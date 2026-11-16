"""LLDB type formatters, DWARF-driven Fortran 2018 memory inspectors, and batch
lockstep trace command for `fvtkhtg` vs C++ VTK 9.7+.

Can be:
1. Imported from Python 3.14 (`uv run python`) by `tools/lldb/verify_lockstep.py`
   for DWARF DIE layout extraction (`DwarfFortranInspector`).
2. Loaded directly inside Apple `/usr/bin/lldb` via:
       command script import tools/lldb/htg_formatters.py
   which registers `vtkBitArray` and `vtkHyperTree` pretty-printers AND the
   `htg_lockstep_trace` batch command used by `verify_lockstep.py`.
"""

import json
import os
import re
import shlex
import struct
import subprocess

try:
    import lldb
except ImportError:
    lldb = None


# ---------------------------------------------------------------------------
# Low-level process memory helpers (used inside `/usr/bin/lldb`)
# ---------------------------------------------------------------------------

def read_process_bytes(process, addr, nbytes):
    if not process or not addr or nbytes <= 0:
        return b""
    err = lldb.SBError()
    data = process.ReadMemory(int(addr), int(nbytes), err)
    if err.Fail() or data is None:
        return b""
    return bytes(data)


def read_u64(process, addr):
    raw = read_process_bytes(process, addr, 8)
    return struct.unpack("<Q", raw)[0] if len(raw) == 8 else 0


def read_i64(process, addr):
    raw = read_process_bytes(process, addr, 8)
    return struct.unpack("<q", raw)[0] if len(raw) == 8 else 0


def read_i32(process, addr):
    raw = read_process_bytes(process, addr, 4)
    return struct.unpack("<i", raw)[0] if len(raw) == 4 else 0


def read_u32(process, addr):
    raw = read_process_bytes(process, addr, 4)
    return struct.unpack("<I", raw)[0] if len(raw) == 4 else 0


def read_u8(process, addr):
    raw = read_process_bytes(process, addr, 1)
    return raw[0] if len(raw) == 1 else 0


def read_u32_array(process, addr, count):
    if count <= 0 or not addr:
        return []
    raw = read_process_bytes(process, addr, 4 * count)
    if len(raw) < 4 * count:
        return []
    return list(struct.unpack(f"<{count}I", raw))


def read_i32_array(process, addr, count):
    if count <= 0 or not addr:
        return []
    raw = read_process_bytes(process, addr, 4 * count)
    if len(raw) < 4 * count:
        return []
    return list(struct.unpack(f"<{count}i", raw))


def read_i64_array(process, addr, count):
    if count <= 0 or not addr:
        return []
    raw = read_process_bytes(process, addr, 8 * count)
    if len(raw) < 8 * count:
        return []
    return list(struct.unpack(f"<{count}q", raw))


# ---------------------------------------------------------------------------
# DWARF Struct & Subprogram Layout Inspector (Flang + GFortran unified)
# ---------------------------------------------------------------------------

_DIE_HEADER_RE = re.compile(r"^0x[0-9a-fA-F]+:(\s+)DW_TAG_(\w+)")
_AT_NAME_RE = re.compile(r'^\s+DW_AT_name\s+\("([^"]+)"\)')
_AT_LINKAGE_RE = re.compile(r'^\s+DW_AT_linkage_name\s+\("([^"]+)"\)')
_AT_BYTE_SIZE_RE = re.compile(r"^\s+DW_AT_byte_size\s+\((0x[0-9a-fA-F]+|\d+)\)")
_AT_TYPE_RE = re.compile(r'^\s+DW_AT_type\s+\(0x[0-9a-fA-F]+\s+"([^"]+)"\)')
_AT_DATA_LOC_RE = re.compile(r"^\s+DW_AT_data_member_location\s+\((0x[0-9a-fA-F]+|\d+)\)")


class DwarfFortranInspector:
    """Parses `dwarfdump` output from `libfvtkhtg.a` to resolve exact member
    byte offsets and `class(...)` conventions across LLVM `flang` & GNU `gfortran`.
    """

    def __init__(self, binary_paths=None):
        self.structs = {}
        self.subprograms = {}
        if isinstance(binary_paths, str):
            binary_paths = [binary_paths]
        for path in (binary_paths or []):
            if path and os.path.exists(path):
                self._parse_dwarfdump(path)

    @classmethod
    def from_dict(cls, data):
        inst = cls(None)
        inst.structs = data.get("structs", {})
        inst.subprograms = data.get("subprograms", {})
        return inst

    def to_dict(self):
        return {"structs": self.structs, "subprograms": self.subprograms}

    def _parse_dwarfdump(self, path):
        out = subprocess.check_output(
            ["/usr/bin/dwarfdump", path], text=True, errors="replace"
        )
        cur_struct = None
        cur_struct_indent = -1
        cur_member = None
        cur_member_indent = -1

        cur_subprog = None
        cur_subprog_indent = -1
        cur_param = None
        cur_param_indent = -1

        for line in out.splitlines():
            m_die = _DIE_HEADER_RE.match(line)
            if m_die:
                indent = len(m_die.group(1))
                tag = m_die.group(2)

                if cur_member is not None and indent <= cur_member_indent:
                    self._commit_member(cur_struct, cur_member)
                    cur_member = None
                if cur_struct is not None and indent <= cur_struct_indent:
                    cur_struct = None

                if cur_param is not None and indent <= cur_param_indent:
                    self._commit_param(cur_subprog, cur_param)
                    cur_param = None
                if cur_subprog is not None and indent <= cur_subprog_indent:
                    self._commit_subprog(cur_subprog)
                    cur_subprog = None

                if tag == "structure_type":
                    cur_struct = {"name": None, "byte_size": 0, "members": {}}
                    cur_struct_indent = indent
                elif tag == "member" and cur_struct is not None and indent == cur_struct_indent + 2:
                    cur_member = {"name": None, "offset": 0, "type": ""}
                    cur_member_indent = indent
                elif tag == "subprogram":
                    cur_subprog = {"name": None, "linkage": None, "self_is_class": None}
                    cur_subprog_indent = indent
                elif (
                    tag == "formal_parameter"
                    and cur_subprog is not None
                    and indent == cur_subprog_indent + 2
                ):
                    cur_param = {"name": None, "type": ""}
                    cur_param_indent = indent
                continue

            if cur_member is not None:
                m = _AT_NAME_RE.match(line)
                if m:
                    cur_member["name"] = m.group(1).lower()
                    continue
                m = _AT_TYPE_RE.match(line)
                if m:
                    cur_member["type"] = m.group(1)
                    continue
                m = _AT_DATA_LOC_RE.match(line)
                if m:
                    cur_member["offset"] = int(m.group(1), 0)
                    continue
            elif cur_struct is not None:
                m = _AT_NAME_RE.match(line)
                if m:
                    sname = m.group(1).split("::")[-1].lower()
                    cur_struct["name"] = sname
                    if sname not in self.structs:
                        self.structs[sname] = cur_struct
                    else:
                        cur_struct = self.structs[sname]
                    continue
                m = _AT_BYTE_SIZE_RE.match(line)
                if m:
                    cur_struct["byte_size"] = max(cur_struct["byte_size"], int(m.group(1), 0))
                    continue

            if cur_param is not None:
                m = _AT_NAME_RE.match(line)
                if m:
                    cur_param["name"] = m.group(1).lower()
                    continue
                m = _AT_TYPE_RE.match(line)
                if m:
                    cur_param["type"] = m.group(1)
                    continue
            elif cur_subprog is not None:
                m = _AT_NAME_RE.match(line)
                if m:
                    cur_subprog["name"] = m.group(1)
                    continue
                m = _AT_LINKAGE_RE.match(line)
                if m:
                    cur_subprog["linkage"] = m.group(1)
                    continue

        if cur_member is not None:
            self._commit_member(cur_struct, cur_member)
        if cur_param is not None:
            self._commit_param(cur_subprog, cur_param)
        if cur_subprog is not None:
            self._commit_subprog(cur_subprog)

    @staticmethod
    def _commit_member(cur_struct, cur_member):
        if cur_struct and cur_struct.get("name") and cur_member.get("name"):
            cur_struct["members"][cur_member["name"]] = {
                "offset": cur_member["offset"],
                "type": cur_member["type"],
            }

    @staticmethod
    def _commit_param(cur_subprog, cur_param):
        if not cur_subprog or cur_subprog["self_is_class"] is not None:
            return
        pname = cur_param.get("name")
        ptype = cur_param.get("type", "")
        if pname in ("self", "grid", "tree", "cursor"):
            is_class = ("__class_" in ptype) or ptype.endswith("*")
            cur_subprog["self_is_class"] = is_class

    def _commit_subprog(self, cur_subprog):
        is_cls = bool(cur_subprog["self_is_class"])
        for key in (cur_subprog.get("name"), cur_subprog.get("linkage")):
            if key:
                entry = {"name": cur_subprog.get("name") or key, "self_is_class": is_cls}
                self.subprograms[key] = entry
                if key.startswith("_"):
                    self.subprograms[key[1:]] = entry

    def resolve_self_addr_at_entry(self, process, frame):
        reg_x0 = frame.FindRegister("x0")
        x0 = reg_x0.GetValueAsUnsigned(0) if reg_x0 and reg_x0.IsValid() else 0
        if not x0:
            return 0
        fname = frame.GetFunctionName() or ""
        sym = frame.GetSymbol().GetName() if frame.GetSymbol() else ""
        is_class = True
        for candidate in (fname, sym, fname.split("::")[-1]):
            if candidate in self.subprograms:
                is_class = self.subprograms[candidate]["self_is_class"]
                break
        return read_u64(process, x0) if is_class else x0

    def member_info(self, struct_name, member_name):
        s = self.structs.get(struct_name.lower())
        if not s:
            return None
        return s["members"].get(member_name.lower())

    def read_scalar(self, process, base_addr, struct_name, member_name, default=0):
        info = self.member_info(struct_name, member_name)
        if not info or not base_addr:
            return default
        addr = base_addr + info["offset"]
        t = info["type"].lower()
        if "kind=1" in t or "bool" in t or "char" in t:
            return read_u8(process, addr)
        if "kind=4" in t or t == "int" or t == "logical":
            return read_i32(process, addr)
        if "*" in t:
            return read_u64(process, addr)
        return read_i64(process, addr)

    def read_allocatable_data_ptr(self, process, base_addr, struct_name, member_name):
        info = self.member_info(struct_name, member_name)
        if not info or not base_addr:
            return 0
        return read_u64(process, base_addr + info["offset"])

    def read_fixed_int_array(self, process, base_addr, struct_name, member_name, count):
        info = self.member_info(struct_name, member_name)
        if not info or not base_addr or count <= 0:
            return []
        addr = base_addr + info["offset"]
        t = info["type"].lower()
        if "kind=8" in t or "int64" in t or "long" in t:
            return read_i64_array(process, addr, count)
        return read_i32_array(process, addr, count)

    def extract_tree_at(self, process, tree_addr):
        if not tree_addr:
            return {}
        nlevels = self.read_scalar(process, tree_addr, "htg_tree_t", "number_of_levels", 0)
        nverts = self.read_scalar(process, tree_addr, "htg_tree_t", "number_of_vertices", 0)
        nnodes = self.read_scalar(process, tree_addr, "htg_tree_t", "number_of_nodes", 0)
        nleaves = self.read_scalar(process, tree_addr, "htg_tree_t", "number_of_leaves", 0)
        bf = self.read_scalar(process, tree_addr, "htg_tree_t", "branch_factor", 0)
        dim = self.read_scalar(process, tree_addr, "htg_tree_t", "dimension", 0)
        nch = self.read_scalar(process, tree_addr, "htg_tree_t", "number_of_children", 0)
        tid = self.read_scalar(process, tree_addr, "htg_tree_t", "tree_index", -1)
        elder_len = self.read_scalar(process, tree_addr, "htg_tree_t", "elder_active_len", 0)

        ptec_ptr = self.read_allocatable_data_ptr(
            process, tree_addr, "htg_tree_t", "parent_to_elder_child"
        )
        vpd_ptr = self.read_allocatable_data_ptr(
            process, tree_addr, "htg_tree_t", "number_of_vertices_per_depth"
        )
        return {
            "tree_addr": int(tree_addr),
            "tree_index": int(tid),
            "branch_factor": int(bf),
            "dimension": int(dim),
            "number_of_children": int(nch),
            "number_of_levels": int(nlevels),
            "number_of_vertices": int(nverts),
            "number_of_nodes": int(nnodes),
            "number_of_leaves": int(nleaves),
            "elder_active_len": int(elder_len),
            "ptec_addr": int(ptec_ptr),
            "vpd_addr": int(vpd_ptr),
            "parent_to_elder_child": read_i64_array(process, ptec_ptr, elder_len),
            "vertices_per_depth": read_i64_array(process, vpd_ptr, nlevels),
        }

    def extract_bit_array_at(self, process, ba_addr):
        if not ba_addr:
            return {"num_bits": 0, "active_bytes": 0, "base_addr": 0, "bytes_hex": ""}
        nbits = self.read_scalar(process, ba_addr, "htg_bit_array_t", "num_bits", 0)
        bcap = self.read_scalar(process, ba_addr, "htg_bit_array_t", "byte_capacity", 0)
        nbytes = (nbits + 7) // 8 if nbits > 0 else 0
        bptr = self.read_allocatable_data_ptr(process, ba_addr, "htg_bit_array_t", "bytes")
        raw = read_process_bytes(process, bptr, nbytes) if nbytes > 0 else b""
        return {
            "num_bits": int(nbits),
            "byte_capacity": int(bcap),
            "active_bytes": int(nbytes),
            "base_addr": int(bptr),
            "bytes_hex": raw.hex(" "),
        }

    def extract_grid_dims_at(self, process, grid_addr):
        if not grid_addr:
            return {}
        dim = self.read_scalar(process, grid_addr, "htg_grid_t", "dimension", 0)
        orient = self.read_scalar(process, grid_addr, "htg_grid_t", "orientation", -1)
        dims = self.read_fixed_int_array(process, grid_addr, "htg_grid_t", "dimensions", 3)
        cdims = self.read_fixed_int_array(process, grid_addr, "htg_grid_t", "cell_dims", 3)
        axes_name = "axes" if self.member_info("htg_grid_t", "axes") else "axis"
        axes = self.read_fixed_int_array(process, grid_addr, "htg_grid_t", axes_name, 2)
        return {
            "grid_addr": int(grid_addr),
            "dimension": int(dim),
            "orientation": int(orient),
            "dimensions": dims,
            "cell_dims": cdims,
            "axes": axes,
        }

    def extract_grid_bit_member_at(self, process, grid_addr, member_name):
        info = self.member_info("htg_grid_t", member_name)
        if not info or not grid_addr:
            return {"num_bits": 0, "active_bytes": 0, "base_addr": 0, "bytes_hex": ""}
        slot_addr = grid_addr + info["offset"]
        ba_addr = read_u64(process, slot_addr) if ("*" in info["type"]) else slot_addr
        return self.extract_bit_array_at(process, ba_addr)

    def extract_cursor_at(self, process, cursor_addr):
        if not cursor_addr:
            return {}
        last_valid = self.read_scalar(process, cursor_addr, "htg_cursor_t", "last_valid_entry", -1)
        level = max(0, last_valid)
        entries_ptr = self.read_allocatable_data_ptr(
            process, cursor_addr, "htg_cursor_t", "entries"
        )
        stack = []
        if entries_ptr and last_valid >= 0:
            entry_struct = self.structs.get("htg_entry_t")
            if entry_struct and entry_struct["byte_size"] > 0:
                esize = entry_struct["byte_size"]
                vid_field = "vertex_id" if "vertex_id" in entry_struct["members"] else "id"
                for lvl in range(last_valid + 1):
                    eaddr = entries_ptr + lvl * esize
                    if vid_field in entry_struct["members"]:
                        vid = self.read_scalar(process, eaddr, "htg_entry_t", vid_field, -1)
                    else:
                        vid = read_i64(process, eaddr)
                    stack.append(int(vid))
            else:
                stack = read_i64_array(process, entries_ptr, last_valid + 1)
        return {
            "cursor_addr": int(cursor_addr),
            "level": int(level),
            "last_valid_entry": int(last_valid),
            "entries_addr": int(entries_ptr),
            "entry_vertex_ids": stack,
        }


# ---------------------------------------------------------------------------
# C++ VTK `SBValue` Inspectors (inside `/usr/bin/lldb`)
# ---------------------------------------------------------------------------

def extract_vtk_bit_array_state(valobj):
    if not valobj or not valobj.IsValid():
        return {"num_bits": 0, "active_bytes": 0, "base_addr": 0, "bytes_hex": ""}
    if valobj.GetType().IsPointerType() or valobj.GetType().IsReferenceType():
        if valobj.GetValueAsUnsigned(0) == 0:
            return {"num_bits": 0, "active_bytes": 0, "base_addr": 0, "bytes_hex": ""}
        valobj = valobj.Dereference()
    obj_member = valobj.GetChildMemberWithName("Object")
    if obj_member and obj_member.IsValid():
        obj_addr = obj_member.GetValueAsUnsigned(0)
        if obj_addr == 0:
            return {"num_bits": 0, "active_bytes": 0, "base_addr": 0, "bytes_hex": ""}
        target = valobj.GetTarget()
        ba_type = target.FindFirstType("vtkBitArray")
        if ba_type and ba_type.IsValid():
            valobj = target.CreateValueFromAddress(
                "vtk_bit_array_live", lldb.SBAddress(obj_addr, target), ba_type
            )
        else:
            valobj = obj_member.GetDynamicValue(lldb.eDynamicDontRunTarget).Dereference()
    process = valobj.GetProcess()
    max_id_val = valobj.GetChildMemberWithName("MaxId")
    max_id = max_id_val.GetValueAsSigned(-1) if max_id_val and max_id_val.IsValid() else -1
    num_bits = max_id + 1 if max_id >= 0 else 0
    active_bytes = (num_bits + 7) // 8 if num_bits > 0 else 0

    base_addr = 0
    buf_member = valobj.GetChildMemberWithName("Buffer")
    if buf_member and buf_member.IsValid():
        buf_obj = buf_member.Dereference() if buf_member.GetType().IsPointerType() else buf_member
        ptr_member = buf_obj.GetChildMemberWithName("Pointer")
        if ptr_member and ptr_member.IsValid():
            base_addr = ptr_member.GetValueAsUnsigned(0)
    raw = read_process_bytes(process, base_addr, active_bytes) if active_bytes > 0 else b""
    return {
        "num_bits": int(num_bits),
        "active_bytes": int(active_bytes),
        "base_addr": int(base_addr),
        "bytes_hex": raw.hex(" "),
    }


def extract_vtk_tree_state(valobj):
    if not valobj or not valobj.IsValid():
        return {}
    if valobj.GetType().IsPointerType() or valobj.GetType().IsReferenceType():
        if valobj.GetValueAsUnsigned(0) == 0:
            return {}
        valobj = valobj.Dereference()
    process = valobj.GetProcess()

    bf = valobj.GetChildMemberWithName("BranchFactor").GetValueAsUnsigned(0)
    dim = valobj.GetChildMemberWithName("Dimension").GetValueAsUnsigned(0)
    nchildren = valobj.GetChildMemberWithName("NumberOfChildren").GetValueAsUnsigned(0)

    datas = valobj.GetChildMemberWithName("Datas")
    datas_ptr = datas.GetChildMemberWithName("__ptr_")
    if datas_ptr and datas_ptr.IsValid():
        datas_obj = datas_ptr.Dereference()
    else:
        datas_obj = datas.GetChildAtIndex(0)
        if datas_obj.GetType().IsPointerType():
            datas_obj = datas_obj.Dereference()

    tree_idx = datas_obj.GetChildMemberWithName("TreeIndex").GetValueAsSigned(-1)
    nlevels = datas_obj.GetChildMemberWithName("NumberOfLevels").GetValueAsUnsigned(0)
    nverts = datas_obj.GetChildMemberWithName("NumberOfVertices").GetValueAsSigned(0)
    nnodes = datas_obj.GetChildMemberWithName("NumberOfNodes").GetValueAsSigned(0)
    gstart = datas_obj.GetChildMemberWithName("GlobalIndexStart").GetValueAsSigned(-1)

    ptec_vec = datas_obj.GetChildMemberWithName("ParentToElderChild")
    ptec_addr = 0
    ptec_len = 0
    begin_field = ptec_vec.GetChildMemberWithName("__begin_")
    end_field = ptec_vec.GetChildMemberWithName("__end_")
    if begin_field and begin_field.IsValid() and end_field and end_field.IsValid():
        b_addr = begin_field.GetValueAsUnsigned(0)
        e_addr = end_field.GetValueAsUnsigned(0)
        ptec_addr = b_addr
        ptec_len = (e_addr - b_addr) // 4 if e_addr >= b_addr else 0
    else:
        ptec_len = ptec_vec.GetNumChildren()
        if ptec_len > 0:
            ptec_addr = ptec_vec.GetChildAtIndex(0).GetLoadAddress()
    ptec = read_u32_array(process, ptec_addr, ptec_len)

    return {
        "tree_index": int(tree_idx),
        "branch_factor": int(bf),
        "dimension": int(dim),
        "number_of_children": int(nchildren),
        "number_of_levels": int(nlevels),
        "number_of_vertices": int(nverts),
        "number_of_nodes": int(nnodes),
        "number_of_leaves": int(nverts - nnodes),
        "global_index_start": int(gstart),
        "parent_to_elder_child": ptec,
        "ptec_addr": int(ptec_addr),
    }


def extract_vtk_cursor_state(valobj):
    if not valobj or not valobj.IsValid():
        return {}
    if valobj.GetType().IsPointerType() or valobj.GetType().IsReferenceType():
        if valobj.GetValueAsUnsigned(0) == 0:
            return {}
        valobj = valobj.Dereference()
    process = valobj.GetProcess()
    level = valobj.GetChildMemberWithName("Level").GetValueAsUnsigned(0)
    last_valid = valobj.GetChildMemberWithName("LastValidEntry").GetValueAsSigned(-1)
    entries_vec = valobj.GetChildMemberWithName("Entries")
    entries_addr = 0
    begin_field = entries_vec.GetChildMemberWithName("__begin_")
    if begin_field and begin_field.IsValid():
        entries_addr = begin_field.GetValueAsUnsigned(0)
    elif entries_vec.GetNumChildren() > 0:
        entries_addr = entries_vec.GetChildAtIndex(0).GetLoadAddress()
    stack = read_i64_array(process, entries_addr, last_valid + 1) if last_valid >= 0 else []
    return {
        "level": int(level),
        "last_valid_entry": int(last_valid),
        "entries_addr": int(entries_addr),
        "entry_vertex_ids": stack,
    }


def extract_vtk_grid_dims_state(valobj):
    if not valobj or not valobj.IsValid():
        return {}
    if valobj.GetType().IsPointerType() or valobj.GetType().IsReferenceType():
        if valobj.GetValueAsUnsigned(0) == 0:
            return {}
        valobj = valobj.Dereference()
    process = valobj.GetProcess()
    dim = valobj.GetChildMemberWithName("Dimension").GetValueAsUnsigned(0)
    orient = valobj.GetChildMemberWithName("Orientation").GetValueAsUnsigned(0)
    dims_addr = valobj.GetChildMemberWithName("Dimensions").GetLoadAddress()
    cdims_addr = valobj.GetChildMemberWithName("CellDims").GetLoadAddress()
    axis_addr = valobj.GetChildMemberWithName("Axis").GetLoadAddress()
    return {
        "dimension": int(dim),
        "orientation": int(orient),
        "dimensions": read_u32_array(process, dims_addr, 3),
        "cell_dims": read_u32_array(process, cdims_addr, 3),
        "axes": read_u32_array(process, axis_addr, 2),
    }


# ---------------------------------------------------------------------------
# In-LLDB Batch Lockstep Stepping Engine (`htg_lockstep_trace`)
# ---------------------------------------------------------------------------

def _run_lldb_cmd(debugger, cmd):
    interp = debugger.GetCommandInterpreter()
    res = lldb.SBCommandReturnObject()
    interp.HandleCommand(cmd, res)
    return res.Succeeded(), (res.GetOutput() or "").strip(), (res.GetError() or "").strip()


def _create_fortran_bp_no_prologue(debugger, target, regex_pattern):
    ok, out, err = _run_lldb_cmd(debugger, f'breakpoint set -K false -r "{regex_pattern}"')
    if not ok:
        raise RuntimeError(f"Failed to set Fortran breakpoint `{regex_pattern}`: {err}")
    bp = target.GetBreakpointAtIndex(target.GetNumBreakpoints() - 1)
    if not bp or not bp.IsValid():
        raise RuntimeError(f"Invalid SBBreakpoint for `{regex_pattern}` ({out})")
    return bp


def _step_statements_and_finish(process, thread, all_bps, max_steps=10):
    saved_states = [(bp, bp.IsEnabled()) for bp in all_bps]
    for bp, _ in saved_states:
        bp.SetEnabled(False)

    stepped_locs = []
    try:
        start_depth = thread.GetNumFrames()
        caller_fp = thread.GetFrameAtIndex(1).GetFP() if start_depth >= 2 else 0

        def _record_loc(f):
            le = f.GetLineEntry()
            if le and le.IsValid() and le.GetLine() > 0:
                fname = le.GetFileSpec().GetFilename() or "?"
                loc = f"{fname}:{le.GetLine()}"
                if not stepped_locs or stepped_locs[-1] != loc:
                    stepped_locs.append(loc)

        _record_loc(thread.GetFrameAtIndex(0))

        returned = False
        for _ in range(max_steps):
            thread.StepOver()
            if process.GetState() != lldb.eStateStopped:
                returned = True
                break
            cur_depth = thread.GetNumFrames()
            f0 = thread.GetFrameAtIndex(0)
            if cur_depth < start_depth or (caller_fp and f0.GetFP() == caller_fp):
                returned = True
                break
            _record_loc(f0)

        if not returned and process.GetState() == lldb.eStateStopped:
            thread.StepOut()
    finally:
        for bp, was_enabled in saved_states:
            bp.SetEnabled(was_enabled)

    return stepped_locs


def _stopped_bp_ids(thread):
    if thread.GetStopReason() != lldb.eStopReasonBreakpoint:
        return set()
    n = thread.GetStopReasonDataCount()
    return {thread.GetStopReasonDataAtIndex(i) for i in range(0, n, 2)}


def lldb_cmd_lockstep_trace(debugger, command, result, _internal_dict):
    """Execute fine-grained C++/Fortran lockstep stepping inside `/usr/bin/lldb`
    and write structured JSON trace to `--out-json` for Python 3.14 verification.
    """
    try:
        tokens = shlex.split(command)
        driver_bin = None
        dwarf_json = None
        out_json = None
        i = 0
        while i < len(tokens):
            if tokens[i] == "--driver":
                driver_bin = tokens[i + 1]
                i += 2
            elif tokens[i] == "--dwarf-json":
                dwarf_json = tokens[i + 1]
                i += 2
            elif tokens[i] == "--out-json":
                out_json = tokens[i + 1]
                i += 2
            else:
                raise ValueError(f"Unknown argument: {tokens[i]}")

        if not driver_bin or not dwarf_json or not out_json:
            raise ValueError("Required: --driver <bin> --dwarf-json <in.json> --out-json <out.json>")

        with open(dwarf_json, "r", encoding="utf-8") as fh:
            inspector = DwarfFortranInspector.from_dict(json.load(fh))

        debugger.SetAsync(False)
        _run_lldb_cmd(debugger, "settings set target.disable-aslr false")

        ok, _, err = _run_lldb_cmd(debugger, f'target create "{driver_bin}"')
        target = debugger.GetSelectedTarget()
        if not ok or not target or not target.IsValid():
            raise RuntimeError(f"Failed to create LLDB target for {driver_bin}: {err}")

        bp_vtk_set_dims = target.BreakpointCreateByRegex(
            r"^vtkHyperTreeGrid::SetDimensions\(int, int, int\)$"
        )
        bp_vtk_subdivide = target.BreakpointCreateByName("vtkHyperTree::SubdivideLeaf")
        bp_vtk_to_child = target.BreakpointCreateByName("vtkHyperTreeGridNonOrientedCursor::ToChild")
        bp_vtk_build_bf = target.BreakpointCreateByName(
            "vtkHyperTree::BuildFromBreadthFirstOrderDescriptor"
        )
        bp_vtk_pure_mask = target.BreakpointCreateByName("vtkHyperTreeGrid::GetPureMask")

        bp_ftn_set_dims = _create_fortran_bp_no_prologue(
            debugger, target, r"fvtkhtg_grid_m.*grid_set_dimensions_ijk_i32$"
        )
        bp_ftn_subdivide = _create_fortran_bp_no_prologue(
            debugger, target, r"fvtkhtg_tree_m.*tree_subdivide_leaf$"
        )
        bp_ftn_to_child = _create_fortran_bp_no_prologue(
            debugger, target, r"fvtkhtg_cursor_m.*cursor_to_child_i32$"
        )
        bp_ftn_build_bf = _create_fortran_bp_no_prologue(
            debugger, target, r"fvtkhtg_tree_m.*tree_build_from_breadth_first_order_descriptor$"
        )
        bp_ftn_pure_mask = _create_fortran_bp_no_prologue(
            debugger, target, r"fvtkhtg_grid_m.*grid_get_pure_mask$"
        )

        bp_step_sync = target.BreakpointCreateByName("lockstep_on_step")

        named_bps = [
            ("vtk_set_dims", bp_vtk_set_dims),
            ("vtk_subdivide", bp_vtk_subdivide),
            ("vtk_to_child", bp_vtk_to_child),
            ("vtk_build_bf", bp_vtk_build_bf),
            ("vtk_pure_mask", bp_vtk_pure_mask),
            ("ftn_set_dims", bp_ftn_set_dims),
            ("ftn_subdivide", bp_ftn_subdivide),
            ("ftn_to_child", bp_ftn_to_child),
            ("ftn_build_bf", bp_ftn_build_bf),
            ("ftn_pure_mask", bp_ftn_pure_mask),
            ("step_sync", bp_step_sync),
        ]
        all_bps = [bp for _, bp in named_bps]

        launch_info = target.GetLaunchInfo()
        launch_info.SetWorkingDirectory(os.path.dirname(os.path.abspath(driver_bin)))
        launch_info.SetLaunchFlags(launch_info.GetLaunchFlags() & ~lldb.eLaunchFlagDisableASLR)
        target.SetLaunchInfo(launch_info)
        err_obj = lldb.SBError()
        process = target.Launch(launch_info, err_obj)
        if err_obj.Fail() or not process or not process.IsValid():
            raise RuntimeError(f"LLDB process launch failed: {err_obj.GetCString()}")

        bp_locations = {name: int(bp.GetNumLocations()) for name, bp in named_bps}
        for name, nloc in bp_locations.items():
            if nloc < 1:
                process.Kill()
                raise RuntimeError(f"Breakpoint `{name}` resolved to 0 locations after launch!")

        events = []

        while process.GetState() == lldb.eStateStopped:
            thread = process.GetSelectedThread()
            hit_ids = _stopped_bp_ids(thread)
            frame = thread.GetFrameAtIndex(0)

            if not hit_ids:
                _, bt, _ = _run_lldb_cmd(debugger, "bt")
                reason = thread.GetStopDescription(256)
                process.Kill()
                raise RuntimeError(f"Unexpected non-breakpoint stop ({reason}):\n{bt}")

            if bp_vtk_set_dims.GetID() in hit_ids:
                this_val = frame.FindVariable("this")
                this_addr = this_val.GetValueAsUnsigned(0)
                this_type = this_val.GetType().GetPointeeType()
                _, fv_out, _ = _run_lldb_cmd(debugger, "frame variable")
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=10)
                vgrid_obj = target.CreateValueFromAddress(
                    "vgrid_live", lldb.SBAddress(this_addr, target), this_type
                )
                state = extract_vtk_grid_dims_state(vgrid_obj)
                events.append(
                    {
                        "kind": "vtk_set_dims",
                        "stepped_locs": locs,
                        "frame_vars": fv_out,
                        "state": state,
                    }
                )

            elif bp_ftn_set_dims.GetID() in hit_ids:
                grid_addr = inspector.resolve_self_addr_at_entry(process, frame)
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=12)
                state = inspector.extract_grid_dims_at(process, grid_addr)
                _, mem_out, _ = _run_lldb_cmd(
                    debugger, f"memory read --size 4 --format d --count 3 0x{grid_addr:x}"
                )
                events.append(
                    {
                        "kind": "ftn_set_dims",
                        "stepped_locs": locs,
                        "memory_read": mem_out,
                        "state": state,
                    }
                )

            elif bp_vtk_subdivide.GetID() in hit_ids:
                this_val = frame.FindVariable("this")
                this_addr = this_val.GetValueAsUnsigned(0)
                this_type = this_val.GetType().GetPointeeType()
                idx_val = frame.FindVariable("index").GetValueAsSigned(-1)
                lvl_var = frame.FindVariable("depth")
                if not lvl_var or not lvl_var.IsValid():
                    lvl_var = frame.FindVariable("level")
                lvl_num = lvl_var.GetValueAsUnsigned(0) if lvl_var and lvl_var.IsValid() else 0
                before_st = extract_vtk_tree_state(this_val)
                _, fv_out, _ = _run_lldb_cmd(debugger, "frame variable")
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=12)
                vtree_live = target.CreateValueFromAddress(
                    "vtree_live", lldb.SBAddress(this_addr, target), this_type
                )
                after_st = extract_vtk_tree_state(vtree_live)
                _, mem_out, _ = _run_lldb_cmd(
                    debugger,
                    f"memory read --size 4 --format u --count {len(after_st['parent_to_elder_child'])} "
                    f"0x{after_st['ptec_addr']:x}",
                )
                events.append(
                    {
                        "kind": "vtk_subdivide",
                        "vertex_id": int(idx_val),
                        "level": int(lvl_num),
                        "stepped_locs": locs,
                        "frame_vars": fv_out,
                        "memory_read": mem_out,
                        "before": before_st,
                        "after": after_st,
                    }
                )

            elif bp_ftn_subdivide.GetID() in hit_ids:
                tree_addr = inspector.resolve_self_addr_at_entry(process, frame)
                before_st = inspector.extract_tree_at(process, tree_addr)
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=14)
                after_st = inspector.extract_tree_at(process, tree_addr)
                _, mem_ptec, _ = _run_lldb_cmd(
                    debugger,
                    f"memory read --size 8 --format d --count {len(after_st['parent_to_elder_child'])} "
                    f"0x{after_st['ptec_addr']:x}",
                )
                _, mem_vpd, _ = _run_lldb_cmd(
                    debugger,
                    f"memory read --size 8 --format d --count {len(after_st['vertices_per_depth'])} "
                    f"0x{after_st['vpd_addr']:x}",
                )
                events.append(
                    {
                        "kind": "ftn_subdivide",
                        "stepped_locs": locs,
                        "memory_read_ptec": mem_ptec,
                        "memory_read_vpd": mem_vpd,
                        "before": before_st,
                        "after": after_st,
                    }
                )

            elif bp_vtk_to_child.GetID() in hit_ids:
                this_val = frame.FindVariable("this")
                this_addr = this_val.GetValueAsUnsigned(0)
                this_type = this_val.GetType().GetPointeeType()
                ichild = frame.FindVariable("ichild").GetValueAsUnsigned(0)
                _, fv_out, _ = _run_lldb_cmd(debugger, "frame variable")
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=10)
                vcur_live = target.CreateValueFromAddress(
                    "vcur_live", lldb.SBAddress(this_addr, target), this_type
                )
                cur_st = extract_vtk_cursor_state(vcur_live)
                events.append(
                    {
                        "kind": "vtk_to_child",
                        "ichild": int(ichild),
                        "stepped_locs": locs,
                        "frame_vars": fv_out,
                        "state": cur_st,
                    }
                )

            elif bp_ftn_to_child.GetID() in hit_ids:
                cur_addr = inspector.resolve_self_addr_at_entry(process, frame)
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=10)
                cur_st = inspector.extract_cursor_at(process, cur_addr)
                mem_entries = ""
                if cur_st.get("entries_addr"):
                    _, mem_entries, _ = _run_lldb_cmd(
                        debugger,
                        f"memory read --size 8 --format d --count {cur_st['level'] + 1} "
                        f"0x{cur_st['entries_addr']:x}",
                    )
                events.append(
                    {
                        "kind": "ftn_to_child",
                        "stepped_locs": locs,
                        "memory_read": mem_entries,
                        "state": cur_st,
                    }
                )

            elif bp_vtk_build_bf.GetID() in hit_ids:
                this_val = frame.FindVariable("this")
                this_addr = this_val.GetValueAsUnsigned(0)
                this_type = this_val.GetType().GetPointeeType()
                desc_val = frame.FindVariable("descriptor")
                nbits = frame.FindVariable("numberOfBits").GetValueAsSigned(-1)
                sidx = frame.FindVariable("startIndex").GetValueAsSigned(-1)
                desc_st = extract_vtk_bit_array_state(desc_val)
                _, fv_out, _ = _run_lldb_cmd(debugger, "frame variable")
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=16)
                vtree_live = target.CreateValueFromAddress(
                    "vtree_bf", lldb.SBAddress(this_addr, target), this_type
                )
                after_st = extract_vtk_tree_state(vtree_live)
                events.append(
                    {
                        "kind": "vtk_build_bf",
                        "nbits": int(nbits),
                        "start_index": int(sidx),
                        "stepped_locs": locs,
                        "frame_vars": fv_out,
                        "descriptor": desc_st,
                        "after": after_st,
                    }
                )

            elif bp_ftn_build_bf.GetID() in hit_ids:
                tree_addr = inspector.resolve_self_addr_at_entry(process, frame)
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=16)
                after_st = inspector.extract_tree_at(process, tree_addr)
                events.append(
                    {
                        "kind": "ftn_build_bf",
                        "stepped_locs": locs,
                        "after": after_st,
                    }
                )

            elif bp_vtk_pure_mask.GetID() in hit_ids:
                this_val = frame.FindVariable("this")
                this_addr = this_val.GetValueAsUnsigned(0)
                this_type = this_val.GetType().GetPointeeType()
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=12)
                vgrid_live = target.CreateValueFromAddress(
                    "vgrid_pm", lldb.SBAddress(this_addr, target), this_type
                )
                vmask_st = extract_vtk_bit_array_state(vgrid_live.GetChildMemberWithName("Mask"))
                vpure_st = extract_vtk_bit_array_state(vgrid_live.GetChildMemberWithName("PureMask"))
                _, mem_pm, _ = _run_lldb_cmd(
                    debugger,
                    f"memory read --size 1 --format x --count {vpure_st['active_bytes']} "
                    f"0x{vpure_st['base_addr']:x}",
                )
                events.append(
                    {
                        "kind": "vtk_pure_mask",
                        "stepped_locs": locs,
                        "memory_read": mem_pm,
                        "mask": vmask_st,
                        "pure_mask": vpure_st,
                    }
                )

            elif bp_ftn_pure_mask.GetID() in hit_ids:
                grid_addr = inspector.resolve_self_addr_at_entry(process, frame)
                locs = _step_statements_and_finish(process, thread, all_bps, max_steps=14)
                fmask_st = inspector.extract_grid_bit_member_at(process, grid_addr, "mask")
                fpure_st = inspector.extract_grid_bit_member_at(process, grid_addr, "pure_mask")
                mem_pm = ""
                if fpure_st["base_addr"] and fpure_st["active_bytes"] > 0:
                    _, mem_pm, _ = _run_lldb_cmd(
                        debugger,
                        f"memory read --size 1 --format x --count {fpure_st['active_bytes']} "
                        f"0x{fpure_st['base_addr']:x}",
                    )
                events.append(
                    {
                        "kind": "ftn_pure_mask",
                        "stepped_locs": locs,
                        "memory_read": mem_pm,
                        "mask": fmask_st,
                        "pure_mask": fpure_st,
                    }
                )

            elif bp_step_sync.GetID() in hit_ids:
                _, fv_snap, _ = _run_lldb_cmd(debugger, "frame variable *snap")
                snap = frame.FindVariable("snap").Dereference()
                step_id = snap.GetChildMemberWithName("step_id").GetValueAsSigned(-1)
                phase = snap.GetChildMemberWithName("phase_name").GetSummary() or ""
                v_nb = snap.GetChildMemberWithName("vtk_bits_nbytes").GetValueAsSigned(0)
                f_nb = snap.GetChildMemberWithName("ftn_bits_nbytes").GetValueAsSigned(0)
                v_hex = ""
                f_hex = ""
                if v_nb > 0:
                    vp = snap.GetChildMemberWithName("vtk_bits_ptr").GetValueAsUnsigned(0)
                    v_hex = read_process_bytes(process, vp, v_nb).hex(" ")
                if f_nb > 0:
                    fp = snap.GetChildMemberWithName("ftn_bits_ptr").GetValueAsUnsigned(0)
                    f_hex = read_process_bytes(process, fp, f_nb).hex(" ")
                events.append(
                    {
                        "kind": "step_sync",
                        "step_id": int(step_id),
                        "phase": phase,
                        "frame_vars": fv_snap,
                        "vtk_bits_nbytes": int(v_nb),
                        "ftn_bits_nbytes": int(f_nb),
                        "vtk_bits_hex": v_hex,
                        "ftn_bits_hex": f_hex,
                    }
                )
                if step_id == 5:
                    for bp in (
                        bp_vtk_set_dims,
                        bp_ftn_set_dims,
                        bp_vtk_subdivide,
                        bp_ftn_subdivide,
                        bp_vtk_to_child,
                        bp_ftn_to_child,
                        bp_vtk_build_bf,
                        bp_ftn_build_bf,
                    ):
                        bp.SetEnabled(False)
                elif step_id == 6:
                    for bp in all_bps:
                        bp.SetEnabled(False)

            process.Continue()

        payload = {
            "exit_state": int(process.GetState()),
            "eStateExited": int(lldb.eStateExited),
            "exit_code": int(process.GetExitStatus()),
            "stdout": process.GetSTDOUT(65536) or "",
            "stderr": process.GetSTDERR(65536) or "",
            "bp_locations": bp_locations,
            "events": events,
        }
        with open(out_json, "w", encoding="utf-8") as fh:
            json.dump(payload, fh, indent=2)

        result.AppendMessage(
            f"[htg_lockstep_trace] Captured {len(events)} lockstep stop events -> {out_json}"
        )
        result.SetStatus(lldb.eReturnStatusSuccessFinishResult)
    except Exception as exc:  # noqa: BLE001
        import traceback

        result.SetError(f"[htg_lockstep_trace ERROR] {exc}\n{traceback.format_exc()}")


# ---------------------------------------------------------------------------
# LLDB Pretty-Printer Summary Callbacks
# ---------------------------------------------------------------------------

def vtk_bit_array_summary(valobj, _internal_dict):
    st = extract_vtk_bit_array_state(valobj)
    return f"num_bits={st['num_bits']} active_bytes={st['active_bytes']} bytes=[{st['bytes_hex'][:47]}]"


def vtk_tree_summary(valobj, _internal_dict):
    st = extract_vtk_tree_state(valobj)
    if not st:
        return "<null vtkHyperTree>"
    return (
        f"tree_id={st['tree_index']} f={st['branch_factor']} d={st['dimension']} "
        f"levels={st['number_of_levels']} verts={st['number_of_vertices']} "
        f"nodes={st['number_of_nodes']} leaves={st['number_of_leaves']} "
        f"ptec={st['parent_to_elder_child'][:12]}"
    )


def __lldb_init_module(debugger, _internal_dict):
    debugger.HandleCommand(
        'type summary add -F htg_formatters.vtk_bit_array_summary '
        '-x "^(VTK_ABI_NAMESPACE_BEGIN::)?vtkBitArray( \\*)?$"'
    )
    debugger.HandleCommand(
        'type summary add -F htg_formatters.vtk_tree_summary '
        '-x "^(VTK_ABI_NAMESPACE_BEGIN::)?vtkHyperTree( \\*)?$"'
    )
    debugger.HandleCommand(
        "command script add -f htg_formatters.lldb_cmd_lockstep_trace htg_lockstep_trace"
    )
