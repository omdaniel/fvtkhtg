#!/usr/bin/env -S uv run python
"""Automated fine-grained `/usr/bin/lldb` lockstep verifier for C++ VTK 9.7+ vs Fortran `fvtkhtg`.

Managed via `uv` and strictly targeting Python >=3.14,<3.15:
    uv run python tools/lldb/verify_lockstep.py <test_lockstep_driver> [--fortran-archive <libfvtkhtg.a>]

Orchestrates `/usr/bin/lldb --no-lldbinit -b -s <batch.lldb>` (keeping Xcode LLDB's embedded
Python framework isolated inside the `/usr/bin/lldb` process) while parsing Fortran 2018 DWARF
DIE layouts (`dwarfdump`) and asserting 100% step-by-step C++ VTK vs Fortran `fvtkhtg` internal
state equivalence in Python 3.14!
"""

import argparse
import json
import os
import subprocess
import sys
import tempfile

# Enforce repository Python >=3.14,<3.15 requirement (`uv` environment)
if sys.version_info[:2] != (3, 14):
    raise RuntimeError(
        f"verify_lockstep.py strictly requires Python >=3.14,<3.15 via `uv run python`, "
        f"got Python {sys.version.split()[0]} ({sys.executable})"
    )

_SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
if _SCRIPT_DIR not in sys.path:
    sys.path.insert(0, _SCRIPT_DIR)

import htg_formatters as fmt  # noqa: E402


UINT32_MAX = 0xFFFFFFFF


def _find_default_fortran_archive(driver_bin: str) -> str | None:
    cur = os.path.dirname(os.path.abspath(driver_bin))
    for _ in range(4):
        candidate = os.path.join(cur, "libfvtkhtg.a")
        if os.path.isfile(candidate):
            return candidate
        parent = os.path.dirname(cur)
        if parent == cur:
            break
        cur = parent
    return None


def _normalize_ptec(seq: list[int]) -> list[int]:
    return [UINT32_MAX if v < 0 else int(v) for v in seq]


def _assert_eq(lhs, rhs, context: str) -> None:
    if lhs != rhs:
        raise AssertionError(
            f"[LLDB LOCKSTEP MISMATCH] {context}\n  VTK (expected): {lhs}\n  FTN (actual):   {rhs}"
        )


def _assert_tree_lockstep(vtk_state: dict, ftn_state: dict, label: str) -> None:
    for key in (
        "tree_index",
        "branch_factor",
        "dimension",
        "number_of_children",
        "number_of_levels",
        "number_of_vertices",
        "number_of_nodes",
        "number_of_leaves",
    ):
        _assert_eq(vtk_state[key], ftn_state[key], f"{label} -> {key}")

    vtk_ptec = _normalize_ptec(vtk_state["parent_to_elder_child"])
    ftn_ptec = _normalize_ptec(ftn_state["parent_to_elder_child"])
    _assert_eq(len(vtk_ptec), len(ftn_ptec), f"{label} -> ParentToElderChild.size()")
    _assert_eq(vtk_ptec, ftn_ptec, f"{label} -> ParentToElderChild[]")

    vpd = ftn_state["vertices_per_depth"]
    _assert_eq(len(vpd), ftn_state["number_of_levels"], f"{label} -> len(VPD)")
    _assert_eq(sum(vpd), ftn_state["number_of_vertices"], f"{label} -> sum(VPD) == NumberOfVertices")


def run_verification(driver_bin: str, fortran_archive: str | None) -> None:
    print(f"[Python] Verified uv Python runtime: {sys.version.split()[0]} ({sys.executable})")

    driver_bin = os.path.abspath(driver_bin)
    if not os.path.isfile(driver_bin):
        raise FileNotFoundError(f"Lockstep driver binary not found: {driver_bin}")

    archive_paths = [
        os.path.abspath(p) for p in (fortran_archive, driver_bin) if p and os.path.exists(p)
    ]
    bridge_cand = os.path.join(os.path.dirname(driver_bin), "libfvtkhtg_test_c_bridge.a")
    if os.path.isfile(bridge_cand):
        archive_paths.append(bridge_cand)

    print(f"[LLDB] Parsing Fortran 2018 DWARF DIE layouts via dwarfdump: {archive_paths}")
    inspector = fmt.DwarfFortranInspector(archive_paths)
    for req_struct in ("htg_bit_array_t", "htg_tree_t", "htg_grid_t", "htg_cursor_t"):
        if req_struct not in inspector.structs:
            raise AssertionError(
                f"Missing Fortran DWARF struct `{req_struct}` in {archive_paths}; "
                f"found: {sorted(inspector.structs.keys())}"
            )
    print(
        "[LLDB] Verified DWARF structs: "
        + ", ".join(
            f"{k} ({inspector.structs[k]['byte_size']}B, {len(inspector.structs[k]['members'])} fields)"
            for k in ("htg_bit_array_t", "htg_tree_t", "htg_grid_t", "htg_cursor_t")
        )
    )

    formatters_py = os.path.join(_SCRIPT_DIR, "htg_formatters.py")

    with tempfile.TemporaryDirectory(prefix="fvtkhtg_lldb_") as tmpdir:
        dwarf_json = os.path.join(tmpdir, "dwarf_layout.json")
        trace_json = os.path.join(tmpdir, "lockstep_trace.json")
        batch_lldb = os.path.join(tmpdir, "run_lockstep.lldb")

        with open(dwarf_json, "w", encoding="utf-8") as fh:
            json.dump(inspector.to_dict(), fh)

        with open(batch_lldb, "w", encoding="utf-8") as fh:
            fh.write("settings set target.disable-aslr false\n")
            fh.write(f'command script import "{formatters_py}"\n')
            fh.write(
                f'htg_lockstep_trace --driver "{driver_bin}" '
                f'--dwarf-json "{dwarf_json}" --out-json "{trace_json}"\n'
            )

        proc = subprocess.run(
            ["/usr/bin/lldb", "--no-lldbinit", "-b", "-s", batch_lldb],
            text=True,
            capture_output=True,
            check=False,
        )
        if proc.returncode != 0 or not os.path.isfile(trace_json):
            raise RuntimeError(
                f"`/usr/bin/lldb -b` failed (rc={proc.returncode}):\n"
                f"--- STDOUT ---\n{proc.stdout}\n--- STDERR ---\n{proc.stderr}"
            )

        with open(trace_json, "r", encoding="utf-8") as fh:
            trace = json.load(fh)

    if trace["stdout"].strip():
        print("[Driver STDOUT]\n" + trace["stdout"].strip())
    if trace["stderr"].strip():
        print("[Driver STDERR]\n" + trace["stderr"].strip(), file=sys.stderr)

    _assert_eq(trace["exit_state"], trace["eStateExited"], "Process final LLDB state")
    _assert_eq(trace["exit_code"], 0, "test_lockstep_driver exit status")

    pending_set_dims: list[dict] = []
    pending_subdivide: list[dict] = []
    pending_to_child: list[dict] = []
    pending_build_bf: list[dict] = []
    pending_pure_mask: list[dict] = []

    counts = {
        "vtk_set_dims": 0,
        "ftn_set_dims": 0,
        "vtk_subdivide": 0,
        "ftn_subdivide": 0,
        "vtk_to_child": 0,
        "ftn_to_child": 0,
        "vtk_build_bf": 0,
        "ftn_build_bf": 0,
        "vtk_pure_mask": 0,
        "ftn_pure_mask": 0,
        "step_sync": 0,
    }

    for ev in trace["events"]:
        kind = ev["kind"]
        counts[kind] += 1

        if kind == "vtk_set_dims":
            pending_set_dims.append(ev)
        elif kind == "ftn_set_dims":
            _assert_eq(bool(pending_set_dims), True, "Paired VTK SetDimensions queued")
            vev = pending_set_dims.pop(0)
            vst, fst = vev["state"], ev["state"]
            for k in ("dimension", "orientation", "dimensions", "cell_dims", "axes"):
                _assert_eq(vst[k], fst[k], f"SetDimensions lockstep `{k}`")
            print(
                f"[LLDB][SetDimensions] VTK({vev['stepped_locs']}) == FTN({ev['stepped_locs']}) -> "
                f"dims={fst['dimensions']} cell_dims={fst['cell_dims']} orient={fst['orientation']} axes={fst['axes']}"
            )

        elif kind == "vtk_subdivide":
            pending_subdivide.append(ev)
        elif kind == "ftn_subdivide":
            _assert_eq(bool(pending_subdivide), True, "Paired VTK SubdivideLeaf queued")
            vev = pending_subdivide.pop(0)
            vid = vev["vertex_id"]
            _assert_tree_lockstep(vev["before"], ev["before"], f"SubdivideLeaf(v={vid}) BEFORE")
            _assert_tree_lockstep(vev["after"], ev["after"], f"SubdivideLeaf(v={vid}) AFTER")
            print(
                f"[LLDB][SubdivideLeaf(v={vid}, lvl={vev['level']})] "
                f"VTK({vev['stepped_locs']}) == FTN({ev['stepped_locs']}) -> "
                f"ptec={ev['after']['parent_to_elder_child']} vpd={ev['after']['vertices_per_depth']}"
            )

        elif kind == "vtk_to_child":
            pending_to_child.append(ev)
        elif kind == "ftn_to_child":
            _assert_eq(bool(pending_to_child), True, "Paired VTK ToChild queued")
            vev = pending_to_child.pop(0)
            ichild = vev["ichild"]
            vst, fst = vev["state"], ev["state"]
            _assert_eq(vst["level"], fst["level"], f"ToChild({ichild}) Level")
            _assert_eq(
                vst["entry_vertex_ids"][: vst["level"] + 1],
                fst["entry_vertex_ids"][: fst["level"] + 1],
                f"ToChild({ichild}) Entry stack vertex IDs",
            )
            print(
                f"[LLDB][ToChild({ichild})] VTK({vev['stepped_locs']}) == FTN({ev['stepped_locs']}) -> "
                f"level={fst['level']} entries={fst['entry_vertex_ids']}"
            )

        elif kind == "vtk_build_bf":
            pending_build_bf.append(ev)
        elif kind == "ftn_build_bf":
            _assert_eq(bool(pending_build_bf), True, "Paired VTK BuildFromBF queued")
            vev = pending_build_bf.pop(0)
            _assert_tree_lockstep(
                vev["after"],
                ev["after"],
                f"BuildFromBF(nbits={vev['nbits']}, start={vev['start_index']}) AFTER",
            )
            print(
                f"[LLDB][BuildFromBF(nbits={vev['nbits']}, start={vev['start_index']})] "
                f"VTK({vev['stepped_locs']}) == FTN({ev['stepped_locs']}) -> "
                f"desc_bytes=[{vev['descriptor']['bytes_hex']}] "
                f"ptec={ev['after']['parent_to_elder_child']} vpd={ev['after']['vertices_per_depth']}"
            )

        elif kind == "vtk_pure_mask":
            pending_pure_mask.append(ev)
        elif kind == "ftn_pure_mask":
            _assert_eq(bool(pending_pure_mask), True, "Paired VTK GetPureMask queued")
            vev = pending_pure_mask.pop(0)
            _assert_eq(vev["mask"]["num_bits"], ev["mask"]["num_bits"], "Grid Mask num_bits")
            _assert_eq(
                vev["mask"]["bytes_hex"], ev["mask"]["bytes_hex"], "Grid Mask MSB-first bytes"
            )
            _assert_eq(
                vev["pure_mask"]["num_bits"], ev["pure_mask"]["num_bits"], "Grid PureMask num_bits"
            )
            _assert_eq(
                vev["pure_mask"]["bytes_hex"],
                ev["pure_mask"]["bytes_hex"],
                "Grid PureMask MSB-first bytes",
            )
            print(
                f"[LLDB][GetPureMask] VTK({vev['stepped_locs']}) == FTN({ev['stepped_locs']}) -> "
                f"nbits={ev['pure_mask']['num_bits']} mask=[{ev['mask']['bytes_hex']}] "
                f"pure_mask=[{ev['pure_mask']['bytes_hex']}]"
            )

        elif kind == "step_sync":
            if ev["vtk_bits_nbytes"] > 0 or ev["ftn_bits_nbytes"] > 0:
                _assert_eq(
                    ev["vtk_bits_nbytes"],
                    ev["ftn_bits_nbytes"],
                    f"Snapshot #{ev['step_id']} bits_nbytes",
                )
                _assert_eq(
                    ev["vtk_bits_hex"],
                    ev["ftn_bits_hex"],
                    f"Snapshot #{ev['step_id']} raw bit bytes",
                )
            print(f"[LLDB][Sync Step #{ev['step_id']} {ev['phase']}] verified *snap via `frame variable`")

    _assert_eq(len(pending_set_dims), 0, "Unpaired SetDimensions stops")
    _assert_eq(len(pending_subdivide), 0, "Unpaired SubdivideLeaf stops")
    _assert_eq(len(pending_to_child), 0, "Unpaired ToChild stops")
    _assert_eq(len(pending_build_bf), 0, "Unpaired BuildFromBF stops")
    _assert_eq(len(pending_pure_mask), 0, "Unpaired GetPureMask stops")

    expected_min = {
        "vtk_set_dims": 1,
        "ftn_set_dims": 1,
        "vtk_subdivide": 3,
        "ftn_subdivide": 3,
        "vtk_to_child": 2,
        "ftn_to_child": 2,
        "vtk_build_bf": 1,
        "ftn_build_bf": 1,
        "vtk_pure_mask": 1,
        "ftn_pure_mask": 1,
        "step_sync": 6,
    }
    for key, min_hits in expected_min.items():
        if counts[key] < min_hits:
            raise AssertionError(
                f"Expected >= {min_hits} LLDB lockstep hits for `{key}`, got {counts[key]} ({counts})"
            )

    print(f"[PASS] LLDB fine-grained lockstep verification passed 100%! Hit summary: {counts}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Verify C++ VTK 9.7+ vs Fortran fvtkhtg lockstep under /usr/bin/lldb."
    )
    parser.add_argument("driver_binary", help="Path to compiled `test_lockstep_driver` executable.")
    parser.add_argument(
        "--fortran-archive",
        default=None,
        help="Optional path to `libfvtkhtg.a` for DWARF DIE layout inspection.",
    )
    args = parser.parse_args()
    archive = args.fortran_archive or _find_default_fortran_archive(args.driver_binary)
    run_verification(args.driver_binary, archive)


if __name__ == "__main__":
    main()
