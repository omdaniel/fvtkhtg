#!/usr/bin/env python3
"""Round trip, part 2: read the VTKHDF files that test_vtkhdf_roundtrip_write
produced, using VTK's own vtkHDFReader (the reader ParaView uses), and compare
every cell with what fvtkhtg said it wrote.

For each line of manifest.txt (<file> <expectation> <step or -1>) this checks:
  * grid: dimensions, branch factor, root indexing, coordinates, interface
    names, number of trees and cells, array names and types, active scalars;
  * every cell, reached by following its child path with a VTK geometry
    cursor: level, bounds (exactly), leaf or refined, masked or not, and
    every cell-data value (exactly);
  * the "standard" arrays against values recomputed here from each cell's
    tree and path, and the changing-mesh arrays against values recomputed
    from the level and bounds VTK itself reports for the cell: both
    independently of fvtkhtg's storage and of the expectation file.

Needs the Python `vtk` package, version 9.7 or later; exits with code 77
(skipped) when it is not available. Uses no other third-party package.

usage: check_vtkhdf_roundtrip.py <directory holding manifest.txt>
"""
import math
import sys
from pathlib import Path

try:
    from vtkmodules.vtkCommonCore import vtkVersion
    from vtkmodules.vtkCommonDataModel import (
        vtkHyperTreeGrid,
        vtkHyperTreeGridNonOrientedGeometryCursor,
    )
    from vtkmodules.vtkIOHDF import vtkHDFReader
except ImportError as exc:
    print(f"SKIP: the Python vtk package is not available ({exc})")
    sys.exit(77)

VTK_VERSION = (
    vtkVersion.GetVTKMajorVersion(),
    vtkVersion.GetVTKMinorVersion(),
    vtkVersion.GetVTKBuildVersion(),
)
if VTK_VERSION < (9, 7, 0):
    print("SKIP: reading VTKHDF HyperTreeGrid files needs VTK 9.7 or later, found "
          + ".".join(str(v) for v in VTK_VERSION))
    sys.exit(77)

# How each fvtkhtg array type must arrive in VTK
VTK_TYPE_NAMES = {
    "REAL64": {"double"},
    "REAL32": {"float"},
    "INT64": {"long long", "long"},
    "INT32": {"int"},
    "BIT": {"unsigned char"},   # the writer has no bit type; bits are written as 0/1 bytes
}


def require(condition, message):
    if not condition:
        raise AssertionError(message)


def path_hash(tree, path):
    h = tree + 1
    for child in path:
        h = h * 29 + child + 1
    return h


def standard_values(name, tree, path):
    """Values the standard arrays must hold, derived only from the cell's position."""
    h = path_hash(tree, path)
    level = len(path)
    require(h < 2 ** 24, "hash too large to be exact in REAL32")
    if name == "Pressure":
        return [0.25 * h + 0.5]
    if name == "Velocity":
        return [float(h), -float(h), 0.5 * level]
    if name == "Id":
        return [h]
    if name == "Count":
        return [level, path[-1] if path else -1]
    if name == "Flag":
        return [1 if h % 3 == 0 else 0]
    return None


def moving_values(name, step, level, bounds):
    """Values the changing-mesh arrays must hold at a time step, derived only
    from what VTK reports for the cell (see case_dynamic in the writer).
    Returns the values and the tolerance they are known to."""
    if name == "Level":
        return [level], 0.0
    if name == "Distance":       # from the cell centre to a circle of radius 0.6 that moves in x
        dx = 0.5 * (bounds[0] + bounds[1]) - (0.9 + 1.1 * step)
        dy = 0.5 * (bounds[2] + bounds[3]) - 1.0
        return [math.hypot(dx, dy) - 0.6], 1.0e-12
    return None


def parse_expectation(path):
    info = {"arrays": [], "cells": [], "interface": None, "active": None}
    for line in path.read_text().splitlines():
        tokens = line.split()
        if not tokens:
            continue
        key = tokens[0]
        if key == "case":
            info["case"] = tokens[1]
        elif key == "dims":
            info["dims"] = tuple(int(t) for t in tokens[1:4])
        elif key == "branch_factor":
            info["branch_factor"] = int(tokens[1])
        elif key == "transposed":
            info["transposed"] = tokens[1] == "1"
        elif key == "interface":
            if tokens[1] == "1":
                info["interface"] = (tokens[2], tokens[3])
        elif key in ("xcoords", "ycoords", "zcoords"):
            count = int(tokens[1])
            info[key] = [float(t) for t in tokens[2:2 + count]]
        elif key == "active_scalars":
            info["active"] = None if tokens[1] == "-" else tokens[1]
        elif key == "has_mask":
            info["has_mask"] = tokens[1] == "1"
        elif key == "array":
            info["arrays"].append((tokens[1], tokens[2], int(tokens[3])))
        elif key == "ncells":
            info["ncells"] = int(tokens[1])
        elif key == "cell":
            values = {}
            position = 12
            for name, type_name, ncomp in info["arrays"]:
                raw = tokens[position:position + ncomp]
                convert = float if type_name.startswith("REAL") else int
                values[name] = [convert(t) for t in raw]
                position += ncomp
            require(position == len(tokens), f"{path.name}: malformed cell line")
            info["cells"].append({
                "tree": int(tokens[1]),
                "path": [] if tokens[2] == "-" else [int(p) for p in tokens[2].split(".")],
                "level": int(tokens[3]),
                "refined": tokens[4] == "1",
                "masked": tokens[5] == "1",
                "bounds": [float(t) for t in tokens[6:12]],
                "values": values,
            })
        else:
            raise AssertionError(f"{path.name}: unknown line '{key}'")
    require(len(info["cells"]) == info["ncells"], f"{path.name}: cell count line disagrees")
    return info


def read_grid(vtkhdf_path, step):
    reader = vtkHDFReader()
    reader.SetFileName(str(vtkhdf_path))
    reader.UpdateInformation()
    if step >= 0:
        require(reader.GetNumberOfSteps() > step,
                f"{vtkhdf_path.name}: file has {reader.GetNumberOfSteps()} steps, need step {step}")
        reader.SetStep(step)
    reader.Update()
    output = reader.GetOutputDataObject(0)
    if isinstance(output, vtkHyperTreeGrid):
        pieces = [output]
    else:
        pieces = [output.GetPartitionAsDataObject(i) for i in range(output.GetNumberOfPartitions())]
    require(len(pieces) == 1 and isinstance(pieces[0], vtkHyperTreeGrid),
            f"{vtkhdf_path.name}: expected one HyperTreeGrid piece")
    return reader, pieces[0]


def check(vtkhdf_path, expect_path, step):
    want = parse_expectation(expect_path)
    reader, htg = read_grid(vtkhdf_path, step)
    where = vtkhdf_path.name + (f" step {step}" if step >= 0 else "")

    # --- the grid as a whole ---------------------------------------------------
    require(tuple(htg.GetDimensions()) == want["dims"], f"{where}: dimensions")
    require(htg.GetBranchFactor() == want["branch_factor"], f"{where}: branch factor")
    require(bool(htg.GetTransposedRootIndexing()) == want["transposed"], f"{where}: root indexing")
    for key, array in (("xcoords", htg.GetXCoordinates()), ("ycoords", htg.GetYCoordinates()),
                       ("zcoords", htg.GetZCoordinates())):
        got = [array.GetTuple1(i) for i in range(array.GetNumberOfTuples())]
        require(got == want[key], f"{where}: {key} {got} != {want[key]}")
    if want["interface"] is None:
        require(not htg.GetHasInterface(), f"{where}: unexpected material interface")
    else:
        require(htg.GetHasInterface(), f"{where}: material interface not enabled by the reader")
        require(htg.GetInterfaceNormalsName() == want["interface"][0], f"{where}: normals name")
        require(htg.GetInterfaceInterceptsName() == want["interface"][1], f"{where}: intercepts name")
    # A grid with a mask must arrive with one. The converse does not hold: for
    # a changing mesh fvtkhdf always stores a mask (all clear when none was
    # given), so only the per-cell mask state below is compared in that case.
    if want["has_mask"]:
        require(htg.GetMask() is not None, f"{where}: the mask is missing from the file")
    require(htg.GetNumberOfCells() == want["ncells"],
            f"{where}: {htg.GetNumberOfCells()} cells in the file, {want['ncells']} expected")
    trees = sorted({cell["tree"] for cell in want["cells"]})
    require(htg.GetNumberOfNonEmptyTrees() == len(trees), f"{where}: number of trees")

    cell_data = htg.GetCellData()
    names = sorted(cell_data.GetArrayName(i) for i in range(cell_data.GetNumberOfArrays()))
    require(names == sorted(name for name, _, _ in want["arrays"]), f"{where}: cell arrays {names}")
    arrays = {}
    for name, type_name, ncomp in want["arrays"]:
        array = cell_data.GetArray(name)
        require(array.GetNumberOfComponents() == ncomp, f"{where}: {name} components")
        require(array.GetNumberOfTuples() == want["ncells"], f"{where}: {name} tuples")
        require(array.GetDataTypeAsString() in VTK_TYPE_NAMES[type_name],
                f"{where}: {name} arrived as {array.GetDataTypeAsString()}, written as {type_name}")
        arrays[name] = array
    scalars = cell_data.GetScalars()
    require((scalars.GetName() if scalars is not None else None) == want["active"],
            f"{where}: active scalars")

    # --- every cell -----------------------------------------------------------
    cursor = vtkHyperTreeGridNonOrientedGeometryCursor()
    seen = set()
    recomputed = 0
    for cell in want["cells"]:
        label = f"{where}: tree {cell['tree']} path {cell['path'] or 'root'}"
        htg.InitializeNonOrientedGeometryCursor(cursor, cell["tree"], False)
        require(cursor.HasTree(), f"{label}: tree is missing from the file")
        for child in cell["path"]:
            require(not cursor.IsLeaf(), f"{label}: an ancestor is a leaf in the file")
            cursor.ToChild(child)
        require(cursor.GetLevel() == cell["level"], f"{label}: level")
        bounds = [0.0] * 6
        cursor.GetBounds(bounds)
        require(bounds == cell["bounds"], f"{label}: bounds {bounds} != {cell['bounds']}")
        require(cursor.IsLeaf() == (not cell["refined"]), f"{label}: leaf/refined")
        require(bool(cursor.IsMasked()) == cell["masked"], f"{label}: mask")

        index = cursor.GetGlobalNodeIndex()
        require(index not in seen, f"{label}: cell index {index} reached twice")
        seen.add(index)
        for name, _, _ in want["arrays"]:
            got = list(arrays[name].GetTuple(index))
            require(got == cell["values"][name], f"{label}: {name} = {got}, expected {cell['values'][name]}")
            independent = standard_values(name, cell["tree"], cell["path"])
            if independent is not None:
                require(got == independent,
                        f"{label}: {name} = {got}, but its position gives {independent}")
                recomputed += 1
            elif step >= 0:
                moving = moving_values(name, step, cursor.GetLevel(), bounds)
                if moving is not None:
                    values, tolerance = moving
                    require(len(got) == len(values)
                            and all(abs(g - v) <= tolerance for g, v in zip(got, values)),
                            f"{label}: {name} = {got}, but the cell VTK reports gives {values}")
                    recomputed += 1
    require(len(seen) == want["ncells"], f"{where}: not every cell of the file was reached")

    print(f"ok  {where}: {want['ncells']} cells in {len(trees)} trees, {len(arrays)} arrays, "
          f"{recomputed} values recomputed independently")
    return want["ncells"]


def main():
    require(len(sys.argv) == 2, "usage: check_vtkhdf_roundtrip.py <directory>")
    directory = Path(sys.argv[1])
    manifest = directory / "manifest.txt"
    require(manifest.is_file(), f"{manifest} not found: run test_vtkhdf_roundtrip_write first")
    total = 0
    checks = 0
    for line in manifest.read_text().splitlines():
        vtkhdf_name, expect_name, step = line.split()
        total += check(directory / vtkhdf_name, directory / expect_name, int(step))
        checks += 1
    require(checks > 0, "manifest is empty")
    print(f"ALL {checks} ROUND-TRIP CHECKS PASSED with VTK "
          + ".".join(str(v) for v in VTK_VERSION) + f" ({total} cells compared)")


if __name__ == "__main__":
    main()
