"""
Tests for as_mesh_raw producing a clean, manifold mesh.

Regression coverage for the fix where as_mesh_raw emitted degenerate triangles (from
degenerate curve edges) -- a triangle with a repeated vertex index contributes a self-edge
that gets counted twice, which made the raw mesh non-manifold (edges shared by 4 faces) and
introduced duplicate faces, so raw pieces could not be closed into watertight solids.  Also
covers honoring keep_all_vertices so the raw faces index the global (unified) vertex set.
"""
import os
from collections import Counter

import numpy as np
import pytest

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _decomp(rel):
    full = os.path.join(_REPO, rel)
    return full if os.path.isdir(full) else None


WHITNEY = _decomp("test/surface/whitney/output_dim_2_comp_0")       # sampled (raw still derivable)
NORDSTRAND = _decomp("test/surface/nordstrands_weird/output_dim_2_comp_0")  # raw, nodal singularities

trimesh = pytest.importorskip("trimesh")


def _max_edge_face_count(mesh):
    ec = Counter()
    for tri in mesh.faces:
        a, b, c = int(tri[0]), int(tri[1]), int(tri[2])
        for u, v in ((a, b), (b, c), (c, a)):
            ec[frozenset((u, v))] += 1
    return max(ec.values(), default=0)


def _num_duplicate_faces(mesh):
    F = np.asarray(mesh.faces)
    return len(F) - len(np.unique(np.sort(F, axis=1), axis=0))


def _num_degenerate_faces(mesh):
    F = np.asarray(mesh.faces)
    return sum(1 for t in F if len({int(t[0]), int(t[1]), int(t[2])}) < 3)


def _assert_clean_manifold(mesh, global_point_count):
    assert _num_degenerate_faces(mesh) == 0, "raw mesh has degenerate triangles"
    assert _num_duplicate_faces(mesh) == 0, "raw mesh has duplicate faces"
    # every edge borders 1 (boundary) or 2 (interior) faces -- never more
    assert _max_edge_face_count(mesh) <= 2, "raw mesh is non-manifold (edge shared by >2 faces)"
    # keep_all_vertices=True: faces index the global/unified vertex set
    assert len(mesh.vertices) == global_point_count, "raw mesh vertices are not the global set"


@pytest.mark.skipif(not WHITNEY, reason="whitney decomposition not present")
def test_raw_mesh_manifold_whitney():
    from bertini_real.surface import Surface
    s = Surface(WHITNEY)
    gp = len(s.extract_points())
    for p in s.separate_into_nonsingular_pieces():
        _assert_clean_manifold(p.as_mesh(smooth=False), gp)


@pytest.mark.skipif(not NORDSTRAND, reason="nordstrand decomposition not present")
def test_raw_mesh_manifold_and_closes_nordstrand():
    from bertini_real.surface import Surface
    s = Surface(NORDSTRAND)
    gp = len(s.extract_points())
    pieces = s.separate_into_nonsingular_pieces()
    for p in pieces:
        _assert_clean_manifold(p.as_mesh(smooth=False), gp)
    # the headline: every raw piece now caps + joins into a watertight solid
    for p in pieces:
        closed = p.as_closed_mesh(smooth=False, resolution=4)
        assert closed is not None and closed.is_watertight
