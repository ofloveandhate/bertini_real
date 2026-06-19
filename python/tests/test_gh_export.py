"""
Tests for the Grasshopper JSON export (Surface.export_gh_json / Curve.export_gh_json and the
helpers they rely on).

Two tiers:
  * fast pure-logic tests with no decomposition data; and
  * fixture-based invariant tests that run against real decompositions under the repo `test/`
    tree when present (skipped otherwise, so the suite still passes without the data).
"""
import json
import os

import numpy as np
import pytest


# the closed vocabulary of curve-type tags the C# side branches on
VOCAB = {"critical", "sphere", "singular", "midslice", "critslice", "unknown", "standalone"}

_REPO = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", ".."))


def _decomp(rel):
    full = os.path.join(_REPO, rel)
    return full if os.path.isdir(full) else None


WHITNEY = _decomp("test/surface/whitney/output_dim_2_comp_0")     # sampled, has a singular curve
SPHERE = _decomp("test/surface/sphere/output_dim_2_comp_0")       # unsampled
EISTUTE = _decomp("test/curve/intersections_of_surfaces/eistute_sphere/output_dim_1_comp_0")
NORDSTRAND = _decomp("test/surface/nordstrands_weird/output_dim_2_comp_0")  # nodal singularities


# --------------------------------------------------------------------------- #
# fast, pure-logic tests (no decomposition data required)
# --------------------------------------------------------------------------- #

def test_mesh_triangles_flattens():
    from bertini_real.surface import _mesh_triangles

    class FakeMesh:
        faces = np.array([[0, 1, 2], [2, 3, 4]])

    out = _mesh_triangles(FakeMesh())
    assert out["triangles"] == [0, 1, 2, 2, 3, 4]
    assert out["triangle_count"] == 2


def test_mesh_triangles_none():
    from bertini_real.surface import _mesh_triangles
    assert _mesh_triangles(None) is None


def test_points_to_xyz_pads_and_truncates():
    from bertini_real.curve import _points_to_xyz
    out = _points_to_xyz([[1.0], [1.0, 2.0], [1.0, 2.0, 3.0, 4.0]])
    assert out == [[1.0, 0.0, 0.0], [1.0, 2.0, 0.0], [1.0, 2.0, 3.0]]


def test_curve_type_for_name():
    from bertini_real.surface import Surface

    s = Surface.__new__(Surface)  # bypass __init__/file IO

    class C:
        def __init__(self, name):
            self.inputfilename = name

    s.critical_curve = C("input_critical_curve")
    s.sphere_curve = C("input_surf_sphere")
    s.critical_point_slices = [C("crit0")]
    s.midpoint_slices = [C("mid0")]
    s.singular_names = ["sing0"]

    assert s._curve_type_for_name("input_critical_curve") == "critical"
    assert s._curve_type_for_name("input_surf_sphere") == "sphere"
    assert s._curve_type_for_name("crit0") == "critslice"
    assert s._curve_type_for_name("mid0") == "midslice"
    assert s._curve_type_for_name("sing0") == "singular"
    assert s._curve_type_for_name("nope") == "unknown"


# --------------------------------------------------------------------------- #
# fixture-based invariant tests (need real decomposition data + trimesh)
# --------------------------------------------------------------------------- #

trimesh = pytest.importorskip("trimesh")


def _check_surface_invariants(contents, surface):
    assert contents["decomposition_type"] == "surface"
    assert contents["vertex_count"] == len(contents["vertices"])
    assert all(len(v) == 3 for v in contents["vertices"])
    assert len(contents["sphere"]["center"]) == 3
    assert contents["sphere"]["radius"] > 0
    assert len(contents["pieces"]) == len(surface.separate_into_nonsingular_pieces())

    vc = contents["vertex_count"]
    for p in contents["pieces"]:
        assert p["mesh_raw"] is not None
        if not contents["is_sampled"]:
            assert p["mesh_smooth"] is None
        for key in ("mesh_raw", "mesh_smooth"):
            m = p[key]
            if m is None:
                continue
            assert len(m["triangles"]) == 3 * m["triangle_count"]
            if m["triangles"]:
                assert 0 <= min(m["triangles"]) and max(m["triangles"]) < vc
        for cv in p["curves"]:
            assert cv["type"] in VOCAB
            if cv["vertex_indices"]:
                assert max(cv["vertex_indices"]) < vc


@pytest.mark.skipif(not WHITNEY, reason="whitney decomposition not present")
def test_surface_export_invariants_whitney(tmp_path):
    from bertini_real.surface import Surface
    s = Surface(WHITNEY)
    contents = json.load(open(s.export_gh_json(str(tmp_path / "w.json"))))
    _check_surface_invariants(contents, s)
    # whitney has a singular curve -- make sure that type is actually emitted
    types = {cv["type"] for p in contents["pieces"] for cv in p["curves"]}
    assert "singular" in types


@pytest.mark.skipif(not WHITNEY, reason="whitney decomposition not present")
def test_index_polyline_matches_points_whitney():
    """to_point_indices mapped through the unified vertices must equal the proven to_points()."""
    from bertini_real.surface import Surface
    s = Surface(WHITNEY)
    pts = s.extract_points()
    for piece in s.separate_into_nonsingular_pieces():
        for cp in piece.edge_pieces():
            idx = cp.to_point_indices()
            via_idx = np.array([pts[i] for i in idx])[:, :3] if idx else np.empty((0, 3))
            via_pts = cp.to_points()
            assert via_idx.shape == via_pts.shape
            if via_pts.size:
                assert np.allclose(via_idx, via_pts)


@pytest.mark.skipif(not WHITNEY, reason="whitney decomposition not present")
def test_piece_as_closed_mesh_and_spread():
    """Rhino-free pipeline: sampled pieces cap+join to watertight solids, and spread moves them."""
    from bertini_real.surface import Surface, spread_pieces
    s = Surface(WHITNEY)
    pieces = s.separate_into_nonsingular_pieces()
    closed = [p.as_closed_mesh(resolution=4) for p in pieces]
    assert all(m is not None and m.is_watertight for m in closed)

    moved = spread_pieces(closed, factor=0.5)
    assert len(moved) == len(closed)
    # spreading must actually displace at least one piece, and leave the originals untouched
    assert any(not np.allclose(a.bounds.mean(axis=0), b.bounds.mean(axis=0))
               for a, b in zip(closed, moved))


@pytest.mark.skipif(not SPHERE, reason="unsampled sphere decomposition not present")
def test_unsampled_surface_has_null_smooth(tmp_path):
    from bertini_real.surface import Surface
    s = Surface(SPHERE)
    contents = json.load(open(s.export_gh_json(str(tmp_path / "s.json"))))
    assert contents["is_sampled"] is False
    _check_surface_invariants(contents, s)
    for p in contents["pieces"]:
        assert p["mesh_smooth"] is None


@pytest.mark.skipif(not NORDSTRAND, reason="nordstrand decomposition not present")
def test_surface_export_singularities(tmp_path):
    pytest.importorskip("bertini")  # tangent-cone directions need the bertini parser
    from bertini_real.surface import Surface
    s = Surface(NORDSTRAND)
    contents = json.load(open(s.export_gh_json(str(tmp_path / "n.json"))))
    sg = contents["singularities"]
    n_pieces = len(contents["pieces"])

    # one direction per singularity, both xyz
    assert len(sg["locations"]) == len(sg["directions"])
    assert all(len(p) == 3 for p in sg["locations"])
    assert all(len(d) == 3 for d in sg["directions"])

    # one parity row per singularity, one entry per piece; a connector joins exactly two pieces
    assert len(sg["parities"]) == len(sg["locations"])
    for row in sg["parities"]:
        assert len(row) == n_pieces
        assert row.count(-1) == 1 and row.count(1) == 1

    # on_pieces is per piece; its singularity ids are in range
    assert len(sg["on_pieces"]) == n_pieces
    n_sing = len(sg["locations"])
    for ids in sg["on_pieces"]:
        for i in ids:
            assert 0 <= i < n_sing


@pytest.mark.skipif(not EISTUTE, reason="eistute_sphere curve decomposition not present")
def test_curve_export_invariants(tmp_path):
    from bertini_real.curve import Curve
    c = Curve(EISTUTE)
    contents = json.load(open(c.export_gh_json(str(tmp_path / "c.json"))))
    assert contents["decomposition_type"] == "curve"
    assert contents["vertex_count"] == len(contents["vertices"])
    assert len(contents["sphere"]["center"]) == 3
    assert contents["sphere"]["radius"] > 0
    vc = contents["vertex_count"]
    for p in contents["curve_pieces"]:
        assert p["type"] in VOCAB
        if p["vertex_indices"]:
            assert max(p["vertex_indices"]) < vc
