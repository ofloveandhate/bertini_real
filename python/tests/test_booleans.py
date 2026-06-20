"""
Tests for mesh_boolean_fold (the ordered boolean fold; Python twin of the Grasshopper
Boolean Piece component).  Uses synthetic primitives, so no decomposition data is needed.
"""
import pytest

trimesh = pytest.importorskip("trimesh")
pytest.importorskip("manifold3d")  # the exact boolean backend


def test_default_is_subtract():
    from bertini_real.surface import mesh_boolean_fold
    box = trimesh.creation.box(extents=(2, 2, 2))           # volume 8
    drill = trimesh.creation.cylinder(radius=0.3, height=4)
    out = mesh_boolean_fold(box, [drill])                   # no signs -> subtract
    assert out.is_watertight
    assert out.volume < box.volume


def test_union_then_subtract_order():
    from bertini_real.surface import mesh_boolean_fold
    box = trimesh.creation.box(extents=(2, 2, 2))
    add = trimesh.creation.box(extents=(2, 2, 2))
    add.apply_translation((1.5, 0, 0))                      # overlaps box -> grows it
    drill = trimesh.creation.cylinder(radius=0.3, height=10)

    out = mesh_boolean_fold(box, [add, drill], signs=[+1, -1])
    assert out.is_watertight
    # the drill pierces the unioned body, so the hole is present in the combined solid
    assert out.volume < trimesh.boolean.union([box, add], engine='manifold').volume


def test_order_matters():
    """subtract-then-union differs from union-then-subtract when the tools overlap."""
    from bertini_real.surface import mesh_boolean_fold
    box = trimesh.creation.box(extents=(2, 2, 2))
    add = trimesh.creation.box(extents=(2, 2, 2))
    add.apply_translation((1.0, 0, 0))
    drill = trimesh.creation.cylinder(radius=0.4, height=10)
    drill.apply_translation((1.0, 0, 0))                    # sits where `add` will be

    union_first = mesh_boolean_fold(box, [add, drill], signs=[+1, -1])
    subtract_first = mesh_boolean_fold(box, [drill, add], signs=[-1, +1])
    # union-then-subtract leaves the hole; subtract-then-union back-fills it -> larger volume
    assert subtract_first.volume > union_first.volume


def test_signs_length_mismatch_raises():
    from bertini_real.surface import mesh_boolean_fold
    box = trimesh.creation.box(extents=(1, 1, 1))
    with pytest.raises(ValueError):
        mesh_boolean_fold(box, [box], signs=[+1, -1])
