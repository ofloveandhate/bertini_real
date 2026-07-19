"""
    :platform: Unix, Windows
    :synopsis: Plot decompositions DOWNSTAIRS -- in the projection plane.

The cellular decomposition is secretly the per-sheet lift of a planar
arrangement: the images of the boundary curves (critical, singular, sphere)
under (pi_0, pi_1) subdivide the projection plane, and the surface's cells lift
that subdivision sheet by sheet.  Apparent crossings are the arrangement's
vertices that live on no single curve; shadow curves are its edges lifted to
sheets that don't own them.  Looking at the arrangement itself -- several curve
decompositions overlaid in the (pi_0, pi_1)-plane, with the event vertices
marked -- is the fastest way to see whether a projection is 'hostile' (images
crossing region interiors) and where a decomposition needs events.

Typical use::

    import bertini_real as br
    from bertini_real.downstairs import plot_downstairs

    crit = br.curve.Curve('output_dim_1_comp_0')     # or any parsed curve(s)
    plot_downstairs([crit], projections=None)        # uses each curve's own pi
"""

# silviana amethyst + Claude
# 2026

import numpy as np

from bertini_real.vertex import VertexType


#: marker style per vertex type, painted in this order (later wins visually)
TYPE_MARKERS = [
    (VertexType.critical, dict(marker='.', color='black', ms=7, ls='none')),
    (VertexType.semicritical, dict(marker='o', mfc='none', color='tab:gray',
                                   ms=7, ls='none')),
    (VertexType.singular, dict(marker='s', mfc='none', color='magenta', ms=8,
                               ls='none')),
    (VertexType.distance_anchor, dict(marker='^', color='tab:green', ms=8,
                                      ls='none')),
    (VertexType.apparent_crossing, dict(marker='x', color='red', ms=11,
                                        mew=2.5, ls='none')),
]


def _real_point(vertex):
    return np.array([c.real if hasattr(c, 'real') else float(c)
                     for c in vertex.point], dtype=float)


def downstairs_coords(vertex, projections):
    """The (pi_0(x), pi_1(x)) coordinates of a vertex under coefficient-vector
    projections (constant term ignored; vertices are dehomogenized points)."""
    p = _real_point(vertex)
    return np.array([float(np.real(np.dot(pi[:len(p)], p)))
                     for pi in projections])


def _projections_for(decomposition, projections):
    if projections is not None:
        return [np.asarray([complex(c).real for c in pi], dtype=float)
                for pi in projections]
    # Decomposition.pi is stored VARIABLE-major by parse_decomposition --
    # pi[variable][projection_index], natural variables only (num_variables in
    # the decomp file already excludes the homogenizing coordinate).  BEWARE:
    # the parser pads pi to two columns regardless of dimension, so a curve
    # decomposition carries a phantom all-zero second projection -- only
    # `dimension` columns are real.  Embedded curves defer to their surface.
    src = decomposition
    if getattr(src, 'dimension', 0) < 2 \
            and getattr(src, 'embedded_into', None) is not None:
        src = src.embedded_into
    mat = np.asarray([[complex(c).real for c in np.ravel(row)]
                      for row in src.pi], dtype=float)
    if getattr(src, 'dimension', 0) < 2 or mat.ndim != 2 or mat.shape[1] < 2:
        raise ValueError(
            "need two projections for a downstairs plot; a curve decomposition "
            "stores only one -- pass projections=(pi0, pi1) explicitly")
    return [mat[:, j] for j in range(2)]


def edge_polyline(curve, edge_index, projections):
    """The downstairs polyline of one edge: left endpoint, samples (if the curve
    was sampled), midpoint, right endpoint."""
    e = curve.edges[edge_index]
    if e[0] == e[1] or e[1] == e[2]:
        return None                     # degenerate
    if getattr(curve, 'sampler_data', None):
        chain = curve.sampler_data[edge_index]
    else:
        chain = e
    pts = [downstairs_coords(curve.vertices[i], projections) for i in chain]
    return np.vstack(pts)


def plot_downstairs(decompositions, *, projections=None, ax=None, labels=None,
                    mark_types=True, legend=True, colors=None):
    """Overlay curve decompositions in the projection plane.

    :param decompositions: iterable of parsed curve decompositions (each with
        ``.vertices``, ``.edges``, optional ``.sampler_data``, and ``.pi``).
    :param projections: optional explicit pair of projection coefficient vectors
        (natural variables only).  Defaults to the first decomposition's own
        ``pi`` -- which works for surfaces' constituent curves; lone curve
        decompositions store a single projection and need this passed.
    :param ax: matplotlib axes to draw on (a new figure by default).
    :param labels: optional per-decomposition labels for the legend.
    :param mark_types: paint per-type vertex markers (critical dots,
        semicritical circles, singular squares, distance-anchor triangles,
        apparent-crossing X's).
    :rtype: the matplotlib axes.
    """
    import matplotlib.pyplot as plt

    decompositions = list(decompositions)
    if not decompositions:
        raise ValueError("nothing to plot")
    pis = _projections_for(decompositions[0], projections)

    if ax is None:
        _fig, ax = plt.subplots(figsize=(9, 9))
    if colors is None:
        colors = [f'C{k % 10}' for k in range(len(decompositions))]

    for k, curve in enumerate(decompositions):
        label = labels[k] if labels else getattr(curve, 'inputfilename', None)
        first = True
        for edge_index in range(len(curve.edges)):
            line = edge_polyline(curve, edge_index, pis)
            if line is None:
                continue
            ax.plot(line[:, 0], line[:, 1], '-', lw=1.3, color=colors[k],
                    label=(label if first else None))
            first = False

    if mark_types:
        seen = set()
        for curve in decompositions:
            for i, vertex in enumerate(curve.vertices):
                if id(vertex) in seen:
                    continue
                seen.add(id(vertex))
                for flag, style in TYPE_MARKERS:
                    if vertex.type & flag:
                        u = downstairs_coords(vertex, pis)
                        ax.plot([u[0]], [u[1]], **style)

    ax.set_xlabel(r'$\pi_0$')
    ax.set_ylabel(r'$\pi_1$')
    ax.set_aspect('equal')
    if legend and labels:
        ax.legend(loc='best', fontsize=9)
    return ax


def save_downstairs(decompositions, path, **kwargs):
    """:func:`plot_downstairs` straight to a file (Agg backend safe)."""
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    ax = plot_downstairs(decompositions, **kwargs)
    ax.figure.savefig(path, dpi=140, bbox_inches='tight')
    plt.close(ax.figure)
    return path
