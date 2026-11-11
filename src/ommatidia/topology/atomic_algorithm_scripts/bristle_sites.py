"""Atomic: choose which cell vertices (facet corners) carry an interommatidial bristle.

In Drosophila there is about one bristle per ommatidium ("An anatomical atlas
of Drosophila melanogaster - the wild-type"), at every other corner of each
hexagon. A hexagonal lattice has two kinds of triangles of centres, pointing
"up" and "down", and adjacent triangles are always of opposite kinds. Taking
all triangles of one kind gives exactly every other corner, and there are about
2 triangles per centre, so about one bristle per ommatidium.
"""

from __future__ import annotations

from collections import deque

import numpy as np

# DD. TRIANGLE_PARITY
# triangle_parity = np.ndarray, shape (T,), int8, values in {0, 1}
# interp. 2-colouring of triangles so that neighbours differ. Exact on a perfect
#         lattice; around 5/7 defects an odd cycle forces some same-colour
#         neighbours (PARITY_CONFLICTS)
triangle_parity0 = np.zeros(0, dtype=np.int8)  # stub
triangle_parity1 = np.array([0, 1, 0, 1], dtype=np.int8)

# DD. BRISTLE_MASK
# bristle_mask = np.ndarray, shape (T,), bool
# interp. True where the cell vertex of that triangle carries a bristle

GREEDY_SWEEPS = 20  # passes of local conflict repair after the BFS colouring


def _neighbour_lists(n_triangles: int, cell_edges: np.ndarray) -> list[list[int]]:
    nbrs: list[list[int]] = [[] for _ in range(n_triangles)]
    for a, b in cell_edges.tolist():
        nbrs[a].append(b)
        nbrs[b].append(a)
    return nbrs


# DD. SEAM_CUTS
# seam_cuts = np.ndarray, shape (E,), bool, aligned with CELL_EDGES
# interp. cell edges whose two corners may share a colour. Every cell with an odd
#         number of sides (5, 7) is an odd cycle that cannot alternate, so a
#         "seam" of such edges must end at it. Short seams pairing nearby
#         defects (or running to an open rim) keep the pattern regular elsewhere.
seam_cuts0 = np.zeros(0, dtype=bool)  # stub


# FD. seam_cuts
# Signature: (E, 2) crossed centre pairs, (N,) bool odd cells, (N,) bool rim cells -> SEAM_CUTS
def seam_cuts(crossed: np.ndarray, odd: np.ndarray, rim: np.ndarray) -> np.ndarray:
    """Pair odd cells greedily by short paths through the lattice (an approximate
    minimum T-join).

    Algorithm:
        for each odd cell v not yet paired:
            breadth-first search over centre-to-centre edges from v until the
                first other unpaired odd cell, or a rim cell
            toggle every cell edge crossing that path (toggling keeps the seam
                count odd at both path ends and even at every other cell)
    """
    nbrs: list[list[tuple[int, int]]] = [[] for _ in range(len(odd))]
    for e, (a, b) in enumerate(crossed.tolist()):
        nbrs[a].append((b, e))
        nbrs[b].append((a, e))
    cut = np.zeros(len(crossed), dtype=bool)
    unpaired = set(np.flatnonzero(odd).tolist())
    for v in np.flatnonzero(odd).tolist():
        if v not in unpaired:
            continue
        unpaired.discard(v)
        came_from = {v: (-1, -1)}
        queue = deque([v])
        target = -1
        while queue:
            x = queue.popleft()
            if x != v and (x in unpaired or rim[x]):
                target = x
                break
            for y, e in nbrs[x]:
                if y not in came_from:
                    came_from[y] = (x, e)
                    queue.append(y)
        if target < 0:
            continue
        unpaired.discard(target)
        x = target
        while x != v:
            x, e = came_from[x]
            cut[e] = not cut[e]
    return cut


# FD. triangle_parity
# Signature: int (T), CELL_EDGES, SEAM_CUTS -> TRIANGLE_PARITY
def triangle_parity(n_triangles: int, cell_edges: np.ndarray,
                    cuts: np.ndarray | None = None) -> np.ndarray:
    """2-colour the triangles so that neighbours differ wherever possible.

    Algorithm:
        breadth-first search over the dual graph minus the seam cuts (bipartite
            once the seams pair every odd cell), giving each newly reached
            triangle the opposite colour of the one it came from
        repeat a few sweeps: flip any triangle that has more same-colour
            neighbours than different-colour ones (each flip strictly removes
            conflicts, so this terminates)
    """
    parity = np.full(n_triangles, -1, dtype=np.int8)
    nbrs = _neighbour_lists(n_triangles, cell_edges)
    keep = np.ones(len(cell_edges), dtype=bool) if cuts is None else ~cuts
    walk = _neighbour_lists(n_triangles, cell_edges[keep])
    for root in range(n_triangles):
        if parity[root] >= 0:
            continue
        parity[root] = 0
        queue = deque([root])
        while queue:
            t = queue.popleft()
            for u in walk[t]:
                if parity[u] < 0:
                    parity[u] = 1 - parity[t]
                    queue.append(u)
    if len(cell_edges) == 0:
        return parity
    a, b = cell_edges[:, 0], cell_edges[:, 1]
    for _ in range(GREEDY_SWEEPS):
        same = parity[a] == parity[b]
        n_same = np.bincount(a[same], minlength=n_triangles) + np.bincount(b[same], minlength=n_triangles)
        degree = np.bincount(a, minlength=n_triangles) + np.bincount(b, minlength=n_triangles)
        flip = 2 * n_same > degree
        if not flip.any():
            break
        # Flip one at a time, re-checking against the current colours, so that
        # two adjacent candidates never undo each other's improvement.
        for t in np.flatnonzero(flip):
            if 2 * sum(parity[u] == parity[t] for u in nbrs[t]) > len(nbrs[t]):
                parity[t] = 1 - parity[t]
    return parity


# FD. count_conflicts
# Signature: TRIANGLE_PARITY, CELL_EDGES -> int
def count_conflicts(parity: np.ndarray, cell_edges: np.ndarray) -> int:
    """Number of adjacent triangle pairs with the same colour."""
    if len(cell_edges) == 0:
        return 0
    return int(np.sum(parity[cell_edges[:, 0]] == parity[cell_edges[:, 1]]))


# FD. select_bristles
# Signature: TRIANGLE_PARITY, BRISTLE_PATTERN, BRISTLE_FREQUENCY, RNG -> BRISTLE_MASK
def select_bristles(parity: np.ndarray, pattern: str, frequency: float,
                    rng: np.random.Generator) -> np.ndarray:
    """Mark the cell vertices that carry a bristle.

    Algorithm:
        ALTERNATING: the triangles of colour 0 (every other corner of each hexagon)
        RANDOM:      each triangle independently with probability BRISTLE_FREQUENCY
        NONE:        none
    """
    if pattern == "ALTERNATING":
        return parity == 0
    if pattern == "RANDOM":
        return rng.random(len(parity)) < frequency
    if pattern == "NONE":
        return np.zeros(len(parity), dtype=bool)
    raise ValueError(f"Unknown bristle pattern {pattern!r}")
