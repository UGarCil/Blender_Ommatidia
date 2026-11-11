"""Atomic: the ordered ring of corners around each facet, from the triangulation."""

from __future__ import annotations

import numpy as np

# DD. CELL_LOOPS
# cell_loops = list[np.ndarray (m,) int64] (one per centre; empty = no facet)
# interp. the triangles (= facet corners, rows of TOPOLOGY.cell_points) around
#         each centre, counter-clockwise about its normal. A facet exists only
#         where the loop closes: rim cells and cells next to a hole get none.
# Invariant: consecutive corners (and the last and first) share a cell edge
cell_loops0 = []                                     # stub
cell_loops1 = [np.array([0, 1, 2, 3, 4, 5]), np.zeros(0, dtype=np.int64)]


# FD. cell_loops
# Signature: TRIANGLES, int -> CELL_LOOPS
def cell_loops(triangles: np.ndarray, n: int) -> list:
    """Walk around each centre from triangle to triangle.

    Algorithm:
        for each triangle (a, b, c) (counter-clockwise), note for corner a:
            "after edge a→b comes this triangle, which leaves along a→c"
        for each centre i: start at any triangle, follow the a→c edge into the
            triangle that starts with it, until back at the start
        a missing triangle (rim or hole) or a repeated edge -> no loop
    """
    nxt = {}
    first = {}
    bad = set()
    for t, (a, b, c) in enumerate(np.asarray(triangles).tolist()):
        for i, j, k in ((a, b, c), (b, c, a), (c, a, b)):
            if (i, j) in nxt:
                bad.add(i)          # two triangles start with the same edge: overlap
            nxt[(i, j)] = (k, t)
            first.setdefault(i, j)
    loops = []
    for i in range(n):
        if i not in first or i in bad:
            loops.append(np.zeros(0, dtype=np.int64))
            continue
        start = first[i]
        j, loop = start, []
        while True:
            step = nxt.get((i, j))
            if step is None or len(loop) > 12:
                loop = []
                break
            j, t = step
            loop.append(t)
            if j == start:
                break
        loops.append(np.asarray(loop, dtype=np.int64))
    return loops
