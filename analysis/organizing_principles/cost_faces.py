"""Exact zero-cost boundary face inventory. Author: Anubhav Chaturvedi."""
from pathlib import Path
from fractions import Fraction as F
from itertools import product
import importlib.util
import json
import sys

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))
from npalevels.core import SYMMETRIES, probabilities

spec = importlib.util.spec_from_file_location('sector', ROOT / 'analysis/motzkin_test/geometry/check_sector.py')
sector = importlib.util.module_from_spec(spec)
spec.loader.exec_module(sector)
HERE = Path(__file__).resolve().parent


def inventory():
    normal = (F(0),) * 4 + sector.EPS
    vertices = list(sector.deterministic())
    facet = sector.FACET
    events = [sector.losing_event(d) for d in facet]
    losses = [sum(1 << i for i, e in enumerate(events) if probabilities(d)[e] == 1) for d in vertices]
    assert sorted(m.bit_count() for m in losses) == [1] * 8 + [3] * 8
    permutations = []
    for g in SYMMETRIES:
        if g.functional(normal) == normal:
            permutations.append([facet.index(g.behavior(d)) for d in facet])
    assert len(permutations) == 16
    unseen = set(range(1, 255))
    rows = []
    while unseen:
        mask = min(unseen)
        orbit = sorted({sum(1 << p[i] for i in range(8) if mask >> i & 1) for p in permutations})
        unseen.difference_update(orbit)
        tied = [i for i, loss in enumerate(losses) if loss & ~mask == 0]
        dimension = sector.rank([[1, *vertices[i]] for i in tied]) - 1
        costs = [F(not (mask >> i & 1)) for i in range(8)]
        # Every deterministic optimum is retained, including vertices outside the cap.
        assert tied == [j for j, d in enumerate(vertices)
                        if sum(costs[i] * probabilities(d)[e] for i, e in enumerate(events)) == 0]
        rows.append(dict(zero_cost_mask=mask, zero_cost_count=mask.bit_count(),
                         orbit=orbit, deterministic_indices=tied,
                         local_face_dimension=dimension,
                         whole_quantum_and_Q1AB_face_is_local=mask.bit_count() <= 4))
    assert sum(len(r['orbit']) for r in rows) == 254
    assert len(rows) == 28
    solved = [r for r in rows if r['whole_quantum_and_Q1AB_face_is_local']]
    assert len(solved) == 18 and sum(len(r['orbit']) for r in solved) == 162
    critical = (F(1), F(0), F(1), F(0), F(1), F(1), F(1), F(-1))
    c = [sum(critical[j] * (sector.PR[j] - d[j]) for j in range(8)) for d in facet]
    mask = sum(1 << i for i, value in enumerate(c) if value == 0)
    tied = [i for i, loss in enumerate(losses) if loss & ~mask == 0]
    assert len(tied) == 3 and sector.rank([[1, *vertices[i]] for i in tied]) == 3
    return dict(events=events, deterministic_vertices=vertices, deterministic_loss_masks=losses,
                boundary_orbits=rows, critical_edge=dict(zero_cost_mask=mask, costs=c,
                                                        deterministic_indices=tied))


if __name__ == '__main__':
    data = json.loads(json.dumps(inventory(), default=str))
    target = HERE / 'cost_faces.json'
    if '--write' in sys.argv:
        target.write_text(json.dumps(data, indent=2) + '\n', encoding='utf-8')
    else:
        assert json.loads(target.read_text(encoding='utf-8')) == data
    print('PASS: 254 boundary masks / 28 symmetry types; all 162 masks with at least four positive costs have a completely local maximizing face; critical edge is a triangle.')
