"""Unused prototype for repeated serialization of integral routing matrices.

Construction and retained memory must be charged in any future integration.
No active collector imports this module.
"""
import math
from numbers import Integral


class MatrixTextCache:
    def __init__(self, costs):
        assert len(costs.shape) == 2 and costs.shape[0] == costs.shape[1]
        self.size = costs.shape[0]
        rows = []
        for row in costs:
            tokens = []
            for value in row:
                assert math.isfinite(value) and value == int(value)
                tokens.append(str(int(value)))
            rows.append(tuple(tokens))
        self.rows = tuple(rows)

    def serialize(self, support):
        support = list(support)
        assert all(isinstance(v, Integral) and not isinstance(v, bool) for v in support)
        assert len(support) >= 3 and len(set(support)) == len(support)
        assert all(0 <= v < self.size for v in support)
        chosen = sorted(support)
        lines = ['NAME: PCCTP_SUPPORT', 'TYPE: TSP', f'DIMENSION: {len(chosen)}',
                 'EDGE_WEIGHT_TYPE: EXPLICIT', 'EDGE_WEIGHT_FORMAT: FULL_MATRIX',
                 'EDGE_WEIGHT_SECTION']
        lines.extend(' '.join(self.rows[a][b] for b in chosen) for a in chosen)
        return '\n'.join(lines+['EOF'])+'\n'
