import math
import cmath
import pytest
from quartic_solver import solve_quartic, ERROR_TOL

def roots_close(r1, r2, tol=ERROR_TOL):
    if len(r1) != len(r2):
        return False
    r1_sorted = sorted(r1, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    r2_sorted = sorted(r2, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    for a, b in zip(r1_sorted, r2_sorted):
        if abs(a - b) > tol:
            return False
    return True

@pytest.mark.parametrize("a, b, c, d, e, expected", [
    # x^4 - 1 = 0 => roots: 1, -1, i, -i
    (1, 0, 0, 0, -1, [1.0, -1.0, 1j, -1j]),
    # x^4 + x^2 - 1 = 0 => roots: approx 0.786, -0.786, 1.272j, -1.272j
    (1, 0, 1, 0, -1, [0.7861513777574233, -0.7861513777574233, 1.27201964951j, -1.27201964951j]),
    # Degenerate: a=0, cubic x^3 - 3x^2 + 2x = 0 => roots: 0, 1, 2
    (0, 1, -3, 2, 0, [0.0, 1.0, 2.0]),
    # x^4 + 2x^2 + 1 = 0 => roots: i, -i, i, -i (double roots)
    (1, 0, 2, 0, 1, [1j, -1j, 1j, -1j]),
    # x^4 + 4x^2 + 4 = 0 => roots: sqrt(2)i, -sqrt(2)i, sqrt(2)i, -sqrt(2)i (double roots)
    (1, 0, 4, 0, 4, [1.41421356237j, -1.41421356237j, 1.41421356237j, -1.41421356237j]),
    # x^4 = 0 => roots: 0, 0, 0, 0
    (1, 0, 0, 0, 0, [0.0, 0.0, 0.0, 0.0]),
])
def test_solve_quartic(a, b, c, d, e, expected):
    result = solve_quartic(a, b, c, d, e)
    assert roots_close(result, expected), f"Failed for a={a}, b={b}, c={c}, d={d}, e={e}: got {result}, expected {expected}"
