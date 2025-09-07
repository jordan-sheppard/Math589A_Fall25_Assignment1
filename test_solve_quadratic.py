import math
import pytest
from cubic_solver import solve_quadratic, ERROR_TOL

def roots_close(r1, r2, tol=ERROR_TOL):
    if len(r1) != len(r2):
        return False
    r1_sorted = sorted(r1, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    r2_sorted = sorted(r2, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    for a, b in zip(r1_sorted, r2_sorted):
        if abs(a - b) > tol:
            return False
    return True

@pytest.mark.parametrize("a, b, c, expected", [
    # Linear case: a = 0, bx + c = 0
    (0, 2, -4, [2.0]),
    (0, 0, 5, []),  # No solution
    (0, 0, 0, []),  # No equation
    # Quadratic with two real roots
    (1, -3, 2, [1.0, 2.0]),
    # Quadratic with one real repeated root
    (1, -2, 1, [1.0, 1.0]),
    # Quadratic with complex roots
    (1, 0, 1, [1j, -1j]),
    # Quadratic with two real roots, negative c
    (1, 1, -6, [2.0, -3.0]),
    # Quadratic with large coefficients
    (1e6, -3e6, 2e6, [1.0, 2.0]),
    # Quadratic with small coefficients
    (1e-6, -3e-6, 2e-6, [1.0, 2.0]),
    # Quadratic with imaginary coefficients 
    (1, -2j, -1, [1j, 1j])
])
def test_solve_quadratic(a, b, c, expected):
    result = solve_quadratic(a, b, c)
    assert roots_close(result, expected), f"Failed for a={a}, b={b}, c={c}: got {result}, expected {expected}"
