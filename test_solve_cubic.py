import math
import cmath
import pytest
from cubic_solver import solve_cubic, ERROR_TOL

def roots_close(r1, r2, tol=ERROR_TOL):
    if len(r1) != len(r2):
        return False
    r1_sorted = sorted(r1, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    r2_sorted = sorted(r2, key=lambda x: (x.real, x.imag) if hasattr(x, 'real') else x)
    for a, b in zip(r1_sorted, r2_sorted):
        if abs(a - b) > tol:
            return False
    return True

@pytest.mark.parametrize("a, b, c, d, expected", [
    # Case 4a: q/h imaginary (should give three complex roots)
    (1, 0, 1, 0, [0, 1j, -1j]),
    # Case 4b: q/h real, q/h < -1 (should give three real roots)
    (1, -6, 11, -6, [1.0, 2.0, 3.0]),
    # Case 4c: q/h real, q/h > 1 (should give one real, two complex roots)
    (1, 0, 0, -1, [1.0, cmath.exp(2j * math.pi / 3), cmath.exp(-2j * math.pi / 3)]),
    # Case 4d: q/h real, |q/h| <= 1 (should give three real roots)
    (1, 0, -3, 2, [2.0, -1.0, -1.0]),
    # Degenerate cubic: triple root (h ≈ 0, p ≈ 0, q ≈ 0)
    (1, 0, 0, 0, [0.0, 0.0, 0.0]),
    # Degenerate cubic: three cube roots of -q/a, shifted (h ≈ 0, p ≈ 0, q ≠ 0)
    (1, 0, 0, 8, [-2.0, 1.0 + math.sqrt(3)*1j, 1.0 - math.sqrt(3)*1j]),
    # p < 0, q/h imaginary (should give three complex roots, explicit p < 0)
    (1, 0, 4, 0, [0.0, 2.0j, -2.0j]),
    # |q/h| <= 1, three real roots (explicit)
    (1, -3, -10, 24, [2.0, -3.0, 4.0]),
    # 1 real and 2 complex roots
    (1, -2, 9, -18, [2.0, -3.0j, 3.0j])
])
def test_solve_cubic(a, b, c, d, expected):
    result = solve_cubic(a, b, c, d)
    assert roots_close(result, expected), f"Failed for a={a}, b={b}, c={c}, d={d}: got {result}, expected {expected}"
