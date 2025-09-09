from cubic_solver import sqrt_trig, solve_cubic, solve_quadratic, is_real

ERROR_TOL = 1e-8

def solve_biquadratic(a, b, c, shift):
    """Solves the equation x = y + shift, where
    ay^4 + by^2 + c = 0
    """
    y_squared_vals = solve_quadratic(a, b, c)
    y_roots = []
    for val in y_squared_vals:
        sqrt_val = sqrt_trig(val)
        y_roots.append(sqrt_val)
        y_roots.append(-sqrt_val)
    x_roots = [y + shift for y in y_roots]
    return x_roots


def solve_quartic(a, b, c, d, e):
    # Degenerate case (a=0) => Really have a cubic equation
    if abs(a) < ERROR_TOL:
        return solve_cubic(b, c, d, e)

    # Depress quartic: x = y - b/(4a)
    p = (8*a*c - 3*b**2)/(8*a**2)
    q = (b**3 - 4*a*b*c + 8*a**2*d)/(8*a**3)
    r = (-3*b**4 + 256*a**3*e - 64*a**2*b*d + 16*a*b**2*c)/(256*a**4)
    shift = -b / (4 * a)
    y_roots = []

    if abs(q) < ERROR_TOL:
        x_roots = solve_biquadratic(1, p, r, shift)
    else:
        # General quartic: solve resolvent cubic
        cubic_roots = solve_cubic(1, -p/2, -r, p*r/2 - q**2/8)
        z0 = cubic_roots[0]  # pick one root
        u2 = 2*z0 - p        # This had issues before - I just had 2sqrt(2z0).
        u = sqrt_trig(u2)

        if abs(u) < ERROR_TOL:
            # u ≈ 0 => biquadratic fallback
            x_roots = solve_biquadratic(1, p, r, shift)
        else:
            y_roots += solve_quadratic(1, u, z0 - q/(2*u))
            y_roots += solve_quadratic(1, -u, z0 + q/(2*u))

            # Shift back
            x_roots = [y - b/(4*a) for y in y_roots]

    # Clamp tiny values to zero for numerical error
    x_roots = [0 if abs(x) < ERROR_TOL else x for x in x_roots]

    return x_roots


def main():
    tests = [
        (1, 0, 0, 0, -1),  # roots of x^4 - 1 = 0 (2 real and 2 complex roots)
        (1, 0, 1, 0, -1),  # roots of x^4 + x^2 - 1 = 0 (2 real and 2 complex roots)
        (0, 1, -3, 2, 0),  # a=0 => cubic: x^3 - 3x^2 + 2x = 0  (roots 0,1,2)
    ]
    for a, b, c, d, e in tests:
        roots = solve_quartic(a, b, c, d, e)
        print(f"solve_quartic({a}, {b}, {c}, {d}, {e}) -> {roots}")

if __name__ == "__main__":
    main()
