from cubic_solver import sqrt_trig, solve_cubic, solve_quadratic
import math, cmath 

ERROR_TOL = 1e-10

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
    """Solve a*x^4 + b*x^3 + c*x^2 + d*x + e = 0.
    Returns a list of 1..4 roots (real numbers or complex numbers).
    If the leading coefficients are zero the function will
    handle lower-degree polynomials automatically.
    """
    ##### Check degenerate cases
    if abs(a) < ERROR_TOL:  # Actully cubic
        return solve_cubic(b, c, d, e)
    
    if abs(e) < ERROR_TOL:  # Can factor out an x
        x1 = 0.
        other_roots = solve_cubic(a, b, c, d)
        return [x1] + other_roots 

    ### Depress the quartic: y^4 + py^2 + qy + r = 0
    ### and x = y + shift
    p = (8*a*c - 3*b**2) / (8 * a**2)
    q = (b**3 - 4*a*b*c + 8*(a**2)*d) / (8*a**3)
    r = (-3*b**4 + 256*(a**3)*e - 64*(a**2)*b*d + 16*a*(b**2)*c) / (256*a**4)
    shift = -b / (4 * a)


    ### Solve degenerate cases of depressed quartic 
    # q = 0 => Biquadratic form y^4 + py^2 + r = 0
    if abs(q) < ERROR_TOL:
        return solve_biquadratic(1, p, r, shift)
    

    ### Solve resolvent cubic and find one real root
    cubic_roots = solve_cubic(1, -p/2, -r, (p*r)/2 - (q**2)/8)
    real_root = None 
    for root in cubic_roots:
        if type(root) is float:
            real_root = root
            break 
        elif type(root) is complex and abs(root.imag) < ERROR_TOL:
            real_root = root.real 
            break 
    if real_root is None:
        raise RuntimeError("No real root of resolvent cubic found. Aborting.")

    ### Solve resulting quadratics to get roots of quartic
    ### y^4 + py^2 + qy + r = (y^2 + sqrt(2*z0)y + alpha)(y^2 - sqrt(2*z0)y + beta)
    ### where alpha = z0 - q/(2*sqrrt(2*z0))
    ### and   beta  = z0 + q/(2*sqrrt(2*z0))
    alpha = real_root - q / (2*sqrt_trig(2*real_root))
    beta = real_root + q / (2*sqrt_trig(2*real_root))
    y_vals_1 = solve_quadratic(1, sqrt_trig(2 * real_root), alpha)
    y_vals_2 = solve_quadratic(1, sqrt_trig(2 * real_root), beta)
    y_roots = y_vals_1 + y_vals_2
    x_roots = [y + shift for y in y_roots]
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
