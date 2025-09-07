import math, cmath

ERROR_TOL = 1e-10

def solve_linear(a, b):
    """Solves the linear equation
    ax + b = 0.
    If a=0, 
    """
    # 1) Not linear (a = 0) => No Equation to Solve at All!
    if abs(a - 0.) < ERROR_TOL:
        return []
    # 2) Linear (a != 0) => Easy to solve for x
    else:
        return [-b / a] 

def solve_quadratic_real_cosine(gamma, shift):
    """Solves the equation
    x = y - shift,
    where
    2y^2 - 1 = gamma,
    and
    |gamma| <= 1.
    This is done using the trig substition
    y = cos(theta)
    leading to the form
    cos(2theta) = gamma.
    This leads to two real roots x_1, x_2, that are not
    repeated unless gamma = 1, -1.
    """
    theta_1 = 0.5 * math.acos(gamma)
    theta_2 = theta_1 + math.pi 

    x_1 = math.cos(theta_1) - shift 
    x_2 = math.cos(theta_2) - shift
    return [x_1, x_2]


def solve_quadratic_real_cosh(gamma, shift):
    """Solves the equation
    x = y - shift,
    where
    2y^2 - 1 = gamma,
    and
    gamma > 1.
    This is done using the trig substition
    y = cos(theta)
    where theta = iu for some real u.
    This leads eading to the form
    cos(2iu) = cosh(2u) = gamma.
    This leads to a repeated real root x_1 = x_2. 
    """
    u = 0.5 * math.acosh(gamma)
    y = math.cosh(u)
    x1 = y - shift
    x2 = -y - shift
    return [x1, x2]


def solve_quadratic_complex_conjugate(gamma, shift):
    """Solves the equation
    x = y - shift,
    where
    2y^2 - 1 = gamma,
    and
    gamma < -1.
    This is done using the trig substition
    y = cos(theta)
    where theta = iu + pi for some real u.
    This leads eading to the form
    cos(2iu + pi) = -cos(2iu) = -cosh(2u) = gamma,
    in other words,
    cosh(2u) = -gamma > 1
    This leads to two roots that are complex conjugates
    x_1 = x_2^*.
    """
    u = 0.5 * math.acosh(-gamma)
    y = 1j * math.sinh(u)
    x1 = -shift + y 
    x2 = -shift - y
    return [x1, x2]


def solve_quadratic(a, b, c):
    """Solves the quadratic equation 
    ax^2 + bx + c = 0
    Returns a list of 1..2 roots (complex if needed)
    """
    ### Step 1: Solve degenerate cases - linear or below
    if abs(a - 0.) < ERROR_TOL:
        return solve_linear(b, c)
    
    ### Step 2: Get constants from depressed quadratic
    shift = b / (2 * a)
    beta = (b**2 - 4 * a * c) / (4 * (a**2))
    gamma = (2 * beta) - 1

    if gamma < -1:
        return solve_quadratic_complex_conjugate(gamma, shift)
    elif gamma > 1:
        return solve_quadratic_real_cosh(gamma, shift)
    else:
        return solve_quadratic_real_cosine(gamma, shift)


def solve_degenerate_cubic_no_p(q, shift):
    """Solves for x, where 
    x = y + shift, 
    and y satisfies the degenerate cubic equation
    y^3 + q = 0.
    """
    if abs(q) < ERROR_TOL:     # We have y^3 = 0 (3 repeated real roots)
        return [shift, shift, shift]
    else:                           # We have a*y^3 = -q, or y^3 = (-q/a)
        rhs = complex(-q)
        mag = abs(rhs)                                      # Complex radius of (-q/a)
        arg = math.atan2(rhs.imag, rhs.real)                # Complex angle of (-q/a)
        cbrt_mag = mag**(1/3)                               # Cube root magnitude of (-q/a)
        cbrt_args = [arg + 2*k*math.pi/3 for k in range(3)] # Cube root angles (there are 3) of (-q/a)
        y_roots = [cbrt_mag * cmath.exp(1j * theta) for theta in cbrt_args]
        x_roots = [y + shift for y in y_roots]
        return x_roots


def solve_degenrate_cubic_no_q(p, shift):
    """Solves for x, where
    x = y + shift
    and y satisfies the degenerate cubic equation
    y^3 + py = 0
    """
    y1 = 0.
    sqrt = cmath.sqrt if p > 0 else math.sqrt
    y2 = sqrt(-p)
    y3 = -sqrt(-p)
    x_roots = [y1 + shift, y2 + shift, y3 + shift]
    return x_roots


def solve_cubic_3_real_roots(p, q, shift):
    """Solves for x when there are 3 distinct 
    real roots
    (that is, when the discriminant is negative).
    
    Uses the formula
    y_k = 2 sqrt(-p/3) * cos(1/3 * acos(3q/2p * sqrt(-3/p)-2k*pi/3))
    and 
    x_k = y_k + shift
    """
    inside_val = 1/3 * math.acos((3*q)/(2*p) * math.sqrt(-3/p))
    y_roots = [
        2 * math.sqrt(-p/3) 
        * math.cos(inside_val - (2 * math.pi * k) / 3)
        for k in range(3)
    ]
    x_roots = [
        root + shift for root in y_roots
    ]
    return x_roots

def solve_cubic_some_complex_roots(p, q, shift):
    """Solves the equation
    x = y + shift
    where 
    y^3 + py + q = 0
    and the discriminant is > 0 by using hyperbolic trig
    substitution.
    """
    # Get first root using hyperbolic trig substitution
    if p > 0:
        inside = math.asinh(
            (3*q)/(2*p) * math.sqrt(3/p)
        ) / 3
        outside = -2 * math.sqrt(p / 3)
        y_1 = outside * math.sinh(inside)
    else:
        inside = math.acosh(
            (3*abs(q))/(2*abs(p)) * math.sqrt(3/abs(p))
        ) / 3
        sign_q = -1 if q < 0 else 1     # q=0 case already handled
        outside = 2 * sign_q * math.sqrt(-p/3)
        y_1 = outside * math.cosh(inside)
    
    # Get other two roots by factoring out (y-y_1):
    # y^3 + py + q = (y-y_1)(ay^2 + by + c)
    # Done using synthetic division, and using a quadratic solver
    # on the quadratic portion
    a = 1
    b = y_1
    c = p + y_1**2
    other_roots = solve_quadratic(a, b, c) 
    y_roots = [y_1] + other_roots 
    x_roots = [y + shift for y in y_roots]
    return x_roots

def solve_cubic(a, b, c, d):
    """Solve a*x^3 + b*x^2 + c*x + d = 0
    Returns list of 1..3 roots (complex if needed).

    Uses a method given by:
    G.C Holmes, "The use of hyperbolic cosines in solving cubic polynomials",
        The Mathematical Gazette, Vol. 86, No. 507 (Nov. 2002), pp. 473-477
    """
    ### Step 1: Check degenerate case - not cubic at all
    if abs(a - 0.) < ERROR_TOL:
        return solve_quadratic(b, c, d)

    ### Step 2: Depress the cubic to y^3 + py + q = 0
    ### where x = y + shift
    # 2a) Normalize to monic: x^3 + A x^2 + B x + C = 0
    A = b / a
    B = c / a
    C = d / a

    shift = -A / 3
    p = B - (A**2 / 3)
    q = ((2 * A**3) / 27) - ((A * B) / 3) + C
    discriminant = (q / 2)**2 + (p / 3)**3

    ####### ------ SOLVE VARIOUS CASES ------ #######
    ### Case A: p = 0 => No linear term in depressed cubic
    if abs(p) < ERROR_TOL:
        print("A")
        return solve_degenerate_cubic_no_p(q, shift)
    
    ### Case B: q = 0 => No constant term in depressed cubic
    elif abs(q) < ERROR_TOL:
        print("B")
        return solve_degenrate_cubic_no_q(p, shift)

    ### Case C: discriminant <= 0
    ### => 3 Distinct Real Roots
    elif discriminant <= 0:
        print("C")
        return solve_cubic_3_real_roots(p, q, shift)
    
    ### Case D: discriminant > 0
    ### => Some complex roots
    else:
        print("D")
        return solve_cubic_some_complex_roots(p, q, shift)
        

def main():
    tests = [
        (1, 0, 0, -1),     # roots of x^3 - 1 = 0 (1 and two complex cube roots)
        (1, -6, 11, -6),   # roots [1.0, 2.0, 3.0]
    ]
    for a, b, c, d in tests:
        roots = solve_cubic(a, b, c, d)
        print(f"solve_cubic({a}, {b}, {c}, {d}) -> {roots}")

if __name__ == "__main__":
    # print(solve_cubic(1, 1, -11, 45))         # Case E
    # print(solve_cubic(1, -5, 19, -15))        # Case B 
    print(solve_cubic(1, -11, 35, -49))         # Case D