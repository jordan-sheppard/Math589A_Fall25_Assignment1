import math, cmath

ERROR_TOL = 1e-7

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


def solve_degenerate_cubic(a, q, shift):
    """Solves for x, where 
    x = y + shift, 
    and y satisfies the degenerate cubic equation
    ay^3 + q = 0.
    """
    if abs(q - 0.) < ERROR_TOL:     # We have a*y^3 = 0 (3 repeated real roots)
        y = 0 
        x = y + shift 
        return [x, x, x]
    else:                           # We have a*y^3 = -q, or y^3 = (-q/a)
        rhs = complex(-q / a)
        mag = abs(rhs)                                      # Complex radius of (-q/a)
        arg = math.atan2(rhs.imag, rhs.real)                # Complex angle of (-q/a)
        cbrt_mag = mag**(1/3)                               # Cube root magnitude of (-q/a)
        cbrt_args = [arg + 2*k*math.pi/3 for k in range(3)] # Cube root angles (there are 3) of (-q/a)
        y_roots = [cbrt_mag * cmath.exp(1j * theta) for theta in cbrt_args]
        x_roots = [y + shift for y in y_roots]
        return x_roots
    

def solve_cubic_imaginary(p:float, check_val:complex, shift:float):
    D = math.sqrt(-p)           # If this is called, p < 0
    phi = math.asinh(check_val.imag)  
    sinh_tmp_val = math.sinh(phi / 3)
    cosh_tmp_val = 1j * math.sqrt(3) * math.cosh(phi / 3)
    y1 = -2 * D * sinh_tmp_val 
    y2 = D * (sinh_tmp_val + cosh_tmp_val)
    y3 = D * (sinh_tmp_val - cosh_tmp_val)
    x_roots = [y1 + shift, y2 + shift, y3 + shift]
    return x_roots 


def solve_cubic_cosine(p:float, check_val:float, shift:float):
    phi = math.acos(check_val)
    sqrt_p = math.sqrt(p)
    y1 = 2 * sqrt_p * math.cos((phi + 0) / 3)
    y2 = 2 * sqrt_p * math.cos((phi - 2*math.pi) / 3)
    y3 = 2 * sqrt_p * math.cos((phi - 4*math.pi) / 3)

    x_roots = [y1 + shift, y2 + shift, y3 + shift]
    return x_roots

def solve_cubic_cosh(p:float, check_val:float, shift:float, negative:bool=False):
    factor = -1 if negative else 1
    phi = math.acosh(check_val)
    cosh_tmp_val = math.cosh(phi / 3)
    sinh_tmp_val = 1j * math.sqrt(3) * math.sinh(phi / 3)
    sqrt_p = math.sqrt(p)
    y1 = factor * 2 * sqrt_p * cosh_tmp_val
    y2 = factor * sqrt_p * (cosh_tmp_val + sinh_tmp_val)
    y3 = factor * sqrt_p * (cosh_tmp_val - sinh_tmp_val)
    x_roots = [y1 + shift, y2 + shift, y3 + shift]
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

    ### Step 2: Depress the cubic to a*y^3 - 3a * p * y + q = 0
    ### where y = x - shift
    shift = -b / (3 * a)
    p = (b**2 - (3 * a * c)) / (9 * a**2)
    q = a*shift**3 + b*shift**2 + c*shift + d
    h = 2 * a * p**(3/2)
    
    #######     ------ SOLVE VARIOUS CASES ------     #######
    ### Case A: h = 0 (that is, p=0; no linear term in depressed cubic)
    ### => Degenerate Cubic
    if abs(h - 0.) < ERROR_TOL:
        print("A")
        return solve_degenerate_cubic(a, q, shift)

    ### Case B: p < 0 (that is, q/h imaginary)
    ### => One real and 2 complex conjugate roots 
    check_val = q / h     # Rest of cases end up using this
    if p < 0:    
        print("B")
        return solve_cubic_imaginary(p, check_val, shift)
    
    else:
        ### Case C: |q/h| <= 1 
        ### => Three real roots
        if abs(check_val) <= 1:
            print("C")
            return solve_cubic_cosine(p, check_val, shift)
        
        ### Case D: q/h < -1 
        ### => One real root, two complex roots
        elif check_val < -1:
            print("D")
            return solve_cubic_cosh(p, -check_val, shift, negative=True)
        
        ### Case E: q/h > 1
        ### => One real root, two complex roots
        else:
            print("E")
            return solve_cubic_cosh(p, check_val, shift)
        


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