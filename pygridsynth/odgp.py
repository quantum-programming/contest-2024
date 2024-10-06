from .mymath import SQRT2, floor, ceil, pow_sqrt2, floorlog
from .ring import ZRootTwo, DRootTwo, LAMBDA


def _solve_ODGP_internal(I, J):
    if I.width < 0 or J.width < 0:
        return []
    elif I.width > 0 and J.width <= 0:
        sol = _solve_ODGP_internal(J, I)
        return [beta.conj_sq2 for beta in sol]
    else:
        (n, _) = (0, 0) if J.width <= 0 else floorlog(J.width, LAMBDA.toReal)
        if n == 0:
            sol = []
            a_min = ceil((I.l + J.l) / 2)
            a_max = floor((I.r + J.r) / 2)
            for a in range(a_min, a_max + 1):
                b_min = ceil(SQRT2() * (a - J.r) / 2)
                b_max = floor(SQRT2() * (a - J.l) / 2)
                for b in range(b_min, b_max + 1):
                    sol.append(ZRootTwo(a, b))
            return sol
        else:
            lambda_n = LAMBDA ** n
            lambda_inv_n = LAMBDA ** -n
            lambda_conj_sq2_n = LAMBDA.conj_sq2 ** n
            sol = _solve_ODGP_internal(I * lambda_n.toReal, J * lambda_conj_sq2_n.toReal)
            sol = [beta * lambda_inv_n for beta in sol]
            return sol


def solve_ODGP(I, J):
    if I.width < 0 or J.width < 0:
        return []

    a = floor((I.l + J.l) / 2)
    b = floor(SQRT2() * (I.l - J.l) / 4)
    alpha = ZRootTwo(a, b)
    sol = _solve_ODGP_internal(I - alpha.toReal, J - alpha.conj_sq2.toReal)
    sol = [beta + alpha for beta in sol]
    sol = [beta for beta in sol if I.within(beta.toReal) and J.within(beta.conj_sq2.toReal)]
    return sol


def solve_ODGP_with_parity(I, J, beta):
    p = beta.parity
    sol = solve_ODGP((I - p) * SQRT2() / 2, (J - p) * (- SQRT2()) / 2)
    sol = [alpha * ZRootTwo(0, 1) + p for alpha in sol]
    return sol


def solve_scaled_ODGP(I, J, k):
    scale = pow_sqrt2(k)
    sol = solve_ODGP(I * scale, -J * scale if k & 1 else J * scale)
    return [DRootTwo(alpha, k) for alpha in sol]


def solve_scaled_ODGP_with_parity(I, J, k, beta):
    if k == 0:
        sol = solve_ODGP_with_parity(I, J, beta.renew_denomexp(0))
        return [DRootTwo.fromZRootTwo(alpha) for alpha in sol]
    else:
        p = beta.renew_denomexp(k).parity
        offset = DRootTwo.fromInt(0) if p == 0 else DRootTwo.pow_sqrt2_inv(k)
        sol = solve_scaled_ODGP(I - offset.toReal, J - offset.conj_sq2.toReal, k - 1)
        return [alpha + offset for alpha in sol]
