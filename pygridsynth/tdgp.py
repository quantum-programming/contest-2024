from .ring import DRootTwo, DOmega
from .region import Interval
from .odgp import solve_scaled_ODGP, solve_scaled_ODGP_with_parity
from .myplot import plot_sol


def solve_TDGP(setA, setB, opG, ellipseA_upright, ellipseB_upright, bboxA, bboxB, k,
               verbose=False, show_graph=False):
    sol_sufficient = []
    sol_x = solve_scaled_ODGP(bboxA.I_x, bboxB.I_x, k + 1)
    sol_y = solve_scaled_ODGP(bboxA.I_y.fatten(bboxA.I_y.width * 1e-4),
                              bboxB.I_y.fatten(bboxB.I_y.width * 1e-4),
                              k + 1)
    if len(sol_x) <= 0 or len(sol_y) <= 0:
        sol_sufficient = []
    else:
        alpha0 = sol_x[0]
        for beta in sol_y:
            dx = DRootTwo.pow_sqrt2_inv(k)
            u0 = opG.inv * DOmega.fromDRootTwoVector(alpha0, beta, k + 1)
            v = opG.inv * DOmega.fromDRootTwoVector(dx, DRootTwo.fromInt(0), k)
            t_A = setA.intersect(u0, v)
            t_B = setB.intersect(u0.conj_sq2, v.conj_sq2)
            if t_A is None or t_B is None:
                continue

            parity = (beta - alpha0).mul_sqrt2_renewing_denomexp(k)
            intA, intB = Interval(*t_A), Interval(*t_B)
            dtA = 10 / max(10, (1 << k) * intB.width)
            dtB = 10 / max(10, (1 << k) * intA.width)
            intA, intB = intA.fatten(dtA), intB.fatten(dtB)
            sol_t = solve_scaled_ODGP_with_parity(intA, intB, 1, parity)
            sol_x = [alpha * dx + alpha0 for alpha in sol_t]
            for alpha in sol_x:
                sol_sufficient.append(DOmega.fromDRootTwoVector(alpha, beta, k))
    sol_transformed = [opG.inv * u for u in sol_sufficient]
    sol = [u for u in sol_transformed if setA.inside(u) and setB.inside(u.conj_sq2)]

    if verbose:
        print(f"sol_sufficient size: {len(sol_sufficient)}, sol size: {len(sol)},")
    if show_graph and len(sol_sufficient) > 0:
        plot_sol([sol_transformed, sol], setA.ellipse, setB.ellipse, None, None,
                 color_list=['limegreen', 'blue'], size_list=[5, 10])

    return sol
