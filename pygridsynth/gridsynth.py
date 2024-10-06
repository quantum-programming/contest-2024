import mpmath
import time

from .mymath import sqrt, solve_quadratic
from .ring import DRootTwo
from .region import Ellipse, ConvexSet
from .to_upright import to_upright_set_pair
from .tdgp import solve_TDGP
from .diophantine import NO_SOLUTION, diophantine_dyadic


class EpsilonRegion(ConvexSet):
    def __init__(self, theta, epsilon):
        self._theta = theta
        self._epsilon = epsilon
        self._d = 1 - epsilon ** 2 / 2
        self._z_x = mpmath.cos(-theta / 2)
        self._z_y = mpmath.sin(-theta / 2)
        D_1 = mpmath.matrix([[self._z_x, -self._z_y], [self._z_y, self._z_x]])
        D_2 = mpmath.matrix([[4 * (1 / epsilon) ** 4, 0], [0, (1 / epsilon) ** 2]])
        D_3 = mpmath.matrix([[self._z_x, self._z_y], [-self._z_y, self._z_x]])
        p = mpmath.matrix([self._d * self._z_x, self._d * self._z_y])
        ellipse = Ellipse(D_1 * D_2 * D_3, p)
        super().__init__(ellipse)

    @property
    def theta(self):
        return self._theta

    @property
    def epsilon(self):
        return self._epsilon

    def inside(self, u):
        cos_similarity = (self._z_x * u.real + self._z_y * u.imag)
        return DRootTwo.fromDOmega(u.conj * u) <= 1 and cos_similarity >= self._d

    def intersect(self, u0, v):
        a = v.conj * v
        b = 2 * v.conj * u0
        c = u0.conj * u0 - 1

        vz = self._z_x * v.real + self._z_y * v.imag
        rhs = self._d - self._z_x * u0.real - self._z_y * u0.imag
        t = solve_quadratic(a.real, b.real, c.real)
        if t is None:
            return None
        t0, t1 = t
        if vz > 0:
            t2 = rhs / vz
            return (t0, t1) if t0 > t2 else (t2, t1)
        elif vz < 0:
            t2 = rhs / vz
            return (t0, t1) if t1 < t2 else (t0, t2)
        else:
            return (t0, t1) if rhs <= 0 else None


class UnitDisk(ConvexSet):
    def __init__(self):
        ellipse = Ellipse(mpmath.matrix([[1, 0], [0, 1]]), mpmath.matrix([0, 0]))
        super().__init__(ellipse)

    def inside(self, u):
        return DRootTwo.fromDOmega(u.conj * u) <= 1

    def intersect(self, u0, v):
        a = v.conj * v
        b = 2 * v.conj * u0
        c = u0.conj * u0 - 1
        return solve_quadratic(a.real, b.real, c.real)


def generate_unitary(sol):
    u, t = sol
    return mpmath.matrix([[u.toComplex, -t.conj.toComplex],
                          [t.toComplex, u.conj.toComplex]])


def generate_target_Rz(theta):
    return mpmath.matrix([[mpmath.exp(- 1.j * theta / 2), 0],
                          [0, mpmath.exp(1.j * theta / 2)]])


def error(theta, sol):
    Rz = generate_target_Rz(theta)
    U = generate_unitary(sol)
    E = U - Rz
    return sqrt(mpmath.fabs(E[0, 0] * E[1, 1] + E[0, 1] * E[1, 0]))


def check(sol, theta):
    print(f"u={sol[0]}, t={sol[1]}")
    Rz = generate_target_Rz(theta)
    U = generate_unitary(sol)
    e = error(theta, sol)
    print(f"{Rz=}")
    print(f"{U=}")
    print(f"{e=}")


def gridsynth(theta, epsilon, verbose=False, show_graph=False):
    epsilon_region = EpsilonRegion(theta, epsilon)
    unit_disk = UnitDisk()
    k = 0

    if verbose:
        start = time.time()
    transformed = to_upright_set_pair(epsilon_region, unit_disk, verbose=verbose, show_graph=show_graph)
    if verbose:
        print(f"to_upright_set_pair: {time.time() - start} s")
        print("------------------")

    while True:
        if verbose:
            print(f"trial {k=}:")
            start = time.time()
        sol = solve_TDGP(epsilon_region, unit_disk, *transformed, k,
                         verbose=verbose, show_graph=show_graph)
        if verbose:
            print(f"solve_TDGP: {time.time() - start} s")
            start = time.time()

        for u in sol:
            xi = 1 - DRootTwo.fromDOmega(u.conj * u)
            t = diophantine_dyadic(xi)
            if t != NO_SOLUTION:
                sol = u, t
                if verbose:
                    print(f"diophantine_dyadic: {time.time() - start} s")
                return sol
        if verbose:
            print(f"diophantine_dyadic: {time.time() - start} s")
            print("------------------")
        k += 1
