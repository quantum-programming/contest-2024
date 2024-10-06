import argparse
import mpmath

from .ring import *
from .gridsynth import *

from .gridsynth import gridsynth, check


def main():
    mpmath.mp.dps = 100
    mpmath.mp.pretty = True

    parser = argparse.ArgumentParser()

    parser.add_argument('theta', type=str)
    parser.add_argument('epsilon', type=str)
    parser.add_argument('--verbose', '-v', action='store_true')
    parser.add_argument('--showgraph', '-g', action='store_true')

    args = parser.parse_args()
    theta = mpmath.mpmathify(args.theta)
    epsilon = mpmath.mpmathify(args.epsilon)

    sol = gridsynth(theta=theta, epsilon=epsilon, verbose=args.verbose, show_graph=args.showgraph)
    print(f"{theta=}, {epsilon=}")
    check(sol, theta)
