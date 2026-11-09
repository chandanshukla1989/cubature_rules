#!/usr/bin/env python3
"""Moller's lower bound on the number of nodes of a cubature rule of odd degree
for a centrally symmetric weight function.

Moller (1979) proved that a cubature rule of degree 2k+1 for a centrally
symmetric integral in n variables needs at least

    N >= 2 * dim P_k^odd              when k is odd,
    N >= 2 * dim P_k^even - 1         when k is even,

where P_k^odd and P_k^even are the spaces of odd, resp. even, polynomials of
total degree at most k in n variables. The two dimensions are counted here by
summing the number of monomials of each admissible total degree.

For degree 5 this reduces to the familiar n^2 + n + 1.

Usage:
    python3 moller_bound.py            # table for small n and odd degrees
    python3 moller_bound.py 5 5        # the bound for n = 5, degree 5

Standard library only.
"""
import sys
from math import comb


def monomials(n, d):
    """Number of monomials in n variables of total degree exactly d."""
    return comb(n + d - 1, d)


def dim_parity(n, k, parity):
    """dim of the span of monomials of total degree <= k whose degree has the
    given parity (0 = even, 1 = odd)."""
    return sum(monomials(n, d) for d in range(k + 1) if d % 2 == parity)


def moller_bound(n, degree):
    """Lower bound on the node count for a degree-`degree` rule in n variables,
    for a centrally symmetric weight. `degree` must be odd."""
    if degree % 2 == 0:
        raise ValueError("Moller's bound in this form is stated for odd degree")
    k = (degree - 1) // 2
    if k % 2 == 1:
        return 2 * dim_parity(n, k, 1)
    return 2 * dim_parity(n, k, 0) - 1


def main():
    if len(sys.argv) == 3:
        n, d = int(sys.argv[1]), int(sys.argv[2])
        k = (d - 1) // 2
        print(f"n = {n}, degree = {d}  (k = {k}, k is {'odd' if k % 2 else 'even'})")
        print(f"Moller lower bound: {moller_bound(n, d)}")
        return

    degrees = [3, 5, 7, 9]
    print("Moller lower bound on the number of nodes, centrally symmetric weight")
    print("(the cube [-1,1]^n is centrally symmetric, so the bound applies to it)\n")
    header = "  n  " + "".join(f"  degree {d:>2}" for d in degrees)
    print(header)
    print("  " + "-" * (len(header) - 2))
    for n in range(2, 9):
        row = f"  {n:<3}" + "".join(f"{moller_bound(n, d):>11}" for d in degrees)
        print(row)
    print("\nThe two rules in this repository, against the bound:")
    for n, d, N in ((5, 5, 31), (4, 7, 52)):
        b = moller_bound(n, d)
        verdict = "ATTAINS THE BOUND (minimal)" if N == b else f"{N - b} above the bound"
        print(f"  cube n={n}, degree {d}: {N} nodes, bound {b}  ->  {verdict}")


if __name__ == "__main__":
    main()
