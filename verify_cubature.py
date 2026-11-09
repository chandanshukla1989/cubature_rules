"""Independent check of a cubature table for the uniform probability measure on [-1,1]^n.

Written from scratch; does not reuse the supplied checker.

Usage: verify_cub.py <table.txt> <n> <degree> [expected_nodes]

Table format: one node per row, n coordinates then the weight, '#' comments ignored.
Targets are exact rationals: for a monomial x^alpha, the average over [-1,1]^n is
0 when any alpha_i is odd, and prod 1/(alpha_i+1) otherwise.

Checks performed:
  * node count, distinctness, strict interiority |x_i| < 1, strict positivity of weights
  * weights sum to 1
  * EVERY monomial of total degree <= degree, with exact rational targets, in
    high-precision decimal arithmetic
  * the first failing degree (degree+1 and beyond), to confirm the rule is not
    accidentally of higher degree than claimed
  * central symmetry of the node set, reported but not required
"""
import itertools
import math
import sys
from decimal import Decimal, localcontext
from fractions import Fraction
from pathlib import Path


def load(path, n):
    rows = []
    for line in Path(path).read_text().splitlines():
        s = line.strip()
        if not s or s.startswith("#"):
            continue
        parts = s.split()
        if len(parts) != n + 1:
            raise SystemExit(f"row has {len(parts)} fields, expected {n+1}")
        rows.append(parts)
    return rows


def monomials(n, deg):
    for alpha in itertools.product(range(deg + 1), repeat=n):
        if sum(alpha) <= deg:
            yield alpha


def exact_target(alpha):
    if any(a % 2 for a in alpha):
        return Fraction(0)
    return Fraction(1, math.prod(a + 1 for a in alpha))


def residuals(rows, n, alphas, prec):
    """Max |computed - exact| over the given monomials."""
    with localcontext() as ctx:
        ctx.prec = prec
        nodes = [[Decimal(v) for v in r[:n]] for r in rows]
        wts = [Decimal(r[n]) for r in rows]
        worst, worst_alpha = Decimal(0), None
        for alpha in alphas:
            t = exact_target(alpha)
            target = Decimal(t.numerator) / Decimal(t.denominator)
            acc = Decimal(0)
            for node, w in zip(nodes, wts):
                term = w
                for x, a in zip(node, alpha):
                    if a:
                        term *= x ** a
                acc += term
            d = abs(acc - target)
            if d > worst:
                worst, worst_alpha = d, alpha
        return worst, worst_alpha


def main():
    path, n, deg = sys.argv[1], int(sys.argv[2]), int(sys.argv[3])
    expected = int(sys.argv[4]) if len(sys.argv) > 4 else None
    rows = load(path, n)
    N = len(rows)
    print(f"file                 {Path(path).name}")
    print(f"dimension            {n}")
    print(f"claimed degree       {deg}")
    print(f"nodes                {N}" + (f"   (expected {expected})" if expected else ""))
    if expected and N != expected:
        raise SystemExit("NODE COUNT MISMATCH")

    with localcontext() as ctx:
        ctx.prec = 140
        coords = [tuple(Decimal(v) for v in r[:n]) for r in rows]
        wts = [Decimal(r[n]) for r in rows]
        print(f"distinct nodes       {len(set(coords))}")
        assert len(set(coords)) == N, "duplicate nodes"
        assert all(w > 0 for w in wts), "a weight is not strictly positive"
        maxabs = max(abs(x) for c in coords for x in c)
        print(f"max |coordinate|     {maxabs}")
        assert maxabs < 1, "a node is on or outside the boundary"
        s = sum(wts)
        print(f"weights sum          {s}")
        print(f"weights sum - 1      {abs(s - 1)}")
        print(f"min weight           {min(wts)}")
        # central symmetry
        cset = set(coords)
        neg = {tuple(-x for x in c) for c in cset}
        sym = cset == neg
        at_origin = sum(1 for c in coords if all(x == 0 for x in c))
        print(f"centrally symmetric  {sym}   (nodes at origin: {at_origin})")

    alphas = list(monomials(n, deg))
    worst, wa = residuals(rows, n, alphas, 140)
    print(f"\nmonomials of degree <= {deg}: {len(alphas)}")
    print(f"max absolute residual        {worst}")
    print(f"attained at exponent vector  {wa}")
    ok = worst < Decimal("1e-100")
    print(f"EXACT TO DEGREE {deg}:         {'YES' if ok else 'NO'}")

    higher = [a for a in monomials(n, deg + 1) if sum(a) == deg + 1]
    w2, wa2 = residuals(rows, n, higher, 140)
    print(f"\nmonomials of degree exactly {deg+1}: {len(higher)}")
    print(f"max absolute residual              {w2}")
    print(f"exact at degree {deg+1} too?          "
          f"{'YES (rule is better than claimed)' if w2 < Decimal('1e-100') else 'no, as expected'}")

    print("\nVERDICT:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
