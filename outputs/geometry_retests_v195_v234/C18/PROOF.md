# Exact reflection-block proof obligations

The operator is the self-adjoint form compression of `-Δ_S²+q` to all real
mean-zero H¹ functions on the unit sphere. Normalized area `dσ/(4π)` is used
throughout, so no area constants are lost in Rayleigh quotients.

## Exact global potential range

Polynomial division verifies `q-q_eff=(x²+y²+z²-1)h`. For reflection-even
`q_eff`, substitute `X=x²,Y=y²,Z=z²`. These coordinates fill the simplex
`X,Y,Z≥0, X+Y+Z=1`. A polynomial in X,Y,Z of degree d is homogenized by
multiplying each degree-k monomial by `(X+Y+Z)^(d-k)`.

Write the resulting polynomial as

```
Σ_{a+b+c=d} β_abc * d!/(a!b!c!) * X^a Y^b Z^c .
```

The nonnegative Bernstein weights sum to 1. The minimum and maximum exact
rational coefficients β therefore bound the potential globally. For C18 the
three coefficients are 1,2,3, giving `[1,3]`. The technique also applies to the
supported even quartic potentials, with d≤2 in squared coordinates. If the
input is not even, this range routine returns a safe coefficient bound and
does not authorize reflection splitting.

## Complete character decomposition

The eight sign-flip maps form the finite group `(Z/2Z)^3`. Its character p in
`{0,1}³` has projector

```
Π_p u(x) = (1/8) Σ_{f∈{0,1}³} (-1)^(p·f) u((-1)^f x).
```

For a polynomial this selects monomials whose coordinate exponents have
parity p. The eight projectors are mutually orthogonal and sum to identity.
Sign flips preserve spherical area and the Dirichlet form. If the exact
potential identity verifies invariance under every flip, the compressed
operator commutes with every projector. Consequently it is the orthogonal
direct sum of the eight character operators. The zero-mean projection
commutes with these maps too.

We use real solid spherical harmonics obtained from homogenized associated
Legendre derivatives. Their on-sphere construction and normalization conventions
are standard; see [NIST DLMF §14.30](https://dlmf.nist.gov/14.30). Rational
normalizations avoid square roots; all masses and integrals are retained
exactly. Each solid polynomial is checked to be harmonic, and each block is
checked for positive diagonal mass and orthogonality.

Let k be the number of odd coordinates in p. At degree l the dimension of
the p-character harmonic subspace is `(l-k)/2+1` when l≥k and l−k is even,
and zero otherwise. Indeed, parity-restricted homogeneous polynomials have
dimension `C(n+2,2)` with n=(l-k)/2; subtracting their r²-multiples of degree
l−2 removes `C(n+1,2)` dimensions. Summing the eight character dimensions
gives `2l+1`. The implementation checks this count through the supported
construction and the tests verify the combined degrees independently.

For p=000 the constant is removed, so the first degree is 2. For every other
character its first degree is k. The corresponding unperturbed floors are
6 for 000, 2 for the three single-odd characters, 6 for the three double-odd
characters, and 12 for 111. Their exact potential lower shift supplies a
valid bound for the entire character, including every degree never assembled.

## Same-character infinite tail

Within each block the harmonic degrees advance by 2. After all degrees ≤L
are kept, the first omitted degree is the next admissible degree of that
character, not necessarily L+1. Let l_last denote the actual last kept degree.
A potential of degree d cannot couple a kept harmonic to a degree beyond
l_last+d. Every same-character harmonic up to that bound is explicitly
included in the coupling columns. Every higher degree remains covered by the
analytic tail bound even though its finite coupling vanishes.

For a spectral test point x smaller than the first omitted Laplace eigenvalue
plus qlo, the true tail operator shifted by x is strictly positive. Bounds
`D_tail+qlo≤T≤D_tail+qhi` and inverse-order monotonicity give

```
S_lower(x) ≤ S_true(x) ≤ S_upper(x),
S_bound(x) = A-xM - Σ_j B_j B_jᵀ /
            [mass_j * (l_j(l_j+1)+q_bound-x)].
```

The qlo comparison is the lower Schur form; qhi is the upper form. No negative
eigenvalues of the lower form prove a block lower spectral endpoint. A
nonpositive direction of the upper form proves an upper endpoint; a Ritz
direction in A−xM can also prove the upper endpoint. These tests include
non-unit masses and use exact rational inertia. Replay rebuilds the complete
form and tail and uses the previous dense inertia routine independently of
the search's banded routine.

## All-character assembly and adaptive acceptance

The global lower endpoint is the minimum of eight valid character lower
bounds. The upper endpoint is the minimum of the computed block upper bounds.
Uncomputed characters MUST remain in the certificate with a verified analytic
floor. Eight entries in fixed character order are required; missing, duplicated,
wrong-input or wrong-scope entries fail verification.

The driver can refine only a character whose lower bound still competes with
the current global upper bound. Every other proof remains in the result.
Increasing a degree uses the next same-character degree. Bit increases and
degree increases are distinguished in the trace; the separate diagnosis
actually runs both experiments.

The requested maximum retained degree also applies during initial coverage.
When a required character would exceed that cap, it remains represented by
its analytic floor and the driver returns an open gap. Wall time is checked
between indivisible exact operations, never by interrupting a proof halfway;
the reported total includes final replay and discloses any overrun.

Acceptance depends on the original potential, full mean-zero scope and exact
requested width bound in the final certificate. Trace, cache counts and wall
time are diagnostic measurements, not substitutes for mathematical proof.

This is rational software verification relying on standard spherical harmonic
completeness, min–max and Schur principles. It is not a proof-assistant kernel
and makes no claim about arbitrary manifolds or non-reflection-invariant
potentials.
