# Mathematical scope and replay argument

All integrals below are expectations under normalized unit-sphere area `dσ/(4π)`.
This common normalization cancels from every Rayleigh quotient. The accepted
operator is the self-adjoint form compression of `-Δ+q` to the real mean-zero
subspace of `H¹(S²)`, not the unrestricted Schrödinger operator.

## Moments and basis

The expectation of `x^a y^b z^c` is zero if any exponent is odd. For all-even
exponents it is

```
(a-1)!! (b-1)!! (c-1)!! / (a+b+c+1)!! .
```

This follows, for example, by writing an isotropic standard Gaussian in R³ as
its independent radius and uniform direction: its coordinate product moments
are the numerator, and its radial moment is the denominator. The convention
`(-1)!!=1` includes zero exponents. No quadrature nodes are used.

For each `l` and `0≤m≤l`, the real and imaginary parts of

```
(x+iy)^m Σ_{k=0}^{floor((l-m)/2)}
  (-1)^k (2l-2k)! z^(l-2k-m) (x²+y²+z²)^k
  / [2^l k! (l-k)! (l-2k-m)!]
```

give rational solid harmonics (the imaginary part is omitted when m=0).
This is the homogenized derivative of the Legendre polynomial. On the sphere
the real components agree with the standard associated-Legendre construction
up to nonzero normalization and a conventional sign; see
[NIST DLMF 14.30.1–3](https://dlmf.nist.gov/14.30).

The implementation checks that each constructed polynomial is homogeneous and
has zero Euclidean Laplacian. It computes Gram entries exactly, checks
orthogonality and positive diagonal mass, and independently computes

```
E[∇_S f · ∇_S g] = E[∇f · ∇g] - E[(x·∇f)(x·∇g)]
```

to check the energy `l(l+1)` against the mass. There are `2l+1` real harmonics
at each degree. Standard completeness of these harmonics in L²(S²) and their
Laplace eigenvalues are analytic dependencies of the certificate rules.

## Full operator and complete tail

Keep every harmonic of degrees `1,…,L`; let the orthogonal tail consist of every
harmonic of every degree `l≥L+1`. The degree-zero constant is absent from BOTH
parts because of the specified mean-zero compression.

Let `A` be the exact finite form, `M` its diagonal mass, and `B` the potential
coupling between retained and tail functions. A degree-d polynomial times a
retained harmonic restricts to a spherical polynomial of degree at most `L+d`.
Thus every finite-to-tail coupling above `L+d` vanishes. The code still bounds
ALL degrees above L; only their coupling columns vanish. All `2l+1` real
harmonics through `L+d` are explicitly enumerated, including zero columns.

The coefficient potential bound supplies globally valid `qlo≤q≤qhi` on S².
On the entire tail its compressed operator T therefore satisfies, in form order,

```
D_tail + qlo ≤ T ≤ D_tail + qhi .
```

For `x < (L+1)(L+2)+qlo`, all three shifted operators are positive. Inversion
reverses form order. Consequently the true Schur complement

```
S_true(x) = A-xM - B(T-x)^(-1)Bᵀ
```

satisfies `S_lower(x) ≤ S_true(x) ≤ S_upper(x)`, where

```
S_lower = A-xM - Σ_j B_j B_jᵀ / [mass_j (l_j(l_j+1)+qlo-x)]
S_upper = A-xM - Σ_j B_j B_jᵀ / [mass_j (l_j(l_j+1)+qhi-x)] .
```

Non-unit tail masses are essential: all basis functions are rationally scaled,
not normalized using square roots.

If `S_lower(x)` has no negative eigenvalues, the full shifted form is
nonnegative and its ground value is at least x. If `S_upper(x)` has a negative
or zero direction, so does the true Schur complement, proving ground value at
most x. A nonpositive direction of `A-xM` also gives a valid Ritz upper bound,
without requiring x below the tail threshold. Both search and replay use exact
rational inertia; the final replay additionally uses the prior dense inertia
implementation instead of the search's banded implementation.

## Sphere reduction and rotations

The polynomial division routine supplies the exact identity

```
q - reduced = (x²+y²+z²-1) * quotient .
```

It never assumes approximate vanishing. The range method chooses the original
or reduced polynomial using coefficient-range width, and the verifier rebuilds
the identity and selection. For `q=x²+y²+z²` the reduced potential is exactly 1,
so the mean-zero ground is exactly 3.

A proper rational rotation R is accepted only when its rows are exactly
orthonormal, its determinant is +1, and `q(x)=p((Rx)₃)` holds as a polynomial
identity. The induced sphere isometry preserves area, mean zero and the
Dirichlet form. Its spectrum therefore equals that of the axisymmetric
potential p(z). The transferred evidence is required to pass the V94 full-sphere
mean-zero verifier, including its separate angular and radial tail bounds.

## Limits

This is rational software replay, not a proof-assistant kernel. It does not
certify arbitrary manifolds, unconstrained higher full-sphere spectra, or a
universal nonlinear interpolation constant. A broad interval remains a valid
enclosure and is not re-labelled an accurate solution.
