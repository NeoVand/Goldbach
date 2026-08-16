# What is known, what is predicted, and why the proof is hard

Every statement below is labeled **[theorem]**, **[heuristic]**, or
**[measurement]**. Keeping those straight is the entire point of this
repository.

## The conjecture

**Binary (strong) Goldbach.** Every even n > 2 is a sum of two primes.
Proposed in the Goldbach–Euler correspondence of 1742. Open.

**Ternary (weak) Goldbach.** Every odd n > 5 is a sum of three primes.
This one is a **[theorem]** — see below. The binary conjecture implies the
ternary trivially (n odd ⇒ n − 3 even); the reverse implication does not hold.

## What is actually proven

- **[theorem]** (Vinogradov, 1937) Every *sufficiently large* odd number is a
  sum of three primes — with an ineffective threshold later made explicit but
  astronomically large.
- **[theorem]** (Helfgott, 2013) Every odd number ≥ 7 is a sum of three
  primes: the ternary conjecture in full, by sharpening both major and minor
  arcs until they met a large computation in the middle.
- **[theorem]** (Chudakov; van der Corput; Estermann, 1937–38) *Almost all*
  even numbers are sums of two primes: the exceptional set up to X has density
  zero.
- **[theorem]** (Montgomery–Vaughan, 1975) The exceptional set is at most
  X^(1−δ) for an effective δ > 0; later work has pushed the exponent below
  3/4. If Goldbach fails, it fails on a thin set.
- **[theorem]** (Chen, 1973) Every sufficiently large even number is p + P₂:
  a prime plus a number with at most two prime factors. This is the
  high-water mark of sieve methods, and it is *exactly* at the parity barrier
  (below).
- **[theorem]** (Linnik 1953; Heath-Brown–Puchta 2002; Pintz–Ruzsa) Every
  sufficiently large even number is a sum of two primes and a bounded number
  of powers of two (13 suffice unconditionally; 8 under refinements).
- **[measurement]** (Oliveira e Silva–Herzog–Pardi, 2014) The conjecture holds
  for every even n ≤ 4·10¹⁸. (This repository independently re-verifies
  n ≤ 2·10⁷ as a by-product of computing every count exactly.)

## The Hardy–Littlewood model

**[heuristic]** Treat "m is prime" as a random event with probability
1/ln m (the prime number theorem's density), *then correct for arithmetic*.
For n = p + q, the naive independent model gives about
∫₂^(n−2) dt / (ln t · ln(n−t)) =: Li₂(n) ordered representations. But
primality is not independent across residue classes: for each odd prime ℓ,
the pair (p, n−p) must avoid ℓ | p and ℓ | n−p, and how constraining that is
depends on whether ℓ divides n. Working out the correction for each ℓ and
multiplying gives the **extended Goldbach conjecture** (Hardy–Littlewood,
1923):

    g(n)  ~  2·C₂ · S(n) · Li₂(n)

    C₂   = ∏ over odd primes (1 − 1/(ℓ−1)²) = 0.66016…   (twin prime constant)
    S(n) = ∏ over odd primes ℓ | n of (ℓ−1)/(ℓ−2)         (singular series)

S(n) ≥ 1 always: divisibility of n only ever *helps*. That is why the comet's
bands sit above the S(n) = 1 baseline, and why dividing by S(n) collapses
them (figure 2). The calibration measurement (figure 3) finds the model
accurate to ~1% at n ~ 10⁷ with an envelope shrinking like a square root —
the fluctuation size the model itself predicts.

## The circle method connection

Everything above has an exact analytic form. Let S(α) = Σ_{p ≤ N} e^{2πipα}.
Then, exactly,

    g(n) = ∫₀¹ |S(α)|² e^{−2πinα} dα.

The Goldbach count is a Fourier coefficient of the primes' power spectrum.
Our FFT instrument is this identity, discretized and made exact (figure 4
shows |S(α)| itself). The spectrum splits into:

- **Major arcs** — α near a rational a/q with small q. There S(α) is large
  and *understood*: |S(a/q)| ≈ π(N)·|μ(q)|/φ(q) (the measured spikes match
  to four decimals, and vanish for non-squarefree q exactly as the Möbius
  function demands). Integrating only the major arcs reproduces precisely
  the Hardy–Littlewood main term 2·C₂·S(n)·Li₂(n).
- **Minor arcs** — everything else. There S(α) exhibits square-root-size
  cancellation (the noise floor in figure 4), which is *believed* to make
  the minor-arc contribution negligible for every n, but is *proven* only
  in an averaged sense.

## Why the known methods stop short

**The L² barrier (why ternary fell and binary didn't).** By Parseval,
∫₀¹ |S(α)|² dα = π(N): the average size of |S|² is under complete control.
For the *ternary* problem the minor-arc error is bounded by
(max over minor arcs |S|) · ∫|S|², so one good *pointwise* bound on the
minor arcs (Vinogradov's) finishes the job. For the *binary* problem the
minor-arc error is ∫ over minor arcs of |S|² itself — Parseval says this can
be as large as ~π(N), which swamps the main term ~ n/ln²n for a *single* n.
To beat it you would need to know that |S|²'s mass on the minor arcs does not
conspire in phase against your particular n. Averaged over n it cannot (that
is exactly the almost-all theorem); for every individual n, nobody knows how.
The gap between "almost all" and "all" is the gap between the L² norm and a
pointwise bound, and it has not moved since 1938.

**The parity barrier (why sieves stop at Chen).** Selberg observed that sieve
methods — which only consume information about divisibility by small primes —
cannot distinguish integers with an odd number of prime factors from those
with an even number. Any sieve argument fine enough to produce p + q (two
primes: parity 1 × parity 1) would equally "produce" configurations it cannot
rule out, so unconditionally it can only ever reach p + P₂ — one factor of
ambiguity left, exactly Chen's theorem. Breaking parity requires an input
beyond divisibility data (bilinear/Type-II information); for Goldbach's
specific additive constraint, no one has found one.

**Why computation cannot finish.** Verification to 4·10¹⁸ plus an
exceptional-set bound of X^0.72-ish still leaves infinitely many candidates
above any checked range. A finite computation settles a finite set; the
conjecture is about all of them. What computation *can* do — what this
repository does — is measure whether the structural model that everyone
believes has any visible crack. (It has none.)

## What our measurements do and do not establish

They establish, up to 2·10⁷ with exact integer arithmetic:
that every even number has a partition; that partition counts match the 1923
model to ~1% with square-root-shrinking spread; that the model's arithmetic
factor explains all visible band structure; that the primes' spectrum has
its predicted major-arc heights and only square-root noise elsewhere; and
that worst cases stay glued to the predicted floor.

They establish **no theorem about any n beyond the ranges computed**, and no
amount of the same evidence could. A proof would need one of: pointwise
control of minor-arc mass for arbitrary individual n; a parity-breaking
input for the additive equation p + q = n; or an idea that belongs to
neither of the two known roads — which, given that both roads have been
walked to their exact ends, is where most number theorists expect it to
come from.

## Sources

- Hardy, Littlewood, *Some problems of 'Partitio Numerorum' III*, Acta Math. 44 (1923)
- Vinogradov, *Representation of an odd number as a sum of three primes* (1937)
- Estermann, *On Goldbach's problem*, Proc. LMS (1938)
- Montgomery, Vaughan, *The exceptional set in Goldbach's problem*, Acta Arith. 27 (1975)
- Chen, *On the representation of a larger even integer as the sum of a prime
  and the product of at most two primes*, Sci. Sinica 16 (1973)
- Heath-Brown, Puchta, *Integers represented as a sum of primes and powers of
  two*, Asian J. Math. 6 (2002)
- Helfgott, *The ternary Goldbach conjecture is true* (2013), arXiv:1312.7748
- Oliveira e Silva, Herzog, Pardi, *Empirical verification of the even
  Goldbach conjecture … up to 4·10¹⁸*, Math. Comp. 83 (2014)
- Friedlander, Iwaniec, *Opera de Cribro*, AMS Colloquium 57 (2010) — the
  parity problem, ch. 16
