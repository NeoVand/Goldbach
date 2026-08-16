# Research log

An ongoing hunt for something genuinely new about Goldbach's conjecture.
Protocol per iteration: **hypothesis → instrument → experiment → control →
verdict**, with the verdict stating plainly what is known theory, what is new
measurement, and what (if anything) is a new phenomenon. Nothing gets called
a discovery until it survives a control and a literature-awareness check.

---

## Iteration 1 — The Riemann zeros are audible in Goldbach data

**Hypothesis.** Fujii's explicit formula (1991) says the Λ-weighted cumulative
Goldbach count obeys Σ_{n≤X} G(n) = X²/2 − 2Σ_ρ X^(ρ+1)/(ρ(ρ+1)) + error.
If real, the rescaled deviation D(X) = (Σ G(n) − X²/2)/X^(3/2) must be a
*superposition of log-periodic tones*: angular frequency γ (each zeta zero's
imaginary part), amplitude 4/|ρ(ρ+1)|. Falsifiable prediction: a periodogram
of measured D in u = log X shows peaks at 14.13, 21.02, 25.01, … with those
exact heights, and the same pipeline on a structureless prime-like set shows
no such tones.

**Instrument.** `goldbach/weighted.py`: Λ-weighted FFT autocorrelation → G(n)
for all n ≤ 2·10⁷ in one shot; cumulative deviation; Hann-windowed
periodogram in log X. Control: Cramér random set (weight ln n with
probability 1/ln n, matching E[Λ] exactly), same pipeline, seed 42.

**Results** (N = 2·10⁷, window X ∈ [10⁴, 2·10⁷], resolution 0.83 rad;
`results/riemann_zeros.json`, `figures/07_riemann_zeros.png`):

| zero γ | peak found at | measured amp | Fujii prediction | agreement |
|---|---|---|---|---|
| 14.135 | 14.049 | 0.01976 | 0.01990 | 0.7% |
| 21.022 | 20.661 | 0.00798 | 0.00903 | 12% |
| 25.011 | 24.793 | 0.00615 | 0.00638 | 4% |
| 30.425 | 30.578 | 0.00421 | 0.00432 | 3% |
| 32.935 | 33.057 | 0.00368 | 0.00368 | 0.1% |

- Mean |location offset| of first five zeros: 0.19 rad (resolution 0.83).
- Template-shift test: sliding the 5-zero template by δ ∈ [−6, 6] (0.05 grid,
  |δ| > 2×resolution excluded), **no shift beats δ = 0**: empirical p = 0.
- γ₁ stands **126×** above the primes' own periodogram noise floor.
- The Cramér control is ~**23× louder** than the primes' floor across the
  band, and its peaks align with nothing.

**Verdict.** Hypothesis confirmed, cleanly. Two readings of what happened:

1. *Known theory, made empirical.* The explicit formula itself is Fujii's;
   detecting the zeros in Goldbach data confirms known mathematics rather
   than revealing unknown mathematics.
2. *The surprise is the control.* The random model — which matches the primes
   in density exactly — produces a deviation signal an order of magnitude
   LOUDER than the real primes. The primes are not "random plus noise"; they
   are *quieter than random*, and essentially all of the residual sound is
   the zeta zeros. In this precise, measurable sense: **the Goldbach error
   is not noise, it is music** — and the instrument (one FFT) hears it at
   amplitudes matching an RH-conditional formula to within ~1% for γ₁ and
   γ₅.

Novelty assessment (honest): the theorem is known; explicit-formula numerics
exist for ψ(x). I am not aware of a published periodogram *detection with a
Cramér control and amplitude-level verification* run on Goldbach counts
specifically — but absence-of-awareness is not absence, and any write-up must
scope the claim as an experimental/methodological contribution, not new
theory. Worth keeping as a paper candidate.

**Next hypotheses queued.**
- H2: pair correlation of Goldbach fluctuations — does cov(g(n), g(n+h))
  follow a singular-series law in h (the additive analogue of prime pair
  correlation)? Instrument exists (FFT of residuals).
- H3: Chebyshev-type bias inside partitions — are the primes used by n's
  partitions equidistributed mod 4, or biased?
- H4: extreme-value law of the dyadic worst cases — Gumbel?
- H5: detection-threshold scaling — how does the number of resolvable zeros
  grow with N? (relates window length to zero density)
