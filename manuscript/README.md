# Sharp fourth-order correlation-minor bounds without filter normal form

The note answers the conjectures immediately following Theorems 1 and 3,
Eqs. (10) and (12), of Peled, Te'eni, Carmi and Cohen, *Scientific Reports*
11, 2849 (2021), DOI 10.1038/s41598-021-82303-3.

## Files

- `note.tex`: complete mathematical exposition and embedded bibliography.
- `note.pdf`: seven-page typeset note.
- `verify.py`: exact rational checks using Python 3's standard library.
- `metadata.json`: author order, affiliations, license, problem identifier,
  and the fixed-revision formal source link.
- `LICENSE.txt`: CC BY 4.0 license notice.
- `verification.json`: exact-arithmetic results printed by `verify.py`.
- `source-manifest.json`: SHA256 and byte count for each accompanying file.

## Build

Use Tectonic 0.17.0:

```sh
tectonic --untrusted note.tex
```

Alternatively use a current LaTeX distribution with `pdflatex`:

```sh
pdflatex -interaction=nonstopmode -halt-on-error note.tex
pdflatex -interaction=nonstopmode -halt-on-error note.tex
```

The bibliography is included in `note.tex`, so BibTeX and external bibliography
files are unnecessary. Packages used are `geometry`, `amsmath`, `amssymb`,
`amsthm`, `booktabs`, `array`, `fontenc`, `lmodern`, and `hyperref`.

## Exact arithmetic

```sh
python3 verify.py
```

The program reconstructs the rational state at noise `1/1000`, checks its
trace and projectors, calculates the real correlation block, and evaluates
all 2,240 third-order minors of the 4-by-16 correlation matrix. Exactly four
entries in one compound row are nonzero. It checks the squared norm against
`9/8192` and the endpoint inequality for the interval `[0,1/200]`.
It also checks the tetrahedral sign identities and the fourth-order inequalities
at `D=5,6,7,8` and noise `0,1/1000`. The script prints a compact JSON summary.
It checks the sharp squared bound `1/110592` and weighted-Gram equality for
the noiseless states at `D=4,...,12`, the squared tail comparison in Appendix
D.5, Eq. (D31), and the positive Gram eigenvalues giving the coherent-design
cardinality obstruction. The general dimension statements follow from the
written proofs, not from these finite ranges.
These finite checks are linked to the all-parameter conclusions by the written
monotonicity and Bernoulli arguments in the note.

The SHA256 of the program's standard output is
`532cb325fdb2368101c380a91a695ffec100ccae5e7f8af6161bf0d033eb23a3`.

## Formal and written coverage

The formal source is
[`CorrelationMinorNormFNFRefutation.lean`](https://github.com/the-omega-institute/trureturing/blob/0f1f822c673c1903dc6b0a209a23dacebca2aa6c/D5/S3/Quantum/Entanglement/Witness/CorrelationMinorNormFNFRefutation.lean).
The declaration
`D5.S3.Quantum.Entanglement.Witness.CorrelationMinorNormFNFRefutation.result`
proves the negation of its universal no-FNF `claim` using the fixed rational
2-by-4 witness. The statement assumes local dimensions at least two and
`D*d <= h^2`, `h <= d^2`, separability, and complete Hilbert--Schmidt
orthonormal Hermitian bases.

The noise interval in Theorem 1, the rank-one compound formulas and Corollary 2,
the companion-conjecture result in Theorem 3, the sharp fourth-order bound
in Theorem 4, and the coherent-design, corrected-tail and equal-dimensional
CCNR Corollaries 5--7
are written proofs. The sharp maximum is `1/(192 sqrt(3))` for every `D>=4`
and every Schatten parameter. The broader design-characterization hypothesis
is distinguished from the proved nonexistence of coherent pairs.
In equal dimensions, CCNR implies the operator bound for `d<=h<=d^2`
without filter normal form; at `h=d^2` this holds for every Schatten parameter.
Li, Poon and Sze's 2011 Theorem 2.4 supplies the prior equal-dimensional
`p=1` comparison.
The original FNF theorems, coherent-design constructions, entropy results,
and discord results are not contradicted by the counterexamples.

The metadata retains the existing problem identifier. The companion conjecture,
the general-parameter design hypothesis, and Appendix D.5 are located explicitly
in the note's source-coverage table. The known covariance inequality is credited
to Zhang et al., Physical Review A 77, 060301(R) (2008), Theorem 1, and
Gittsovich et al., Physical Review A 78, 052319 (2008), Proposition IV.2.

## Authors and license

Haobo Ma — ChronoAI Pte Ltd, Singapore.

Wenlin Zhang — National University of Singapore; ChronoAI Pte Ltd, Singapore.

Copyright 2026 Haobo Ma and Wenlin Zhang.
Creative Commons Attribution 4.0 International (CC BY 4.0; `cc-by-4.0`).
See `LICENSE.txt` and https://creativecommons.org/licenses/by/4.0/.
