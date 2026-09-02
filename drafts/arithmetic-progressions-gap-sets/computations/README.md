# Computational certificates

This directory contains the computational certificate material for
*Arithmetic Progressions in Gap Sets of Numerical Semigroups*.

The manuscript source remains one level above this directory:

```text
../main.tex
../references.bib
```

The certificate package is organized as follows:

```text
computations/
├── README.md
├── CERTIFICATE.txt
├── SHA256SUMS
├── generate_certificate.sh
├── src/
│   ├── plateau_certificate.c
│   ├── plateau_certificate_independent.c
│   ├── direct_ap_search.c
│   ├── staircase_criterion.c
│   └── staircase_chunk.c
└── logs/
    ├── plateau_certificate.log
    ├── plateau_certificate_independent.log
    ├── direct_ap_tests.log
    ├── direct_ap_extremizer9.log
    ├── direct_ap_A8.log
    ├── direct_ap_table_12.log
    ├── staircase_criterion_1e8.log
    ├── staircase_chunk_sample.log
    └── staircase_chunk_1e8.log
```

The five `.c` files in `src/` are the computational source code. The retained
execution logs are in `logs/`. The script `generate_certificate.sh` does **not**
create the source files: it validates the retained logs and regenerates only
the concise summary `CERTIFICATE.txt`.

## File map

- `src/plateau_certificate.c` — primary finite Kunz optimization for $A=9$ and
  $m=8,9,10$.
- `src/plateau_certificate_independent.c` — independent reordered
  branch-and-bound implementation of the same optimization.
- `src/direct_ap_search.c` — direct construction of gap sets and direct
  computation of arithmetic progressions; it verifies $M(8)=48$, the small
  table through $A=12$, the multiplicity-$9$ extremizer, and regression
  examples.
- `src/staircase_criterion.c` — exact composite-loss criterion through
  $A=10^8$.
- `src/staircase_chunk.c` — the same criterion on an inclusive interval, used as
  a separate full-range cross-check.
- `generate_certificate.sh` — checks all PASS markers, compares the two
  full-range criterion outputs, and regenerates `CERTIFICATE.txt`.
- `CERTIFICATE.txt` — short human-readable summary generated from the logs.
- `SHA256SUMS` — checksums for the manuscript and computational package.

The retained execution logs are

```text
logs/plateau_certificate.log
logs/plateau_certificate_independent.log
logs/direct_ap_tests.log
logs/direct_ap_A8.log
logs/direct_ap_table_12.log
logs/direct_ap_extremizer9.log
logs/staircase_criterion_1e8.log
logs/staircase_chunk_sample.log
logs/staircase_chunk_1e8.log
```

Compiled binaries are temporary and belong in `build/`; they are not part of
the certificate package.

## Mathematical scope

The plateau programs impose exactly the ordinary and wrapped Kunz
inequalities and the Run-Lemma consequences displayed as equations
(13)--(17) in the manuscript. These orbit conditions are necessary for

\[
\operatorname{AP}(S)\leq 9,
\]

not sufficient. This relaxation nevertheless gives a valid upper bound,
because every relevant semigroup belongs to the enumerated set.

The direct program is logically independent: during branch-and-bound it builds
the genuine gap prefixes already forced by assigned Kunz coordinates and prunes
a branch only when those actual gaps already contain an arithmetic progression
longer than A.  Any surviving complete candidate is checked again from its full
gap set.  No Run-Lemma constraint is used. The
criterion programs do not enumerate numerical semigroups; they factor each
candidate multiplicity and test the exact inequality from Proposition 6.1.
In their failure records, `p` always means $P(A+1)$; no `q` is printed because
a failure means that every prime divisor $q\mid m$ failed the criterion.

## Compilation

Run these commands from this `computations/` directory:

```sh
mkdir -p build
cc -O3 -std=c11 -Wall -Wextra -Wpedantic \
  src/plateau_certificate.c -o build/plateau_certificate
cc -O3 -std=c11 -Wall -Wextra -Wpedantic \
  src/plateau_certificate_independent.c -o build/plateau_certificate_independent
cc -O3 -std=c11 -Wall -Wextra -Wpedantic \
  src/direct_ap_search.c -o build/direct_ap_search
cc -O3 -std=c11 -Wall -Wextra -Wpedantic \
  src/staircase_criterion.c -o build/staircase_criterion
cc -O3 -std=c11 -Wall -Wextra -Wpedantic \
  src/staircase_chunk.c -o build/staircase_chunk
```

The compiler should emit no warnings.

## Reproduction commands

```sh
build/plateau_certificate > logs/plateau_certificate.log
build/plateau_certificate_independent > logs/plateau_certificate_independent.log

build/direct_ap_search tests > logs/direct_ap_tests.log
build/direct_ap_search extremizer9 > logs/direct_ap_extremizer9.log
build/direct_ap_search a8 > logs/direct_ap_A8.log
build/direct_ap_search table 12 > logs/direct_ap_table_12.log

build/staircase_criterion 100000000 > logs/staircase_criterion_1e8.log
build/staircase_chunk 1 100000 > logs/staircase_chunk_sample.log
build/staircase_chunk 1 100000000 > logs/staircase_chunk_1e8.log

sh generate_certificate.sh
```

The complete criterion runs allocate approximately 200 MB. Their running
times depend on the machine.

## Expected results

The two plateau implementations must independently return

```text
m=8   maximum=44   vector=(7,9,2,9,4,4,9)
m=9   maximum=54   vector=(3,3,6,6,9,9,9,9)
m=10  maximum=53   vector=(9,9,9,8,6,4,4,3,1)
```

Each maximum is unique within the stated relaxation. The direct search gives
$M(8)=48$. For the multiplicity-$9$ vector, it gives

```text
genus=54  F=80  AP=9  maximizing differences=3,9
```

The direct table is

```text
A     1  2  3   4   5   6   7   8   9   10   11   12
M(A)  1  4  7  16  20  36  42  48  54  100  110  144
```

Both complete criterion programs must report (1,449,077,567) tested pairs
and precisely the four failures

```text
(A,m)=(3,4), (8,9), (9,9), (9,10).
```

Each principal program exits with a nonzero status if its advertised
assertions fail. In addition, `generate_certificate.sh` refuses to produce a
summary unless the required logs exist and agree.

## Scope of the theorem through \(4\cdot10^{18}\)

No program here exhausts every $A\leq4\cdot10^{18}$. The local computation
ends at $10^8$. Above that point, Theorem 6.2 uses the published computation
of consecutive-prime gaps together with the mathematical prime-gap criterion.
The external prime-gap tables are not reproduced in this package.

## Robustness notes

- Kunz coordinates are indexed by `1,...,m-1`; the case `i+j==m` is skipped.
- Wrapped inequalities use target `i+j-m` and the carry `+1`.
- The direct checker never marks `0` as a gap and scans every
  `1 <= d <= F`.  Its `direct_AP_checks` counter includes direct checks on
  forced partial gap sets as well as checks of surviving complete candidates.
- Gap storage is dynamically allocated through `F`; no fixed-width bit mask
  is used, so positions above `127` are handled correctly.
- The criterion programs enforce $A\leq10^8$. Products are evaluated in
  `uint64_t`, and the smallest-prime-factor table is valid throughout this
  range.
- Allocation failures produce a nonzero exit status.

## Checksums

After regenerating all logs and `CERTIFICATE.txt`, run

```sh
sha256sum \
  ../main.tex ../references.bib README.md CERTIFICATE.txt \
  generate_certificate.sh src/*.c \
  logs/plateau_certificate.log logs/plateau_certificate_independent.log \
  logs/direct_ap_tests.log logs/direct_ap_A8.log logs/direct_ap_table_12.log \
  logs/direct_ap_extremizer9.log logs/staircase_criterion_1e8.log \
  logs/staircase_chunk_sample.log logs/staircase_chunk_1e8.log \
  > SHA256SUMS
```

Verify the package with

```sh
sha256sum -c SHA256SUMS
```

Checksums establish byte identity; they do not replace source inspection or
independent execution.
