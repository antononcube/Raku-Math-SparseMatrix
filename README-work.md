# Math::SparseMatrix


[![Actions Status](https://github.com/antononcube/Raku-Math-SparseMatrix/actions/workflows/linux.yml/badge.svg)](https://github.com/antononcube/Raku-Math-SparseMatrix/actions)
[![Actions Status](https://github.com/antononcube/Raku-Math-SparseMatrix/actions/workflows/macos.yml/badge.svg)](https://github.com/antononcube/Raku-Math-SparseMatrix/actions)
[![Actions Status](https://github.com/antononcube/Raku-Math-SparseMatrix/actions/workflows/windows.yml/badge.svg)](https://github.com/antononcube/Raku-Math-SparseMatrix/actions)

[![](https://raku.land/zef:antononcube/Math::SparseMatrix/badges/version)](https://raku.land/zef:antononcube/Math::SparseMatrix)
[![License: Artistic-2.0](https://img.shields.io/badge/License-Artistic%202.0-0298c3.svg)](https://opensource.org/licenses/Artistic-2.0)


Raku package for sparse matrix algorithms:

- Implements (some of) the algorithms described (and spelled-out in FORTRAN) in the book 
"Sparse Matrix Technology" by S. Pissanetzky, [SP1].

- Provides convenient interface to accessing sparse matrix elements, rows, column, and sub-matrices. 

-----

## Motivation

Sparse Matrix Algebra (SMA) is a "must have" for many computational workflows. 
Here is a (non-exhaustive) list given in the order of _my_ preferences:

- Recommendation Systems (RS)
    - I make recommenders often during Exploratory Data Analysis (EDA).
    - For me, RS are "first order regression."
    - I also specialize in the making of RS.
    - I implemented a Raku recommender without SMA, 
      ["ML::StreamsBlendingRecommender"](https://github.com/antononcube/Raku-ML-StreamsBlendingRecommender), [AAp1], 
      but it is too slow for "serious" datasets.
      - Still useful; see [AAv1].
- Latent Semantic Analysis (LSA)
    - LSA is one my favorite Unsupervised Machine Learning (ML) workflows.
    - That means that this SMA package should have algorithms facilitating the programming of:
        - Singular Value Decomposition (SVD)
        - Non-Negative Matrix Factorization (NNMF)
- Graphs
    - There is a natural (representation) connection between sparse matrices and graphs.
    - Many graph algorithms can leverage (fast) SMA.
    - So far (2024-09-25) the algorithms in "Graph", [AAp2], do not use SMA and that is feature- and speed limiting.
- Optimization
    - For large scale optimization problems using SMA is a must.
        - Since their constraints are given with sparse matrices.
- Partial Differential Equations (PDE) solving

-----

## Usage examples

Here is a _simple_ sparse matrix in Compressed Sparse Row (CSR) format:

```perl6
use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

my $nrow = 5;
my $ncol = 8;
my $density = 0.2;
my $tol = 0.01;

my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
say $matrix1;
```

Here it is "pretty printed": 

```perl6
say "-" x 50;
$matrix1.print;
say "-" x 50;
```

Here `10` is added to all elements:

```perl6
my $matrix2 = $matrix1.add(10); 

say "-" x 50;
$matrix2.print;
say "-" x 50;
```

Here is the dot-product of the original matrix with its transpose:

```perl6
my $matrix3 = $matrix1.dot($matrix1.transpose); 

say "-" x 50;
$matrix3.print;
say "-" x 50;
```

-----

## Special features

Here are few features that other SMA packages typically do not provide.

### Named rows and columns

It is very convenient to have named rows and columns that are respected (or preserved)
in the typical SMA operations.

Here is an example:

```perl6
my $smat = Math::SparseMatrix.new(sparse-matrix => $matrix1, row-names => 'a' .. 'e', column-names => 'A' .. 'H'); 

say "-" x 50;
$smat.print;
say "-" x 50;
```

### Implicit value

-----

## Implementation notes

- The most important operation is matrix-vector multiplication.
    - It has to be as fast as possible.
    - There are two Dot Product implementations:
        - Direct
        - Symbolic-&-numeric
    - The direct one is 20-50% faster.
- It seems it a good idea to provide for some operations a _symbolic_ (or sparse matrix elements pattern) method.
    - For example:
        - `add-pattern` / `add`
        - `dot-pattern` / `dot-numeric`
- It is important to have access methods / operators.
    - All three are used in the accessor implementation: `AT-POS`, `postcircumfix:<[ ]>`, `postcircumfix:<[; ]>` .
- The core sparse matrix classes do not have named rows, columns, and dimensions. 
- That is done via suitable role. (Not implemented yest.) 

-----

## References

### Books

[SP1] Sergio Pissanetzky, Sparse Matrix Technology, Academic Pr (January 1, 1984), ISBN-10: 0125575807, ISBN-13: 978-0125575805.

### Packages

[AAp1] Anton Antonov,
[ML::StreamsBlendingRecommender Raku package](https://github.com/antononcube/Raku-ML-StreamsBlendingRecommender), 
(2021-2024),
[GitHub/antononcube](https://github.com/antononcube).

[AAp2] Anton Antonov,
[Graph Raku package](https://github.com/antononcube/Raku-Graph),
(2024),
[GitHub/antononcube](https://github.com/antononcube).

### Videos

[AAv1] Anton Antonov,
["TRC 2022 Implementation of ML algorithms in Raku"](https://youtu.be/efRHfjYebs4?si=J5P8pK1TgGSxdlmD&t=193),
(2022),
[YouTube/antononcube](https://www.youtube.com/@AAA4prediction).