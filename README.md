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
```
# Math::SparseMatrix(:specified-elements(8), :dimensions((5, 8)), :density(0.2))
```

Here it is "pretty printed": 

```perl6
$matrix1.print;
```
```
# –––––––––––––––––––––––––––––––––––––––––––
#     0    1    2    3    4    5    6    7   
# ––┼––––––––––––––––––––––––––––––––––––––––
# 0 │ .    .    0.58 .    0.6  .    .    0.19
# 1 │ .    .    .    .    .    .    .    .   
# 2 │ .    .    .    .    .    .    0.31 .   
# 3 │ .    .    .    0.79 .    .    .    .   
# 4 │ 0.03 .    0.73 .    .    .    .    0.9
```

Here `10` is multiplied with all elements:

```perl6
my $matrix2 = $matrix1.multiply(10); 

$matrix2.print;
```
```
# –––––––––––––––––––––––––––––––––––
#     0   1   2   3   4   5   6   7  
# ––┼––––––––––––––––––––––––––––––––
# 0 │ .   .   5.8 .   6   .   .   1.9
# 1 │ .   .   .   .   .   .   .   .  
# 2 │ .   .   .   .   .   .   3.1 .  
# 3 │ .   .   .   7.9 .   .   .   .  
# 4 │ 0.3 .   7.3 .   .   .   .   9
```

Here is the dot-product of the original matrix with its transpose:

```perl6
my $matrix3 = $matrix1.dot($matrix1.transpose); 

$matrix3.print;
```
```
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#     0                  1                  2                  3                  4                 
# ––┼–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# 0 │ 0.7324999999999999 .                  .                  .                  0.5943999999999999
# 1 │ .                  .                  .                  .                  .                 
# 2 │ .                  .                  0.0961             .                  .                 
# 3 │ .                  .                  .                  0.6241000000000001 .                 
# 4 │ 0.5943999999999999 .                  .                  .                  1.3437999999999999
```

-----

## Special features

Here are few features that other SMA packages typically do not provide.

### Named rows and columns

It is very convenient to have named rows and columns that are respected (or preserved)
in the typical SMA operations.

Here is an example:

```perl6
my $smat = Math::SparseMatrix.new($matrix1, row-names => 'a' .. 'e', column-names => 'A' .. 'H'); 

$smat.print;
```
```
# –––––––––––––––––––––––––––––––––––––––––––
#     A    B    C    D    E    F    G    H   
# ––┼––––––––––––––––––––––––––––––––––––––––
# a │ .    .    0.58 .    0.6  .    .    0.19
# b │ .    .    .    .    .    .    .    .   
# c │ .    .    .    .    .    .    0.31 .   
# d │ .    .    .    0.79 .    .    .    .   
# e │ 0.03 .    0.73 .    .    .    .    0.9
```


Here is the dot-product of that matrix with its transpose:

```perl6
my $smat2 = $smat.dot($smat.transpose); 

$smat2.print;
```
```
# ––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
#     a                  b                  c                  d                  e                 
# ––┼–––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––––
# a │ 0.7324999999999999 .                  .                  .                  0.5943999999999999
# b │ .                  .                  .                  .                  .                 
# c │ .                  .                  0.0961             .                  .                 
# d │ .                  .                  .                  0.6241000000000001 .                 
# e │ 0.5943999999999999 .                  .                  .                  1.3437999999999999
```

### Implicit value

The sparse matrices can have an _implicit value_ that is different from 0.
For example, adding a number to a sparse matrix produces another sparse matrix 
with different implicit value:

```perl6
my $matrix3 = $matrix1.add(10);
```
```
# Math::SparseMatrix(:specified-elements(8), :dimensions((5, 8)), :density(0.2))
```

```perl6
$matrix3.implicit-value
```
```
# 10
```

Here is the pretty print:

```perl6
$matrix3.print(:iv)
```
```
# –––––––––––––––––––––––––––––––––––––––––––––––––––
#     0     1     2     3     4     5     6     7    
# ––┼––––––––––––––––––––––––––––––––––––––––––––––––
# 0 │ .     .     10.58 .     10.6  .     .     10.19
# 1 │ .     .     .     .     .     .     .     .    
# 2 │ .     .     .     .     .     .     10.31 .    
# 3 │ .     .     .     10.79 .     .     .     .    
# 4 │ 10.03 .     10.73 .     .     .     .     10.9
```

**Remark:** Currently, the implicit values are ignored in `dot`.

-----

## Design

### General

- There should be a "main" class, `Math::SpareMatrix` that:
  - Provides the SMA functionalities
  - Delegates to concrete sparse matrix classes that are based on different representation formats
  - Can have named rows, columns, and dimensions
  - Gives access to sparse matrix elements, rows, columns, and sub-matrices
- The default or "main" core sparse matrix class should use Compressed Sparse Row (CSR) format.
- Also, a class using Dictionary Of Keys (DOK) format should be provided.
- The core sparse matrix classes do not have named rows, columns, and dimensions.
- Ideally, a class using `NativeCall` should be implemented at some point.
  - It looks like this is "a must", since the CSR and DOK classes are fairly slow.
  - Both "plain C" and macOS [Accelerate](https://developer.apple.com/accelerate/) implementations should be made.
- The _most important operation_ is Matrix-Vector Dot Product.
  - The current design is to use one-row or one-column matrices for the vectors.
  - Dense vectors are (of course) also supported

### Object-Oriented Programming (OOP) architecture

- The OOP [Decorator Design Pattern](https://en.wikipedia.org/wiki/Decorator_pattern) is used to organize the SMA functionalities.
- In that pattern:
  - The _Component_ is played by the class [`Math::SparseMatrix::Abstract`](./lib/Math/SparseMatrix/Abstract.rakumod).
  - The _ConcreteComponent_ is played by the classes:
    - [`Math::SparseMatrix::CSR`](./lib/Math/SparseMatrix/CSR.rakumod)
    - [`Math::SparseMatrix::DOK`](./lib/Math/SparseMatrix/DOK.rakumod)
  - The concrete component classes provide the core SMA operations.
  - The _Decorator_ is played by [`Math::SparseMatrix`](./lib/Math/SparseMatrix.rakumod).
    - That is a "top level", interface class.
    - Allows access using named rows and columns.
    - "Hides" the actual component class used.

Here is a corresponding diagram:

```mermaid
classDiagram
    class Abstract["Math::SparseMatrix::Abstract"] {
        <<abstract>>
        +value-at()
        +row-at()
        +column-at()
        +row-slice()
        +AT-POS()
        +print()
        +transpose()
        +add()
        +multiply()
        +dot()
    }
    
    class CSR["Math::SparseMatrix::CSR"] {
        @row-ptr
        @col-index
        @values
        nrow
        ncol
        implicit-value
    }
    
    class DOK["Math::SparseMatrix::DOK"] {
        %adjacency-map
        nrow
        ncol
        implicit-value
    }
    
    class SparseMatrix["Math::SparseMatrix"] {
        Abstract core-matrix
        +AT-POS()
        +print()
        +transpose()
        +add()
        +multiply()
        +dot()
    }
    
    CSR --> Abstract : implements
    DOK --> Abstract : implements
    SparseMatrix --> Abstract : Hides actual component class
    SparseMatrix *--> Abstract
```

### Implementation details

- Again, the most important operation is Matrix-Vector Dot Product.
    - It has to be as fast as possible.
    - There are two Dot Product implementations for CSR:
        - Direct
        - Symbolic-&-numeric
    - (Currently) the direct one is 20-50% faster.
- It seems it is a good idea to provide for some operations _symbolic_ (or sparse matrix elements pattern) methods.
    - For example:
        - `add-pattern` / `add`
        - `dot-pattern` / `dot-numeric`
- It is important to have access methods / operators.
    - All three are used in the accessor implementation: `AT-POS`, `postcircumfix:<[ ]>`, `postcircumfix:<[; ]>` .


-----

## Performance

- Performance of CSR and DOK sparse matrices is not good: between 40 to 150 times slower than Wolfram Language.
    - (Using the same matrices, of course.)
- It is somewhat surprising that DOK is faster than CSR.

-----

## Acknowledgements

Thanks to [@lizmat](https://github.com/lizmat) and [@timo](https://github.com/timo) for helping figuring out the proper use of `postcircumfix:<[]>` and `postcircumfix:<[; ]>`
in order to have the named rows and columns functionalities.

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