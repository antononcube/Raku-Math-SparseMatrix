#!/usr/bin/env raku
use v6.d;

use Math::SparseMatrix;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::Utilities;


my $nrow = 5;
my $ncol = 8;
my $density = 0.35;
my $tol = 0.01;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
say $matrix1;
say "-" x 100;
$matrix1.print;
say "-" x 100;

say "=" x 100;
say "Result [2,3]:";
say "-" x 100;
my $result = $matrix1[2;3];
say (:$result);

say "=" x 100;
say "Result [^3,2]:";
say "-" x 100;
$result = $matrix1[^3;2];
$result.print();

say "=" x 100;
say "Result [*;1..2]:";
say "-" x 100;
$result = $matrix1[*;1..2];
$result.print();

say "=" x 100;
say "Result [^3;2..4]:";
say "-" x 100;
$result = $matrix1[^3;2..4];
$result.print();

