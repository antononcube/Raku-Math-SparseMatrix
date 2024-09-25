#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::CSR;


my $nrow = 5;
my $ncol = 8;
my $density = 0.15;
my $tol = 0.01;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
say $matrix1;
say "-" x 100;
$matrix1.print;
say "-" x 100;

say "=" x 100;
say "Result :";
say "-" x 100;
my $result = $matrix1.transpose;
$result.print();

say "=" x 100;
note $result.verify():pairs;
note $result.row-ptr;
