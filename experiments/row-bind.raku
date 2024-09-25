#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::CSR;


my $nrow = 5;
my $ncol = 8;
my $density = 0.25;
my $tol = 0.01;

say "=" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
say $matrix1;
say "-" x 100;
$matrix1.print;
say "-" x 100;

say "=" x 100;
say "Matrix 2:";
my $matrix2 = generate-random-sparse-matrix($nrow, $ncol+2, :$density, :$tol);
say $matrix2;
say "-" x 100;
$matrix2.print;
say "-" x 100;


#my $result = $matrix1.row-bind($matrix2);
my $result = $matrix1.column-bind($matrix2);
say "=" x 100;
say "Result :";
say "-" x 100;
$result.print();

say "=" x 100;
note $matrix1.transpose.verify():pairs;
