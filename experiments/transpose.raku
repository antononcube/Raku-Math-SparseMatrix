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
#say "Non-zero values 1: ", $matrix1.rules;
#say "Non-zero values 1: ", $matrix1.values;
say "Non-zero values 1: ", $matrix1.values.elems;
$matrix1.print;
say "-" x 100;

my $result = $matrix1.transpose;
say "-" x 100;
say "Result :";
$result.print();

say "=" x 100;
note $result.verify():pairs;
note $result.row-ptr;
