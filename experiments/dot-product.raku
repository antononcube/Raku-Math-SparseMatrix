#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::CSR;


my $nrow = 5;
my $ncol = 8;
my $density = 0.2;
my $tol = 0.01;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
#say "Non-zero values 1: ", $matrix1.rules;
#say "Non-zero values 1: ", $matrix1.values;
say "Non-zero values 1: ", $matrix1.values.elems;
$matrix1.print;
say "-" x 100;

say "-" x 100;
say "Matrix 2:";
my $matrix2 = generate-random-sparse-matrix($ncol, $nrow, :$density, :$tol);
#say "Non-zero values 2: ", $matrix2.values;
say "Non-zero values 2: ", $matrix2.values.elems;
$matrix2.print;
say "-" x 100;


my $pattern = $matrix1.dot-pattern($matrix2);

say "-" x 100;
say "Pattern :";
$pattern.print();

my $result = $matrix1.dot($matrix2);
say "-" x 100;
say "Result :";
$result.print();

my $result2 = $matrix1.dot-numeric($matrix2);
say "-" x 100;
say "Result2 :";
$result2.print();
