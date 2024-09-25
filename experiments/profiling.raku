#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::CSR;


my $nrow = 10_000;
my $ncol = 400;
my $density = 0.001;
my $tol = 0.2;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
#say "Non-zero values 1: ", $matrix1.rules;
#say "Non-zero values 1: ", $matrix1.values;
say "Non-zero values 1: ", $matrix1.values.elems;
#$matrix1.print;
say "-" x 100;

say "-" x 100;
say "Matrix 2:";
#my $matrix2 = generate-random-sparse-matrix($ncol, 2, 0.9);

my $tstart = now;
my $matrix2 = $matrix1.transpose;
my $tend = now;
say "Transpose time: {$tend - $tstart} seconds.";

#say "Non-zero values 2: ", $matrix2.values;
say "Non-zero values 2: ", $matrix2.values.elems;
#$matrix2.print;
say "-" x 100;



$tstart = now;
my $result = $matrix1.dot($matrix2);
$tend = now;

say "Multiplication time dot        : {$tend - $tstart} seconds.";

$tstart = now;
$result = $matrix1.dot-numeric($matrix2);
$tend = now;

say "Multiplication time dot-numeric: {$tend - $tstart} seconds.";