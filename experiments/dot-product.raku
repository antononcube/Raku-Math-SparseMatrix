#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::Utilities;

my $nrow = 5;
my $ncol = 8;
my $density = 0.2;
my $tol = 0.01;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol);
say $matrix1;
say "-" x 100;
$matrix1.print;
say "-" x 100;

say $matrix1.wl();
say "-" x 100;

say "-" x 100;
say "Matrix 2:";
my $matrix2 = generate-random-sparse-matrix($ncol, $nrow, :$density, :$tol);
say $matrix2;
say "-" x 100;
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
