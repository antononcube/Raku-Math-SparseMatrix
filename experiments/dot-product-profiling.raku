#!/usr/bin/env raku
use v6.d;

#use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

my $nrow = 50;
my $ncol = 8_000;
my $density = 0.2;
my $nnz = ($nrow * $ncol * $density).round;
my $tol = 0.0;
my $type = 'CSR';
my $n = 10;

my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):decorated;
my $matrix2 = generate-random-sparse-matrix(nrow => $ncol, ncol => $nrow, :$density, :$tol, :$type):decorated;

say '=' x 100;
say 'Raku implementation';
say '=' x 100;

say (:$matrix1);
say (:$matrix2);
my $tstart = now;
for ^$n {
    $matrix1.dot($matrix2)
}
my $tend = now;

say "Total time : { $tend - $tstart }";
say "Mean time  : { ($tend - $tstart) / $n }";
say '◻️';

say '=' x 100;
say 'Raku/C/NativeCall implementation';
say '=' x 100;

$matrix1 = $matrix1.to-adapted;
$matrix2 = $matrix2.to-adapted;

say (:$matrix1);
say (:$matrix2);

my $tstart2 = now;
for ^$n {
    $matrix1.dot($matrix2)
}
my $tend2 = now;

say "Total time : { $tend2 - $tstart2 }";
say "Mean time  : { ($tend2 - $tstart2) / $n }";
say '◻️';
