#!/usr/bin/env raku
use v6.d;

# use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

my $nrow = 1_000;
my $ncol = 10_000;
my $density = 0.002;
my $tol = 0.0;
my $type = 'CSR';

#----------------------------------------------------------------------------------------------------
say "-" x 100;
say "Matrix 1:";
my $tstart = now;
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):!decorated;
my $tend = now;
say (:$matrix1);
say "Creation time: {$tend - $tstart} seconds.";
say "Non-zero values 1: ", $matrix1.explicit-length;
#----------------------------------------------------------------------------------------------------

say "-" x 100;

#spurt 'test-sparse-matrix.wl', $matrix1.to-wl;

#----------------------------------------------------------------------------------------------------
say "-" x 100;
say "Matrix 2:";
#my $matrix2 = generate-random-sparse-matrix($ncol, 2, 0.9);

$tstart = now;
my $matrix2 = $matrix1.transpose;
$tend = now;
say (:$matrix2);
say "Transpose time: {$tend - $tstart} seconds.";

say "Non-zero values 2: ", $matrix2.explicit-length;

#----------------------------------------------------------------------------------------------------
say "-" x 100;

#$tstart = now;
#my $result = $matrix1.dot($matrix2);
#$tend = now;
#say "Multiplication time dot        : {$tend - $tstart} seconds.";

$tstart = now;
my $result2 = $matrix1.dot-numeric($matrix2);
$tend = now;
say "Multiplication time dot-numeric: {$tend - $tstart} seconds.";

#----------------------------------------------------------------------------------------------------
say "-" x 100;

$tstart = now;
my $result3 = $matrix1.values.rotor($nrow, :partial).map({ sum($_ <<*>> $_) });
#note (:$result3);
$tend = now;

say "Simple dot products time       : {$tend - $tstart} seconds.";

#----------------------------------------------------------------------------------------------------
say "-" x 100;

my $matrix3 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):decorated;

$tstart = now;
my @rowSums = $matrix3.row-sums;
$tend = now;

say "Generic row sums time : {$tend - $tstart} seconds.";

$tstart = now;
@rowSums = $matrix3.core-matrix.row-sums;
$tend = now;

say "CSR row sums time     : {$tend - $tstart} seconds.";


my $matrix3a = $matrix3.clone;
$matrix3a.to-adapted;

$tstart = now;
@rowSums = $matrix3a.core-matrix.row-sums;
$tend = now;

say "Adapted row sums time : {$tend - $tstart} seconds.";

#----------------------------------------------------------------------------------------------------
say "-" x 100;

$tstart = now;
my $matrix3to10 = $matrix3.top-k-elements-matrix(10);
$tend = now;

say (:$matrix3to10);
say "CSR top-K elements matrix time : {$tend - $tstart} seconds.";

$tstart = now;
my $matrix3a10 = $matrix3a.top-k-elements-matrix(10);
$tend = now;

say (:$matrix3a10);
say "Adapted top-K elements matrix time : {$tend - $tstart} seconds.";

say 'explicit-length    : ', $matrix3a10.core-matrix.explicit-length;
say 'explicit-positions : ', $matrix3a10.core-matrix.explicit-positions;
say 'explicit-values    : ', $matrix3a10.core-matrix.explicit-values;