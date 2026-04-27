#!/usr/bin/env raku
use v6.d;

#use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

my $nrow = 10_000;
my $ncol = 1200;
my $density = 0.002;
my $tol = 0.0;
my $type = 'CSR';
my $k = 100;
my $tolerance = 1e-4;

say "-" x 100;
say "Matrix:";
my $tstart = now;
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type).to-adapted;
my $tend = now;
say (:$matrix1);
say "Creation time: { $tend - $tstart } seconds.";
say "Non-zero values: ", $matrix1.explicit-length;
say "Fill in: ", $matrix1.explicit-length / $matrix1.rows-count / $matrix1.columns-count;
say "-" x 100;

$tstart = now;
my ($u, $s, $v) = $matrix1.svd($k, :$tolerance);
$tend = now;
say "SVD time: { $tend - $tstart } seconds.";


.say for |(:$u, :$s, :$v);
