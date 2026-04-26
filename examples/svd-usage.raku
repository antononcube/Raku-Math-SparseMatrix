#!/usr/bin/env raku
use v6.d;

# Corresponding Wolfram Language code
# m = {{0.5, 0.5, 0.5, 0, 0}, {1., 2., 4., 8., 0}, {1.5, 4.5, 13.5, 40.5, 121.5}, {0, 8., 32., 128., 512.}, {0, 0, 62.5, 312.5, 1562.5}};
# {u, s, v} = SingularValueDecomposition[m, 3];
# MatrixForm /@ Map[Round[#, 0.0001] &, {u, s, v}, {-1}]

use Math::SparseMatrix;

my @m = [0.5,0.5,0.5,0,0],[1.0,2.0,4.0,8.0,0],[1.5,4.5,13.5,40.5,121.5],[0,8.0,32.0,128.0,512.0],[0,0,62.5,312.5,1562.5],[0,8.0,4.0,2.0,1.0];

#my $matrix = Math::SparseMatrix::NativeAdapter.new(dense-matrix => @m, row-names => 'ABCDE'.comb, column-names => 'abcde'.comb);
my $matrix = Math::SparseMatrix.new(dense-matrix => @m, row-names => 'ABCDEF'.comb, column-names => 'abcde'.comb).to-adapted;

#====================================================================================================
say '=' x 100;
say (:$matrix);

my ($u, $s, $v) = $matrix.svd(3);

my $tol = 0.0001;

#====================================================================================================
say '=' x 100;

say (:$u);
$u.round($tol).print;

#====================================================================================================
say '=' x 100;

say (:$s);
$s.round($tol).print;

#=====================================================================================================
say '=' x 100;

say (:$u);
$v.round($tol).print;