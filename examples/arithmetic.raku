#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix :ALL;
use Math::SparseMatrix::Utilities;

my $nrow = 5;
my $ncol = 8;
my $density = 0.25;
my $tol = 0.001;
my $type = 'CSR';

# Here we create a decorated sparse matrix
my $smat1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):decorated;
my $smat2 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):decorated;

# Here that matrix is given row-names and column-names:
$smat1 =
        Math::SparseMatrix.new(
                core-matrix => $smat1,
                row-names => 'a' .. 'e',
                column-names => 'A' .. 'H'
        );

$smat2 =
        Math::SparseMatrix.new(
                core-matrix => $smat2,
                row-names => 'a' .. 'e',
                column-names => 'A' .. 'H'
                );

say "\n" x 2;
# pretty print the matrices
$smat1.print;
say '-' x 100;
$smat2.print;

say "\n" x 2;
say 'Sum of a matrix and a scalar:';
my $scalar = 10;
($smat1 + $scalar).print(:iv);

say "implicit value : ", ($smat1 + $scalar).implicit-value;

say "\n" x 2;
say 'Sum of the matrices:';
($smat1 + $smat2).print;

say "\n" x 2;
say 'Multiplication of the matrices:';
($smat1 * $smat2).print;

say "\n" x 2;
say 'Multiplication by scalar and sum:';
($smat1 + 10000 * $smat2).print;


say "\n" x 2;
say 'Dot product:';
($smat1 • $smat2.transpose + 100).print(:iv);
