#!/usr/bin/env raku
use v6.d;

use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

use Data::Importers;
use Data::Summarizers;
use Data::Reshapers;

my $nrow = 5;
my $ncol = 8;
my $density = 0.25;
my $tol = 0.001;
my $type = 'CSR';

# Here we create a decorated sparse matrix
my $smat = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type):decorated;

# Here that matrix is given row-names and column-names:
$smat =
        Math::SparseMatrix.new(
                core-matrix => $smat,
                row-names => ('a' .. 'e').pairs.Hash.invert.Hash.deepmap(*.Int),
                column-names => ('A' .. 'H').pairs.Hash.invert.Hash.deepmap(*.Int)
        );

say "\n" x 2;
# Here it is a "pretty print" of it
$smat.print;

say "\n" x 2;
# Here is matrix's transpose
$smat.transpose.print;

say "\n" x 2;
# Here is row slice specified by using row names
my $smat2 = $smat.row-slice(<a c d>);
$smat2.print;

say "\n" x 2;
# Here is row slice specified by using row names
my $smat2a = $smat.row-slice(<d c a>);
$smat2a.print;

say "\n" x 2;
# Here is a column sub-matrix specified with a column name
my $smat3 = $smat.column-at('B');
$smat3.print;

say "\n" x 2;
# Alternative specification using postcircumfix:<[; ]>
my $smat4 = $smat[*;'B'];
$smat4.print;

say "\n" x 2;

my $smat5 = $smat['a'..'c';'B'..'D'];
$smat5.print;

say "\n" x 2;
my $smat6 = $smat5.multiply(100);
$smat6.print;

say "\n" x 2;
my $smat7 = $smat.dot(10 xx $smat.columns-count);
$smat7.print;

say "\n" x 2;
my $smat8 = $smat.dot($smat.transpose);
$smat8.round(0.02).print;

say $smat8;
say $smat8.Hash;

say $smat8.to-wl;