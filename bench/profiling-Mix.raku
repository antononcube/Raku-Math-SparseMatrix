#!/usr/bin/env raku
use v6.d;

#use lib <. lib>;
use Math::SparseMatrix;
use Math::SparseMatrix::Utilities;

my $nrow = 1000;
my $ncol = 10_000;
my $density = 0.002;
my $tol = 0.0;
my $type = 'DOK';

say "-" x 100;
say "Matrix 1:";
my $tstart = now;
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, :$density, :$tol, :$type);
my $tend = now;
say (:$matrix1);
say "Creation time: { $tend - $tstart } seconds.";
say "Non-zero values 1: ", $matrix1.explicit-length;
say "-" x 100;

#spurt 'test-sparse-matrix.wl', $matrix1.wl;

say "-" x 100;
say "Matrix 2:";
#my $matrix2 = generate-random-sparse-matrix($ncol, 2, 0.9);

$tstart = now;
my $matrix2 = $matrix1.transpose;
$tend = now;
say (:$matrix2);
say "Transpose time: { $tend - $tstart } seconds.";

say "Non-zero values 2: ", $matrix2.explicit-length;
say "-" x 100;


# Convert
my %adj-map1 = $matrix1.adjacency-list.nodemap(*.Mix);
say (:%adj-map1);

my %adj-map2 = $matrix2.transpose.adjacency-list.nodemap(*.Mix);
say (:%adj-map2);

sub mix-dot-product(%adj-map1, %adj-map2) {
    my %result;
    for %adj-map1.kv -> $k1, %row1 {
        my @rowInds1 = %row1.keys;
        for %adj-map2.kv -> $k2, %row2 {
            %result{$k1}{$k2} = (%row2{@rowInds1} <<*>> %row1.values).sum
        }
    }
}

sub dot-product(Math::SparseMatrix::DOK:D $matrix1, Math::SparseMatrix::DOK:D $matrix2) {
    my %adj-map1 = $matrix1.adjacency-list.nodemap(*.Mix);
    my %adj-map2 = $matrix2.transpose.adjacency-list.nodemap(*.Mix);
    my %result;
    for %adj-map1.kv -> $k1, %row1 {
        my @rowInds1 = %row1.keys;
        for %adj-map2.kv -> $k2, %row2 {
            %result{$k1}{$k2} = (%row2{@rowInds1} <<*>> %row1.values).sum
        }
    }
    Math::SparseMatrix::DOK.new(adjacency-map => %result, nrow => $matrix1.nrow, ncol => $matrix2.ncol);
}


$tstart = now;
my $result = mix-dot-product(%adj-map1, %adj-map2);
$tend = now;

say "Multiplication time mix-dot-product : { $tend - $tstart } seconds.";


$tstart = now;
my $result2 = dot-product($matrix1, $matrix2);
$tend = now;

say "Multiplication time mix-dot-product : { $tend - $tstart } seconds.";

my $result3 = $matrix1.dot($matrix2);

say $result3.eqv($result2);

$result2.print;
say '-' x 120;
$result3.print;
