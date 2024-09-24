#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use SparseMatrixCSR;

sub generate-random-sparse-matrix(UInt:D $nrow, UInt:D $ncol, Numeric:D $density) returns SparseMatrixCSR {
    my $n = ($nrow * $ncol * $density).Int;

    my @rules = ((^$nrow).pick($n) Z (^$ncol).pick($n)) Z=> (round(rand, 0.2) xx $n);

    return SparseMatrixCSR.new(
            :@rules,
            :$nrow,
            :$ncol
            );
}


my $nrow = 10_000;
my $ncol = 400;
my $density = 0.02;

say "-" x 100;
say "Matrix 1:";
my $matrix1 = generate-random-sparse-matrix($nrow, $ncol, $density);
#say "Non-zero values 1: ", $matrix1.rules;
#say "Non-zero values 1: ", $matrix1.values;
say "Non-zero values 1: ", $matrix1.values.elems;
#$matrix1.print;
say "-" x 100;

say "-" x 100;
say "Matrix 2:";
#my $matrix2 = generate-random-sparse-matrix($ncol, 2, 0.9);
my $matrix2 = $matrix1.transpose;
#say "Non-zero values 2: ", $matrix2.values;
say "Non-zero values 2: ", $matrix2.values.elems;
#$matrix2.print;
say "-" x 100;



my $tstart = now;
my $result = $matrix1.dot($matrix2);
my $tend = now;

say "Multiplication time: {$tend - $tstart} seconds.";