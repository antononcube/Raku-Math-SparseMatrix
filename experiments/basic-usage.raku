#!/usr/bin/env raku
use v6.d;

use lib <. lib>;
use Math::SparseMatrix::CSR;

# Example tuples: (row, col, value)
my @entries = (
(0, 0, 1),
(0, 2, 2),
(1, 1, 3),
(2, 0, 4),
(2, 2, 5),
(3, 2, 7)
);

my @rules = @entries.map({ Pair.new($_.head(2).List, $_.tail) });

say (@rules.raku);

# Create a row-major sparse matrix
my $matrix = Math::SparseMatrix::CSR.new(:@rules);

say "Matrix: ", $matrix.raku;

say "Matrix rules: ", $matrix.rules.raku;
say '-' x 100;
$matrix.print;
say '-' x 100;


# First two rows
say '-' x 100;
$matrix.head(2).print;
say '-' x 100;

# Multiply with a dense vector
my @vector = [1, 2, 3];
my @result = $matrix.dot(@vector);
say "Matrix * Vector: ", @result;

# Transpose the matrix
my $transposed = $matrix.transpose();
say "Transposed Matrix: ", $transposed.raku;
$transposed.print;


# Multiply two matrices
my $product = $matrix.dot($transposed);
say "Product Matrix: ", $product.raku;

say "Product Matrix pretty: ", $product.print;