unit module Math::SparseMatrix::Utilities;

use Math::SparseMatrix::CSR;

#=====================================================================
proto sub generate-random-sparse-matrix(|) is export {*}

multi sub generate-random-sparse-matrix(
        UInt:D $nrow,
        UInt:D $ncol,
        Numeric:D :d(:$density) = 0.01,
        Numeric:D :$tol = 0.001) {
    return generate-random-sparse-matrix(:$nrow, :$ncol, :$density, :$tol);
}

multi sub generate-random-sparse-matrix(
        UInt:D :$nrow,
        UInt:D :$ncol,
        Numeric:D :d(:$density) = 0.01,
        Numeric:D :$tol = 0.001)
{
    my $n = ($nrow * $ncol * $density).Int;
    my @rules = ((^$nrow) X (^$ncol)).pick($n) Z=> (round(rand, $tol) xx $n);
    return Math::SparseMatrix::CSR.new(:@rules, :$nrow, :$ncol);
}
