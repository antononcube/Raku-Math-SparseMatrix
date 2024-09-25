unit module Math::SparseMatrix;

use Math::SparseMatrix::CSR;

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

#multi sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, **@indexes) {
#    die 'The indexes are expected to be non-negative integers.'
#    unless @indexes.all ~~ Int:D && min(@indexes) ≥ 0;
#
#    my @mats = @indexes.map({ $mat.row-at($_) });
#    my $res = @mats.head;
#    for @mats.tail(*-1) -> $m {
#        $res = $res.row-bind($m)
#    }
#    return $res;
#}