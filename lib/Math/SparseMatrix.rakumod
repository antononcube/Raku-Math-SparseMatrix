unit module Math::SparseMatrix;

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

#=====================================================================
#proto sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {*}
multi sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {
    return $mat.row-slice(@indexes);
}

multi sub postcircumfix:<[; ]>(Math::SparseMatrix::CSR:D $mat, @indexes) is export {
    return do given (@indexes[0], @indexes[1]) {
        when $_.head ~~ Int && $_.tail ~~ Int {
            $mat.value-at(@indexes[0], @indexes[1]);
        }
        when $_.head ~~ Range && $_.tail ~~ Int {
            $mat.row-slice($_.head).column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Int {
            $mat.column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Range {
            $mat.transpose.row-slice($_.tail).transpose
        }
        when $_.head ~~ Range && $_.tail ~~ Range {
            $mat.row-slice($_.head).transpose.row-slice($_.tail).transpose
        }
        default {
            die "Cannot process the given range"
        }
    }
}