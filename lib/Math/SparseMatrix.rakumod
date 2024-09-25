unit module Math::SparseMatrix;

use Math::SparseMatrix::CSR;

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