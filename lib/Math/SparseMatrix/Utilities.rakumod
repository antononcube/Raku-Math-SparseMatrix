unit module Math::SparseMatrix::Utilities;

use Math::SparseMatrix;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::DOK;

#=====================================================================
# Random sparse matrix
#=====================================================================
proto sub generate-random-sparse-matrix(|) is export {*}

multi sub generate-random-sparse-matrix(
        UInt:D $nrow,
        UInt:D $ncol,
        Numeric:D :$density = 0.01,
        Numeric:D :$tol = 0.001,
        :$type = Whatever,
        Bool:D :d(:$decorated) = True) {
    return generate-random-sparse-matrix(:$nrow, :$ncol, :$density, :$tol, :$type, :$decorated);
}

multi sub generate-random-sparse-matrix(
        UInt:D :$nrow,
        UInt:D :$ncol,
        Numeric:D :$density = 0.01,
        Numeric:D :$tol = 0.001,
        :$type is copy = Whatever,
        Bool:D :d(:$decorated) = True) {

    my $n = ($nrow * $ncol * $density).Int;
    my @rules = do if $tol {
        (($nrow.rand.floor, $ncol.rand.floor) => round(rand, $tol)) xx $n
    } else {
        (($nrow.rand.floor, $ncol.rand.floor) => rand) xx $n
    }

    my $core-matrix = do given $type {
        when Whatever {
            Math::SparseMatrix::CSR.new(:@rules, :$nrow, :$ncol)
        }
        when ($_ ~~ Str:D) && $_.lc eq 'csr' {
            Math::SparseMatrix::CSR.new(:@rules, :$nrow, :$ncol)
        }
        when ($_ ~~ Str:D) && $_.lc eq 'dok' {
            Math::SparseMatrix::DOK.new(rules => @rules.Hash, :$nrow, :$ncol)
        }
        default {
            die 'The argument $type is expected to be one of "CSR", "DOK", or Whatever.'
        }
    }
    return $decorated ?? Math::SparseMatrix.new($core-matrix) !! $core-matrix;
}

#=====================================================================
# Element-wise operation
#=====================================================================
sub elementwise-operation(@matrix1, @matrix2, &op) is export {
    die "Incompatible matrices" unless @matrix1.elems == @matrix2.elems && @matrix1[0].elems == @matrix2[0].elems;
    my @result;
    for @matrix1.kv -> $i, @row {
        for @row.kv -> $j, $elem {
            @result[$i][$j] = op($elem, @matrix2[$i][$j]);
        }
    }
    return @result;
}

#=====================================================================
# Dot Product for dense matrices (arrays)
#=====================================================================
sub dense-dot-product(@matrix1, @matrix2) is export {
    die "Incompatible matrices" unless @matrix1[0].elems == @matrix2.elems;
    my @result;
    for @matrix1.kv -> $i, @row {
        for @matrix2[0].keys -> $j {
            @result[$i][$j] = [+] @row Z* @matrix2.map({ $_[$j] });
        }
    }
    return @result;
}