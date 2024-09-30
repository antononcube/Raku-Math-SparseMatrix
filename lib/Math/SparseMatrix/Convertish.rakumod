use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::DOK;

role Math::SparseMatrix::Convertish {

    #=================================================================
    # To CSR
    #=================================================================
    #| To CSR sparse matrix.
    proto method to-csr(Math::SparseMatrix::Abstract:D $m) {*}

    multi method to-csr(Math::SparseMatrix::CSR:D $m) {
        return $m;
    }

    multi method to-csr(Math::SparseMatrix::DOK:D $m) {
        my @values;
        my @col-index;
        my @row-ptr = 0;
        my $nnz = 0;

        for ^$m.nrow -> $row {
            for %.adjacency-list{$row}.kv -> $col, $value {
                @values.push($value);
                @col-index.push($col);
                $nnz++;
            }
            @row-ptr.push($nnz);
        }

        return Math::SparseMatrix::CSR.new(
                :@values,
                :@col-index,
                :@row-ptr,
                nrow => $m.nrow,
                ncol => $m.ncol,
                implicit-value => $m.implicit-value
                );
    }

    #=================================================================
    # To DOK
    #=================================================================
    #| To DOK sparse matrix.
    proto method to-dok(Math::SparseMatrix::Abstract:D $m) {*}

    multi method to-dok(Math::SparseMatrix::DOK:D $m) {
        return $m;
    }

    multi method to-dok(Math::SparseMatrix::CSR:D $m) {
        my %adjacency-map;
        for ^$m.nrow -> $row {
            my $start = $m.row-ptr[$row];
            my $end = $m.row-ptr[$row + 1];
            for $start ..^ $end -> $idx {
                %adjacency-map{$row}{$m.col-index[$idx]} = $m.values[$idx];
            }
        }
        return Math::SparseMatrix::DOK.new(
                :%adjacency-map,
                :nrow($m.nrow),
                :ncol($m.ncol),
                :implicit-value($m.implicit-value)
                );
    }
}