use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::CSC;
use Math::SparseMatrix::DOK;
use Math::SparseMatrix::NativeAdapter;

role Math::SparseMatrix::Convertish {

    #=================================================================
    # To CSR
    #=================================================================
    #| To CSR sparse matrix.
    proto method to-csr(Math::SparseMatrix::Abstract:D $m) {*}

    multi method to-csr(Math::SparseMatrix::CSR:D $m) {
        return $m;
    }

    multi method to-csr(Math::SparseMatrix::CSC:D $m) {
        return self.to-csr(self.to-dok($m));
    }

    multi method to-csr(Math::SparseMatrix::DOK:D $m) {
        my @values;
        my @col-index;
        my @row-ptr = 0;
        my $nnz = 0;

        for ^$m.nrow -> $row {
            for $m.adjacency-map{$row}.kv -> $col, $value {
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
    # To CSC
    #=================-===============================================
    #| To CSC sparse matrix.
    proto method to-csc(Math::SparseMatrix::Abstract:D $m) {*}

    multi method to-csc(Math::SparseMatrix::CSC:D $m) {
        return $m;
    }

    multi method to-csc(Math::SparseMatrix::CSR:D $m) {
        my $mt = $m.transpose;
        return
                Math::SparseMatrix::CSC.new(
                        values => $mt.values,
                        col-ptr => $mt.row-ptr,
                        row-index => $mt.col-index,
                        nrow => $mt.ncol,
                        ncol => $mt.nrow,
                        implicit-value => $mt.implicit-value
                )
    }

    multi method to-csc(Math::SparseMatrix::DOK:D $m) {
       return self.to-csc(self.to-csr($m))
    }

    #=================================================================
    # To DOK
    #=================================================================
    #| To DOK sparse matrix.
    proto method to-dok(|) {*}

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
                nrow => $m.nrow,
                ncol => $m.ncol,
                implicit-value => $m.implicit-value
                );
    }

    multi method to-dok(Math::SparseMatrix::CSC:D $m) {
        my %adjacency-map;
        for ^$m.ncol -> $col {
            my $start = $m.col-ptr[$col];
            my $end = $m.col-ptr[$col + 1];
            for $start ..^ $end -> $idx {
                %adjacency-map{$m.row-index[$idx]}{$col} = $m.values[$idx];
            }
        }
        return Math::SparseMatrix::DOK.new(
                :%adjacency-map,
                nrow => $m.nrow,
                ncol => $m.ncol,
                implicit-value => $m.implicit-value
                );
    }

    #=================================================================
    # To NativeAdapter
    #=================================================================
    proto method to-adapted(| --> Math::SparseMatrix::NativeAdapter:D) {*}

    multi method to-adapted(Math::SparseMatrix::NativeAdapter:D $m) {
        return $m;
    }
    multi method to-adapted(Math::SparseMatrix::CSR:D $m --> Math::SparseMatrix::NativeAdapter:D) {
        return Math::SparseMatrix::NativeAdapter.new(
                values => $m.values,
                row-ptr => $m.row-ptr,
                col-index => $m.col-index,
                nrow => $m.nrow,
                ncol => $m.ncol,
                implicit-value => $m.implicit-value,
                );
    }

    multi method to-adapted(Math::SparseMatrix::DOK:D $m --> Math::SparseMatrix::NativeAdapter:D) {
        return self.to-adapted($m.to-csr);
    }
}