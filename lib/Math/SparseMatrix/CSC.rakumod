use Math::SparseMatrix::Abstract;

# This class, for Compressed Sparse Column (CSC) representation,
# is for used to make more convenient the interfacing with other packages or tools.
# I.e. it is not "operational", it is to do conversions with.

class Math::SparseMatrix::CSC is Math::SparseMatrix::Abstract {
    has @.values;
    has @.row-index;
    has @.col-ptr;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value is rw = 0;

    #=================================================================
    # Creators
    #=================================================================
    multi method new(:@values! where @values.all ~~ Numeric:D,
                     :@row-index, :@col-ptr,
                     UInt:D :$nrow, UInt:D :$ncol,
                     Numeric:D :$implicit-value = 0) {
        self.bless(:@values, :@row-index, :@col-ptr, :$nrow, :$ncol, :$implicit-value);
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose() {
        my @col-ptr = 0 xx ($!nrow + 1);

        # Each row of the original matrix becomes a column of its transpose.
        for @!row-index -> $row {
            @col-ptr[$row + 1]++;
        }
        for ^$!nrow -> $col {
            @col-ptr[$col + 1] += @col-ptr[$col];
        }

        my @next = @col-ptr.clone;
        my @values;
        my @row-index;
        for ^$!ncol -> $col {
            for @!col-ptr[$col] ..^ @!col-ptr[$col + 1] -> $idx {
                my $dest = @next[@!row-index[$idx]]++;
                @values[$dest] = @!values[$idx];
                @row-index[$dest] = $col;
            }
        }

        return Math::SparseMatrix::CSC.new(
                :@values,
                :@row-index,
                :@col-ptr,
                nrow => $!ncol,
                ncol => $!nrow,
                implicit-value => $!implicit-value
                );
    }
}
