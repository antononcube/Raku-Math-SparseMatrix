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
}