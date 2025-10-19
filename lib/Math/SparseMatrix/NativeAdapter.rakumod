use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::Native;
use NativeHelpers::Array;

class Math::SparseMatrix::NativeAdapter is Math::SparseMatrix::Abstract {

    has Math::SparseMatrix::Native::CSRStruct:D $.csr-struct
            is required
            handles(
            "nrow", "ncol", "nnz",
            "value-at",
            "adjacency-lists", "columns-count", "explicit-length",
            "density", "dimensions", "implicit-value" => "implicit_value",
            "rows-count",
            "rules", "tuples",
            "row-sums", "row-maxes", "column-sums", "column-maxes"
            );

    # In principle we could use these handles:
    #   "clone", "row-at", "row-ptr" => "row_ptr", "col-index" => "col_index", "values"
    # but they give C-array and it is better to "standard" arrays.

    #=================================================================
    # Creators
    #=================================================================
    multi method new(Math::SparseMatrix::Native::CSRStruct:D $csr-struct) {
        self.bless(:$csr-struct);
    }

    multi method new(Math::SparseMatrix::Native::CSRStruct:D :$csr-struct!) {
        self.bless(:$csr-struct);
    }

    multi method new(:@values! where @values.all ~~ Numeric:D,
                     :@col-index, :@row-ptr,
                     UInt:D :$nrow, UInt:D :$ncol,
                     Numeric:D :$implicit-value = 0) {
        my $csr-struct =
                Math::SparseMatrix::Native::CSRStruct.new(
                :@values,
                col_index => @col-index,
                row_ptr => @row-ptr,
                :$nrow,
                :$ncol,
                implicit_value => $implicit-value);
        self.bless(:$csr-struct);
    }

    multi method new(:@rules! where @rules.all ~~ Pair:D,
                     :$nrow is copy = @rules.map(*.key[0]).max + 1,
                     :$ncol is copy = @rules.map(*.key[1]).max + 1,
                     Numeric:D :$implicit-value = 0) {
        my $csr-struct =
                Math::SparseMatrix::Native::CSRStruct.new(
                        :@rules,
                        :$nrow,
                        :$ncol,
                        implicit_value => $implicit-value);
        self.bless(:$csr-struct);
    }

    multi method new(:@dense-matrix! where @dense-matrix ~~ List:D && @dense-matrix.all ~~ List:D,
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     Numeric:D :$implicit-value = 0) {
        my $csr-struct =
                Math::SparseMatrix::Native::CSRStruct.new(
                        :@dense-matrix,
                        :$nrow,
                        :$ncol,
                        implicit_value => $implicit-value);
        self.bless(:$csr-struct);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        Math::SparseMatrix::NativeAdapter.new(csr-struct => $!csr-struct.clone);
    }

    #=================================================================
    # Access
    #=================================================================

    method row-ptr() { return copy-to-array($!csr-struct.row_ptr, $!csr-struct.nrow + 1)}
    method col-index() { return copy-to-array($!csr-struct.col_index, $!csr-struct.nnz)}
    method values() { return copy-to-array($!csr-struct.values, $!csr-struct.nnz)}

    # Delegated
    #method elems(::?CLASS:D:) {...}

    # Delegated
    #method value-at(Int:D $row, Int:D $col) {...}

    method row-at(Int:D $i --> Math::SparseMatrix::Abstract) {
        my $csr-struct = $!csr-struct.row-at($i);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    method row-slice(*@indexes) {
        my $csr-struct = $!csr-struct.row-slice(@indexes);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    method column-at(Int:D $col --> Math::SparseMatrix::Abstract) {
        my $csr-struct = $!csr-struct.column-at($col);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    method AT-POS(*@index) {
        my $csr-struct = $!csr-struct.AT-POS(@index);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    method !row-map(UInt:D $i) returns Hash {
        my %row;
        for $!csr-struct.row-ptr[$i] ..^ $!csr-struct.row-ptr[$i + 1] -> $j {
            %row{$!csr-struct.col-index[$j]} = $!csr-struct.values[$j];
        }
        return %row;
    }

    #=================================================================
    # Tuples
    #=================================================================
    # Same as Math::SparseMatrix::Abstract.tuples
    #| Tuples (or triplets)
    #method tuples(Bool:D :d(:$dataset) = False) {...}

    #=================================================================
    # Equivalence
    #=================================================================
    method eqv(Math::SparseMatrix::NativeAdapter:D $other --> Bool:D) {
        return $!csr-struct.eqv($other.csr-struct);
    }

    #=================================================================
    # Info
    #=================================================================

    # Those methods are delegated.

    #=================================================================
    # Row-bind
    #=================================================================
    method row-bind(Math::SparseMatrix::NativeAdapter:D $other --> Math::SparseMatrix::NativeAdapter:D) {
        my $csr-struct = $!csr-struct.row-bind($other.csr-struct);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Column-bind
    #=================================================================
    method column-bind(Math::SparseMatrix::NativeAdapter:D $other --> Math::SparseMatrix::NativeAdapter:D) {
        my $csr-struct = $!csr-struct.row-bind($other.csr-struct);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Slicing
    #=================================================================
    method head(Int $n = 1 --> Math::SparseMatrix::NativeAdapter:D) {
        die "The first element is expected to be a positive integer."
        unless $n ~~ Int:D && $n > 0;

        return $n == 1 ?? self[0] !! self.row-slice(^$n);
    }

    #=================================================================
    # Dense array
    #=================================================================
    #| (Dense) array of arrays representation.
    #| C<:$implicit-value> -- Implicit value to use.
    method Array(:i(:iv(:$implicit-value)) is copy = Whatever) {
        return $!csr-struct.Array(:$implicit-value);
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose(--> Math::SparseMatrix::Abstract) {
        my $csr-struct = $!csr-struct.transpose;
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Matrix-vector multiplication
    #=================================================================
    #| Dot product of a sparse matrix with dense vector
    multi method dot(@vector --> Array) {
        return $!csr-struct.dot(@vector);
    }

    #=================================================================
    # Matrix-matrix multiplication
    #=================================================================
    #| Dot product of two sparse matrices
    multi method dot(Math::SparseMatrix::Abstract:D $other --> Math::SparseMatrix::Abstract) {
        my $csr-struct = $!csr-struct.dot($other.csr-struct);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition
    multi method add(Numeric:D $a --> Math::SparseMatrix::Abstract:D) {
        my $csr-struct = $!csr-struct.add($a);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #| Numeric addition of two matrices
    multi method add(Math::SparseMatrix::Abstract $other --> Math::SparseMatrix::Abstract:D) {
        my $csr-struct = $!csr-struct.add($other.csr-struct);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication
    multi method multiply(Numeric:D $a --> Math::SparseMatrix::NativeAdapter:D) {
        my $csr-struct = $!csr-struct.multiply($a);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #| Matrix-matrix element-wise multiplication
    multi method multiply(Math::SparseMatrix::NativeAdapter:D $other --> Math::SparseMatrix::NativeAdapter:D) {
        my $csr-struct = $!csr-struct.multiply($other.csr-struct);
        return Math::SparseMatrix::NativeAdapter.new(:$csr-struct);
    }

    #=================================================================
    # Row sums and maxes
    #=================================================================
    # Handled by delegation

    #=================================================================
    # Unitize
    #=================================================================
    #| Unitize the sparse matrix
    #| C<:$clone> -- Whether to operate in-place.
    method unitize(Bool:D :$clone = True) {
        if $clone {
            return self.clone.unitize(:!clone);
        }
        self.csr-struct.unitize;
        return self;
    }

    #=================================================================
    # Clip
    #=================================================================
    #| Clip the sparse matrix
    #| C<:$v-min> -- The new min value.
    #| C<:$v-max> -- The new max value.
    #| C<:$clone> -- Whether to operate in-place.
    method clip(Numeric:D :min(:$v-min)!, Numeric:D :max(:$v-max)!, Bool:D :$clone = True) {
        if $clone {
            return self.clone.clip(:$v-min, :$v-max, :!clone);
        }
        self.csr-struct.clip(:$v-min, :$v-max);
        return self;
    }

    #=================================================================
    # Round
    #=================================================================
    #| Round the sparse matrix
    #| C<:$scale> -- Scale to round to.
    #| C<:$clone> -- Whether to operate in-place.
    method round(Numeric:D $scale = 1, Bool:D :$clone = True) {
        if $clone {
            return self.clone.round($scale, :!clone);
        }
        self.csr-struct.round($scale);
        return self;
    }

    #=================================================================
    # Top-k elements matrix
    #=================================================================
    #| Modify the sparse matrix or give a new sparse matrix with the largest, top-k elements only.
    #| C<$k> -- Number of the top elements.
    #| C<:$clone> -- Whether to operate in-place.
    method top-k-elements-matrix(UInt:D $k, Bool:D :$clone = True) {
        note 'In place top-K elements matrix modification is not implemented yet.' unless $clone;
        my $csr-struct = self.csr-struct.top-k-elements-matrix($k);
        return self.new(:$csr-struct);
    }

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {
        # Delegate to Math::SparseMatrix::CSR
        return self.to-csr.print;
    }

    #=================================================================
    # Representation
    #=================================================================
    method to-csr(-->Math::SparseMatrix::CSR:D) {
        return Math::SparseMatrix::CSR.new(
                values => self.values,
                row-ptr => self.row-ptr,
                col-index => self.col-index,
                nrow => self.nrow,
                ncol => self.ncol,
                implicit-value => self.implicit-value,
                );
    }

    # As Math::SparseMatrix::CSR.wl
    #| Wolfram Language (WL) representation
    method to-wl() {
        return self.to-csr.to-wl;
    }

    # As Math::SparseMatrix::CSR.Hash
    #| To Hash
    multi method Hash(::?CLASS:D:-->Hash) {
        return self.csr-struct.Hash;
    }

    #| To string
    multi method Str(::?CLASS:D:-->Str) {
        return self.csr-struct.Str;
    }

    #| To gist
    multi method gist(::?CLASS:D:-->Str) {
        return self.csr-struct.gist
    }
}
