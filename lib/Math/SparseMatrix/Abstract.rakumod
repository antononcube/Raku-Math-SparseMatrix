class Math::SparseMatrix::Abstract {
    
    #=================================================================
    # Clone
    #=================================================================
    method clone() {...}
    
    #=================================================================
    # Access
    #=================================================================
    method elems(::?CLASS:D:) {...}
    
    method value-at(Int:D $row, Int:D $col) {...}

    method row-at(Int:D $i --> Math::SparseMatrix::Abstract) {...}

    method row-slice(*@indexes) {...}

    method column-at(Int:D $col --> Math::SparseMatrix::Abstract) {...}

    method AT-POS(*@index) {...}

    method !row-map(UInt $i) returns Hash {...}

    #=================================================================
    # Tuples
    #=================================================================
    # Same as Math::SparseMatrix::Abstract.tuples
    #| Tuples (or triplets)
    method tuples(Bool:D :d(:$dataset) = False) {...}

    #=================================================================
    # Equivalence
    #=================================================================
    method eqv(Math::SparseMatrix::Abstract:D $other --> Bool:D) {...}

    #=================================================================
    # Info
    #=================================================================
    method adjacency-lists() {...}

    method column-indices() {...}

    method columns-count() {...}

    method density() {...}

    method dimensions() {...}

    method explicit-length() {...}

    method explicit-positions() {...}

    method explicit-values() {...}

    method row-pointers() {...}

    method rows-count() {...}

    #=================================================================
    # Row-bind
    #=================================================================
    method row-bind(Math::SparseMatrix::Abstract:D $other --> Math::SparseMatrix::Abstract:D) {...}

    #=================================================================
    # Column-bind
    #=================================================================
    method column-bind(Math::SparseMatrix::Abstract:D $other --> Math::SparseMatrix::Abstract:D) {...}

    #=================================================================
    # Slicing
    #=================================================================
    method head(Int $n = 1 --> Math::SparseMatrix::Abstract:D) {...}

    #=================================================================
    # Dense array
    #=================================================================
    #| (Dense) array of arrays representation.
    #| C<:$implicit-value> -- Implicit value to use.
    method Array(:i(:iv(:$implicit-value)) is copy = Whatever) {...}

    #=================================================================
    # Transpose
    #=================================================================
    method transpose(--> Math::SparseMatrix::Abstract) {...}

    #=================================================================
    # Matrix-vector multiplication
    #=================================================================
    #| Dot product of a sparse matrix with dense vector
    multi method dot(@vector --> Array) {...}

    #=================================================================
    # Matrix-matrix multiplication
    #=================================================================
    #| Dot product of two sparse matrices
    multi method dot(Math::SparseMatrix::Abstract:D $other --> Math::SparseMatrix::Abstract) {...}

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition
    multi method add(Numeric:D $a --> Math::SparseMatrix::Abstract:D) {...}

    #| Numeric addition of two matrices
    multi method add(Math::SparseMatrix::Abstract $other --> Math::SparseMatrix::Abstract:D) {...}

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication
    multi method multiply(Numeric:D $a --> Math::SparseMatrix::Abstract:D) {...}

    #| Matrix-matrix element-wise multiplication
    multi method multiply(Math::SparseMatrix::Abstract:D $other --> Math::SparseMatrix::Abstract:D) {...}

    #=================================================================
    # Unitize
    #=================================================================
    #| Unitize the sparse matrix
    #| C<:$clone> -- Whether to operate in-place.
    method unitize(Bool:D :$clone = True) {...}

    #=================================================================
    # Clip
    #=================================================================
    #| Clip the sparse matrix
    #| C<:$v-min> -- The new min value.
    #| C<:$v-max> -- The new max value.
    #| C<:$clone> -- Whether to operate in-place.
    method clip(Numeric:D :min(:$v-min)!, Numeric:D :max(:$v-max)!, Bool:D :$clone = True) {...}

    #=================================================================
    # Round
    #=================================================================
    #| Round the sparse matrix
    #| C<:$scale> -- Scale to round to.
    method round(Numeric:D $scale = 1, Bool:D :$clone = True) {...}

    #=================================================================
    # Conjugate
    #=================================================================
    #| Round the sparse matrix
    method conjugate(Bool:D :$clone = True) {...}

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {...}

    #=================================================================
    # Representation
    #=================================================================
    #| Convert to adjacency maps.
    method to-adjacency-list(%rules) {
        my %hash-of-hashes;
        for %rules.kv -> $key, $value {
            my ($row, $col) = $key.words;
            %hash-of-hashes{$row}{$col} = $value;
        }
        return %hash-of-hashes;
    }

    #| From adjacency map to dictionary of keys
    method to-rules(%adjacency-map) {
        my %rules;
        for %adjacency-map.kv -> $i, %row {
            for %row.kv -> $j, $val {
                %rules{($i, $j).Str} = $val;
            }
        }
        return %rules;
    }

    # As Math::SparseMatrix::CSR.wl
    #| Wolfram Language (WL) representation
    method to-wl() {...}

    # As Math::SparseMatrix::CSR.Hash
    #| To Hash
    multi method Hash(::?CLASS:D:-->Hash) {...}

    #| To string
    multi method Str(::?CLASS:D:-->Str) {
        return self.gist;
    }

    #| To gist
    multi method gist(::?CLASS:D:-->Str) {...}
}