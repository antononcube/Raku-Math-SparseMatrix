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
    # Row and column sums and maxes
    #=================================================================
    method !row-op(&op, Bool:D :p(:$pairs) = False) {
        # Probably more effective implementations can be provided by the core matrix classes.
        # But this is universal and quick to implement.
        my @sums = do for ^self.nrow -> $i {
            self.core-matrix.row-at($i).tuples>>.tail.&op
        }

        if $pairs {
            my @rn = self.row-names.defined ?? self.row-names !! (^self.nrow);
            return (@rn.Array Z=> @sums.Array).Hash;
        }
        return @sums;
    }

    #| Row sums for a sparse matrix.
    method row-sums(Bool:D :p(:$pairs) = False) {
        my &asum = -> @a { @a.sum + self.implicit-value * (self.ncol - @a.elems) };
        self!row-op(&asum, :$pairs);
    }

    #| Column sums for a sparse matrix.
    method column-sums(Bool:D :p(:$pairs) = False) {
        # Quick to implement.
        return self.transpose.row-sums(:$pairs);
    }

    #| Row maxes for a sparse matrix.
    method row-maxes(Bool:D :p(:$pairs) = False) {
        my &amax = -> @a { @a.Array.push(self.implicit-value).max };
        self!row-op(&amax, :$pairs);
    }

    #| Column for a sparse matrix.
    method column-maxes(Bool:D :p(:$pairs) = False) {
        # Quick to implement.
        return self.transpose.row-maxes(:$pairs);
    }

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
    #| C<:$clone> -- Whether to operate in-place.
    method round(Numeric:D $scale = 1, Bool:D :$clone = True) {...}

    #=================================================================
    # Top-k elements matrix
    #=================================================================
    #| Give a new sparse matrix with the largest, top-k elements only
    #| C<$k> -- Number of the top elements.
    #| C<:$clone> -- Whether to operate in-place.
    method top-k-elements-matrix(UInt:D $k, Bool:D :$clone = True) {...}

    #=================================================================
    # Apply elementwise
    #=================================================================
    #| Apply a function to the elements of sparse matrix.
    #| C<&func> -- Function to apply.
    #| C<:$skip-implicit-value> -- Should application to the implicit value be skipped or not?
    #| C<:$clone> -- Whether to operate in-place.
    method apply-elementwise(&func, Bool:D :$skip-implicit-value = False, Bool:D :$clone = True) {...}

    #=================================================================
    # Conjugate
    #=================================================================
    #| Conjugate the sparse matrix
    method conjugate(Bool:D :$clone = True) {...}

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {...}

    #=================================================================
    # Representation
    #=================================================================
    #| Convert to adjacency maps.
    method to-adjacency-map(%rules) {
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