use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;

class Math::SparseMatrix::DOK is Math::SparseMatrix::Abstract {
    has %.adjacency-map;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value is rw = 0;

    #=================================================================
    # Creators
    #=================================================================
    multi method new(:%adjacency-map!,
                     :$nrow is copy = Whatever,
                     :$ncol is copy = Whatever,
                     Numeric:D :$implicit-value = 0) {
        if $nrow.isa(Whatever) { $nrow = %adjacency-map.keys.max }
        if $ncol.isa(Whatever) { $ncol = %adjacency-map.map({ $_.keys }).map(*.Slip).max }
        self.bless(:%adjacency-map, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:@rules! where @rules.all ~~ Pair:D,
                     :$nrow is copy = @rules.map(*.key[0]).max + 1,
                     :$ncol is copy = @rules.map(*.key[1]).max + 1,
                     Numeric:D :$implicit-value = 0) {
        # There should be correctness check, etc.
        self.new(rules => @rules.Hash, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:%rules!,
                     :$nrow is copy = %rules.map(*.key.words.head).max + 1,
                     :$ncol is copy = %rules.map(*.key.words.tail).max + 1,
                     Numeric:D :$implicit-value = 0) {
        # There should be correctness check, etc.
        my %adjacency-map = self.to-adjacency-map(%rules);
        self.bless(:%adjacency-map, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:@dense-matrix! where @dense-matrix ~~ List:D && @dense-matrix.all ~~ List:D,
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     Numeric:D :$implicit-value = 0) {
        my %adjacency-map;
        for @dense-matrix.kv -> $row, @cols {
            for @cols.kv -> $col, $val {
                if $val != $implicit-value && $row < $nrow && $col < $ncol {
                    %adjacency-map{$row}{$col} = $val;
                }
            }
        }
        self.bless(:%adjacency-map, :$nrow, :$ncol, :$implicit-value);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix::DOK.new(
                adjacency-map => %!adjacency-map.clone,
                :$!nrow,
                :$!ncol,
                :$!implicit-value
                );
    }


    #=================================================================
    # Access
    #=================================================================
    method elems(::?CLASS:D:) {
        return $!nrow;
    }

    method value-at(Int:D $row, Int:D $col) {
        return %!adjacency-map{$row}{$col} // $!implicit-value;
    }

    method row-at(Int:D $i --> Math::SparseMatrix::DOK) {
        return Math::SparseMatrix::DOK.new(
                adjacency-map => %!adjacency-map{$i}:exists ?? {0 => %!adjacency-map{$i}} !! %(),
                nrow => 1,
                :$!ncol,
                :$!implicit-value
                );
    }

    method row-slice(*@indexes) {
        die 'The indexes are expected to be non-negative integers.'
        unless (@indexes.all ~~ Int:D) && min(@indexes) ≥ 0;

        my @mats = @indexes.map({ self.row-at($_) });
        my $res = @mats.head;
        for @mats.tail(*- 1) -> $m {
            $res = $res.row-bind($m)
        }
        return $res;
    }

    method column-at(Int:D $col --> Math::SparseMatrix::DOK) {
        # return self.transpose.row-at($col).transpose;
        my %res-adjacency-map = %!adjacency-map.map({ $_.value{$col}:exists ?? ($_.key => %(0 => $_.value{$col})) !! Empty });
        return Math::SparseMatrix::DOK.new(
                adjacency-map => %res-adjacency-map,
                :$!nrow,
                ncol => 1,
                :$!implicit-value
                );
    }

    method AT-POS(*@index) {
        if @index.elems == 1 {
            return self.row-at(@index.head);
        }
        die "Only one index is expected.";
    }

    #=================================================================
    # Rules and tuples
    #=================================================================
    method rules() {
        return self.to-rules(%!adjacency-map);
    }

    # Same as Math::SparseMatrix::DOK.tuples
    #| Tuples (or triplets)
    method tuples(Bool:D :d(:$dataset) = False) {
        my @res = self.rules.map({ [|$_.key, $_.value].List }).Array;
        if $dataset {
            @res = @res.map({ <i j x>.Array Z=> $_.Array })>>.Hash.Array
        }
        return @res;
    }

    #=================================================================
    # Equivalence
    #=================================================================
    method eqv(Math::SparseMatrix::DOK:D $other, Numeric:D :$tol = 1e-14 --> Bool:D) {
        return False unless $!nrow == $other.nrow && $!ncol == $other.ncol;
        for %!adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $value {
                return False unless abs($value - ($other.adjacency-map{$row}{$col} // $other.implicit-value)) <= $tol;
            }
        }
        for $other.adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $value {
                return False unless abs($value - (%!adjacency-map{$row}{$col} // $!implicit-value)) <= $tol;
            }
        }
        return True;
    }

    #=================================================================
    # Info
    #=================================================================
    method adjacency-lists() {
        return do for ^$!nrow -> $r { %!adjacency-map{$r}.keys.sort // ()};
    }

    method column-indices() {
        return %!adjacency-map.map({ $_.keys }).map(*.Slip).unique;
    }

    method columns-count() {
        return $!ncol;
    }

    method density() {
        return %!adjacency-map.values>>.elems.sum / ($!nrow * $!ncol);
    }

    method dimensions() {
        return ($!nrow, $!ncol);
    }

    method explicit-length() {
        return %!adjacency-map.values>>.elems.sum;
    }

    method explicit-positions() {
        return self.rules.keys;
    }

    method explicit-values() {
        return %!adjacency-map.values>>.values.map(*.Slip);
    }

    method row-pointers() {
        return Empty;
    }

    method rows-count() {
        return $!nrow;
    }

    #=================================================================
    # Row-bind
    #=================================================================
    method row-bind(Math::SparseMatrix::DOK:D $other --> Math::SparseMatrix::DOK:D) {
        die 'The number of columns of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.ncol;

        die 'The the implicit value of the argument is expected to be same as the explicit value of the object.'
        unless $!implicit-value == $other.implicit-value;

        my %res-adjacency-map = %!adjacency-map.clone;
        my %res-adjacency-map2 = $other.adjacency-map.map({
            ($_.key + $!nrow) => $_.value
        });
        %res-adjacency-map = %res-adjacency-map , %res-adjacency-map2;
        return Math::SparseMatrix::DOK.new(
                adjacency-map => %res-adjacency-map,
                nrow => $!nrow + $other.nrow,
                :$!ncol,
                :$!implicit-value
                );
    }

    #=================================================================
    # Column-bind
    #=================================================================
    method column-bind(Math::SparseMatrix::DOK:D $other --> Math::SparseMatrix::DOK:D) {
        die 'The number of rows of the argument is expected to be equal to the number of rows of the object.'
        unless $!nrow == $other.nrow;

        # Not very effective, but quick to implement
        return self.transpose.row-bind($other.transpose).transpose;
    }

    #=================================================================
    # Slicing
    #=================================================================
    method head(Int $n = 1 --> Math::SparseMatrix::DOK:D) {
        my %resRules = self.rules.grep({ $_.key.words.head < $n });
        return Math::SparseMatrix::DOK(rules => %resRules);
    }

    #=================================================================
    # Dense array
    #=================================================================
    #| (Dense) array of arrays representation.
    #| C<:$implicit-value> -- Implicit value to use.
    method Array(:i(:iv(:$implicit-value)) is copy = Whatever) {
        if $implicit-value.isa(Whatever) { $implicit-value = $!implicit-value }
        my @matrix;
        for ^$!nrow -> $i {
            my @row;
            for ^$!ncol -> $j {
                @row.push(%!adjacency-map{$i}{$j} // $implicit-value);
            }
            @matrix.push(@row);
        }
        return @matrix;
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose(--> Math::SparseMatrix::DOK) {
        my %transposed;
        for %!adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $value {
                %transposed{$col}{$row} = $value;
            }
        }
        self.new(adjacency-map => %transposed, :nrow($!ncol), :ncol($!nrow), :$!implicit-value);
    }

    #=================================================================
    # Matrix-vector multiplication
    #=================================================================
    #| Dot product of a sparse matrix with dense vector
    multi method dot(@vector --> Array) {
        die "Vector length does not match matrix column count."
        unless @vector.elems == $!ncol;

        my @result = 0 xx $!nrow;
        for %.rules.kv -> $key, $value {
            my ($i, $j) = $key.words>>.Int;
            @result[$i] += $value * @vector[$j];
        }
        return @result;
    }

    #=================================================================
    # Matrix-matrix multiplication
    #=================================================================
    #| Dot product of two sparse matrices
    multi method dot(Math::SparseMatrix::DOK $other --> Math::SparseMatrix::DOK) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.nrow;

        my %result;
        for %!adjacency-map.kv -> $i, %row {
            for %row.kv -> $j, $val {
                for $other.adjacency-map{$j}.kv -> $k, $other-val {
                    %result{$i}{$k} += $val * $other-val;
                }
            }
        }
        return Math::SparseMatrix::DOK.new(
                adjacency-map => %result,
                :$!nrow,
                ncol => $other.ncol,
                implicit-value => $!implicit-value * $other.implicit-value
                );
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition
    multi method add(Numeric:D $a --> Math::SparseMatrix::DOK:D) {
        my %res-adjacency-map;
        for %.adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                %res-adjacency-map{$row}{$col} = $val + $a;
            }
        }
        return Math::SparseMatrix::DOK.bless(
                adjacency-map => %res-adjacency-map,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $a
                );
    }

    #| Numeric addition of two matrices
    multi method add(Math::SparseMatrix::DOK $other --> Math::SparseMatrix::DOK:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        my %res-adjacency-map;
        for %!adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                %res-adjacency-map{$row}{$col} = $val + ($other.adjacency-map{$row}{$col} // $other.implicit-value);
            }
        }

        for $other.adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                %res-adjacency-map{$row}{$col} //= $val + $!implicit-value;
            }
        }

        return Math::SparseMatrix::DOK.bless(
                adjacency-map => %res-adjacency-map,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $other.implicit-value
                );
    }

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication
    multi method multiply(Numeric:D $a --> Math::SparseMatrix::DOK:D) {
        my %res-adjacency-map;
        for %.adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                %res-adjacency-map{$row}{$col} = $val * $a;
            }
        }
        return Math::SparseMatrix::DOK.bless(
                adjacency-map => %res-adjacency-map,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value * $a
                );
    }

    #| Matrix-matrix element-wise multiplication
    multi method multiply(Math::SparseMatrix::DOK:D $other --> Math::SparseMatrix::DOK:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        my %res-adjacency-map;
        for %!adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                if $other.adjacency-map{$row}{$col}:exists {
                    %res-adjacency-map{$row}{$col} = $val * $other.adjacency-map{$row}{$col};
                } else {
                    %res-adjacency-map{$row}{$col} = $val * $other.implicit-value;
                }
            }
        }

        for $other.adjacency-map.kv -> $row, %cols {
            for %cols.kv -> $col, $val {
                unless %!adjacency-map{$row}{$col}:exists {
                    %res-adjacency-map{$row}{$col} = $.implicit-value * $val;
                }
            }
        }
        return Math::SparseMatrix::DOK.bless(
                adjacency-map => %res-adjacency-map,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value * $other.implicit-value
                );
    }

    #=================================================================
    # Unitize
    #=================================================================
    #| Unitize the sparse matrix
    #| C<:$clone> -- Whether to operate in-place.
    method unitize(Bool:D :$clone = True) {
        if $clone {
            return self.clone.unitize(:!clone);
        }
        %!adjacency-map = %!adjacency-map.map({ $_.map({ $_.key => 1 }) });
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
        %!adjacency-map = %!adjacency-map.map({
            $_.map({
                $_.key => do if $_.value < $v-min { $v-min }
                elsif $v-max < $_.value { $v-max }
                else { $_.value }
            })
        });
        return self;
    }

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {
        my $max-length = %!adjacency-map.values>>.values.map(*.Slip).map(*.Str.chars).max // 1;
        my $default = $show-implicit-value ?? $!implicit-value !! '.';
        $max-length = max($max-length, $default.Str.chars);

        my @rows;
        for ^$!nrow -> $i {
            my @row = do for ^$!ncol -> $j {
                my $formatted = do if %!adjacency-map{$i}{$j}:exists {
                    %!adjacency-map{$i}{$j} ;
                } else {
                    $default;
                }
                sprintf "%-*s ", $max-length, $formatted;
            }
            @rows.push: @row
        }

        if $echo {
            for @rows -> @row {
                say @row.join(' ');
            }
        }

        return @rows>>.Array.Array;
    }

    #=================================================================
    # Representation
    #=================================================================
    #| To CSR sparse matrix.
    method to-csr() {
        my @values;
        my @col-index;
        my @row-ptr = 0;
        my $nnz = 0;

        for ^$!nrow -> $row {
            for %.adjacency-map{$row}.kv -> $col, $value {
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
                :$!nrow,
                :$!ncol,
                :$!implicit-value
                );
    }

    # As Math::SparseMatrix::CSR.wl
    #| Wolfram Language (WL) representation
    method wl() {
        my $rules = self.rules.map({ "\{{ $_.key.words.head + 1 },{ $_.key.words.tail + 1 }\}->{ $_.value }" }).join(',');
        return "SparseArray[\{$rules\}, \{{ $!nrow }, { $!ncol }\}, { $!implicit-value }]"
    }


    # As Math::SparseMatrix::CSR.Hash
    #| To Hash
    multi method Hash(::?CLASS:D:-->Hash) {
        return
                {
                    specified-elements => self.explicit-length,
                    dimensions => ($!nrow, $!ncol),
                    default => $!implicit-value,
                    density => self.density,
                    elements => self.rules,
                };
    }

    #| To string
    multi method Str(::?CLASS:D:-->Str) {
        return self.gist;
    }

    #| To gist
    multi method gist(::?CLASS:D:-->Str) {
        return 'Math::SparseMatrix::DOK' ~ (specified-elements => self.explicit-length, dimensions => ($!nrow, $!ncol),
                                            density => self.density).List.raku;
    }
}