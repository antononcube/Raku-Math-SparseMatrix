use Math::SparseMatrix::Abstract;

class Math::SparseMatrix::DOK is Math::SparseMatrix::Abstract {
    has %.rules;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value is rw = 0;

    #=================================================================
    # Creators
    #=================================================================
    multi method new(:@rules! where @rules.all ~~ Pair:D,
                     :$nrow is copy = @rules.map(*.key[0]).max + 1,
                     :$ncol is copy = @rules.map(*.key[1]).max + 1,
                     Numeric:D :$implicit-value = 0) {
        # There should be correctness check, etc.
        self.bless(rules => @rules.Hash, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:%rules!,
                     :$nrow is copy = %rules.map(*.key.words.head).max + 1,
                     :$ncol is copy = %rules.map(*.key.words.tail).max + 1,
                     Numeric:D :$implicit-value = 0) {
        # There should be correctness check, etc.
        self.bless(:%rules, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:@dense-matrix! where @dense-matrix ~~ List:D && @dense-matrix.all ~~ List:D,
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     Numeric:D :$implicit-value = 0) {
        my %rules;
        for @dense-matrix.kv -> $row, @cols {
            for @cols.kv -> $col, $val {
                if $val != $implicit-value && $row < $nrow && $col < $ncol {
                    %rules.push(Pair.new(($row, $col), $val));
                }
            }
        }

        self.bless(:%rules, :$nrow, :$ncol, :$implicit-value);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix::DOK.new(
                rules => %!rules.clone,
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
        return %!rules{($row, $col).Str} // $!implicit-value;
    }

    method row-at(Int:D $i --> Math::SparseMatrix::DOK) {
        # Can this be optimized?
        my %row = %!rules.grep({ $_.key.words.head == $i }).map({ (0, $_.key.words.tail) => $_.value });
        return Math::SparseMatrix::DOK.new(
                rules => %row,
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
        return self.transpose.row-at($col).transpose;
    }

    method AT-POS(*@index) {
        if @index.elems == 1 {
            return self.row-at(@index.head);
        }
        die "Only one index is expected.";
    }

    method !row-map(UInt $i) returns Hash {
        my %row = %!rules.grep({ $_.key.head == $i });
        return %row;
    }

    #=================================================================
    # Tuples
    #=================================================================
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
    method eqv(Math::SparseMatrix::DOK:D $other --> Bool:D) {
        if $!nrow != $other.nrow || $!ncol != $other.ncol ||
                %!rules.elems != $other.rules.elems ||
                $!implicit-value != $other.implicit-value {
            return False;
        }

        return self.rules eqv $other.rules;
    }

    #=================================================================
    # Info
    #=================================================================
    method adjacency-lists() {
        my @adj-lists;
        for @.row-ptr.kv -> $i, $ptr {
            my $next_ptr = $i == @.row-ptr.end ?? @.values.elems !! @.row-ptr[$i + 1];
            my @list;
            for $ptr ..^ $next_ptr -> $j {
                @list.push(@.col-index[$j])
            }
            @adj-lists.push(@list)
        }
        return @adj-lists;
    }

    method column-indices() {
        return %!rules.map({ $_.key.words[1] });
    }

    method columns-count() {
        return $!ncol;
    }

    method density() {
        return %!rules.elems / ($!nrow * $!ncol);
    }

    method dimensions() {
        return ($!nrow, $!ncol);
    }

    method explicit-length() {
        return %!rules.elems;
    }

    method explicit-positions() {
        return self.rules.keys;
    }

    method explicit-values() {
        return %!rules.values;
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

        my %resRules = %!rules.clone;
        my %resRules2 = $other.rules.clone.map({
            my ($i, $j) = $_.key.words>>.Int;
            ($i + $!nrow, $j) => $_.value
        });
        %resRules = %resRules, %resRules2;
        return Math::SparseMatrix::DOK.new(
                rules => %resRules,
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
        my %resRules = %!rules.grep({ $_.key.words.head < $n });
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
                @row.push(%.rules{($i, $j).Str} // $implicit-value);
            }
            @matrix.push(@row);
        }
        return @matrix;
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose(--> Math::SparseMatrix::DOK) {
        return Math::SparseMatrix::DOK.new(
                rules => %!rules.map({ $_.key.words.reverse => $_.value }).Hash,
                nrow => self.ncol,
                ncol => self.nrow
                );
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
    multi method dot(Math::SparseMatrix::DOK:D $other --> Math::SparseMatrix::DOK) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.nrow;

        my %result;
        for %.rules.kv -> $key, $value {
            my ($i, $j) = $key.words>>.Int;
            for 0 ..^ $other.ncol -> $k {
                %result{($i, $k).Str} += $value * ($other.rules{($j, $k).Str} // $other.implicit-value);
            }
        }

        return Math::SparseMatrix::DOK.new(
                :rules(%result),
                :$!nrow,
                :ncol($other.ncol),
                :implicit-value($!implicit-value * $other.implicit-value)
                );
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition
    multi method add(Numeric:D $a --> Math::SparseMatrix::DOK:D) {
        my %resRules = %!rules.map({ $_.key => $_.value + $a });
        return Math::SparseMatrix::DOK.bless(
                rules => %resRules,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $a
                );
    }

    #| Numeric addition of two matrices
    multi method add(Math::SparseMatrix::DOK $other --> Math::SparseMatrix::DOK:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        my %resRules;
        for %!rules -> $r1 {
            for $other.rules -> $r2 {
                if $r1.key eqv $r2.key {
                    %resRules{$r1.key.Str} = $r1.value + $r2.value
                } else {
                    if %resRules{$r1.key.Str}:!exists { %resRules{$r1.key.Str} = $r1.value + $other.implicit-value }
                    if %resRules{$r2.key.Str}:!exists { %resRules{$r2.key.Str} = $r2.value + $!implicit-value }
                }
            }
        }

        return Math::SparseMatrix::DOK.bless(
                rules => %resRules,
                :$!nrow,
                :$!ncol,
                :$!implicit-value
                );
    }

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication
    multi method multiply(Numeric:D $a --> Math::SparseMatrix::DOK:D) {
        my %resRalues = %!rules.map({ $_.key => $_.value * $a });
        return Math::SparseMatrix::DOK.bless(
                rules => %resRalues,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value * $a
                );
    }

    #| Matrix-matrix element-wise multiplication
    multi method multiply(Math::SparseMatrix::DOK:D $other --> Math::SparseMatrix::DOK:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        my %resRules;
        my %r2 = $other.rules;
        my @keys = [|%!rules.keys, |%r2.keys].unique;

        my $use-iv = !($!implicit-value == 0 || $other.implicit-value == 0);

        for @keys -> $pos {
            if (%!rules{$pos}:exists) && (%r2{$pos}:exists) {
                %resRules{$pos} = %!rules{$pos} * %r2{$pos};
            } elsif (%r2{$pos}:exists) && $use-iv {
                %resRules{$pos} = $!implicit-value * %r2{$pos};
            } elsif (%!rules{$pos}:exists) && $use-iv {
                %resRules{$pos} = %!rules{$pos} * $other.implicit-value;
            }
        }

        return Math::SparseMatrix::DOK.bless(
                rules => %resRules,
                :$!nrow,
                :$!ncol,
                :$!implicit-value
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
        %!rules = %!rules.key X=> 1;
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
        %!rules = %!rules.map({
            $_.key => do if $_.value < $v-min { $v-min }
            elsif $v-max < $_.value { $v-max }
            else { $_.value }
        });
        return self;
    }

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {
        my $max-length = %.rules.values.map(*.Str.chars).max // 1;
        my $default = $show-implicit-value ?? $!implicit-value !! '.';
        $max-length = max($max-length, $default.Str.chars);

        my @rows;
        for ^$!nrow -> $i {
            my @row = do for ^$!ncol -> $j {
                my $formatted = do if %!rules{($i, $j).Str}:exists {
                    %!rules{($i, $j).Str} ;
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