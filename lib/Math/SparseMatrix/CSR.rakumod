use Math::SparseMatrix::Abstract;

class Math::SparseMatrix::CSR is Math::SparseMatrix::Abstract {
    has @.values;
    has @.col-index;
    has @.row-ptr;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value is rw = 0;

    #=================================================================
    # Creators
    #=================================================================
    multi method new(:@values! where @values.all ~~ Numeric:D,
                     :@col-index, :@row-ptr,
                     UInt:D :$nrow, UInt:D :$ncol,
                     Numeric:D :$implicit-value = 0) {
        self.bless(:@values, :@col-index, :@row-ptr, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:@rules! where @rules.all ~~ Pair:D,
                     :$nrow is copy = @rules.map(*.key[0]).max + 1,
                     :$ncol is copy = @rules.map(*.key[1]).max + 1,
                     Numeric:D :$implicit-value = 0) {
        my $msg = 'Each key of the argument @rules is expected to be a list of two non negative integers.';

        my %temp;
        for @rules -> $rule {
            die $msg unless ($rule.key ~~ Positional:D) && $rule.key.elems == 2;

            my ($row, $col) = $rule.key;

            die $msg unless ($row ~~ Int:D) && $row ≥0 && ($col ~~ Int:D) && $col ≥ 0;

            %temp{$row}{$col} = $rule.value;
        }

        my @values;
        my @col-index;
        my @row-ptr = 0;

        for ^$nrow -> $row {
            if %temp{$row}:exists {
                for %temp{$row}.keys.sort -> $col {
                    @values.push: %temp{$row}{$col};
                    @col-index.push: $col.Int;
                }
            }
            @row-ptr.push: @values.elems;
        }

        self.bless(:@values, :@col-index, :@row-ptr, :$nrow, :$ncol, :$implicit-value);
    }

    multi method new(:@dense-matrix! where @dense-matrix ~~ List:D && @dense-matrix.all ~~ List:D,
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     Numeric:D :$implicit-value = 0) {
        my @values;
        my @col-index;
        my @row-ptr = 0;

        for @dense-matrix.kv -> $row, @cols {
            for @cols.kv -> $col, $val {
                if $val != $implicit-value && $row < $nrow && $col < $ncol {
                    @values.push: $val;
                    @col-index.push: $col;
                }
            }
            @row-ptr.push: @values.elems;
        }

        if $nrow > @dense-matrix.elems {
            @row-ptr.append( @row-ptr.tail xx ($nrow - @dense-matrix.elems))
        }
        self.bless(:@values, :@col-index, :@row-ptr, :$nrow, :$ncol, :$implicit-value);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix::CSR.new(
                values => @!values.clone,
                col-index => @!col-index.clone,
                row-ptr => @!row-ptr.clone,
                :$!nrow,
                :$!ncol,
                :$!implicit-value
                );
    }

    #=================================================================
    # Verify
    #=================================================================
    method verify(Bool:D :$pairs = False) {
        my %res =
                row-ptr-elems => @!row-ptr.elems == $!nrow + 1,
                row-ptr-head => @!row-ptr.head == 0,
                row-ptr-tail => @!row-ptr.tail == @!values.elems,
                values-elems => @!values.elems == @!col-index.elems;
        return $pairs ?? %res !! [&&] %res.values;
    }

    #=================================================================
    # Access
    #=================================================================
    method elems(::?CLASS:D:) {
        return $!nrow;
    }

    method value-at(Int:D $row, Int:D $col) {
        for @!row-ptr[$row] ..^ @!row-ptr[$row + 1] -> $ip {
            return @!values[$ip] if @!col-index[$ip] == $col;
        }
        return $!implicit-value;
    }

    method row-at(Int:D $row --> Math::SparseMatrix::CSR) {
        my @values;
        my @row-ptr;
        my @col-index;

        for @!row-ptr[$row] ..^ @!row-ptr[$row + 1] -> $ip {
            @values[$ip - @!row-ptr[$row]] = @!values[$ip];
            @col-index[$ip - @!row-ptr[$row]] = @!col-index[$ip];
        }
        @row-ptr = [0, @values.elems];
        return Math::SparseMatrix::CSR.new(
                :@values,
                :@col-index,
                :@row-ptr,
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
        for @mats.tail(*-1) -> $m {
            $res = $res.row-bind($m)
        }
        return $res;
    }

    method column-at(Int:D $col --> Math::SparseMatrix::CSR) {
        # Not effective, but very quick to implement.
        return self.transpose.row-at($col).transpose;
    }

    method AT-POS(*@index) {
        if @index.elems == 1 {
            return self.row-at(@index.head);
        }
        die "Only one index is expected.";
    }

    method !row-map(UInt $i) returns Hash {
        my %row;
        for @!row-ptr[$i] ..^ @!row-ptr[$i + 1] -> $j {
            %row{@!col-index[$j]} = @!values[$j];
        }
        return %row;
    }

    #=================================================================
    # Rules and tuples
    #=================================================================
    #| Rules in the for C<(row, column) => value>.
    method rules() {
        my @rules;
        for @.row-ptr.kv -> $i, $ptr {
            my $next_ptr = $i == @.row-ptr.end ?? @.values.elems !! @.row-ptr[$i + 1];
            for $ptr ..^ $next_ptr -> $j {
                # We show all explicit values
                # i.e. using this check is wrong: if @.values[$j] != $!implicit-value { ...
                @rules.push(Pair.new(($i, @.col-index[$j]), @.values[$j]))
            }
        }
        return @rules;
    }

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
    method eqv(Math::SparseMatrix::CSR:D $other, Numeric:D :$tol = 1e-14 --> Bool:D) {
        if $!nrow != $other.nrow || $!ncol != $other.ncol ||
                @!values.elems != $other.values.elems ||
                $!implicit-value != $other.implicit-value {
            return False;
        }

        # Ineffective, but quick to implement
        #return self.rules.Hash eqv $other.rules.Hash;
        # This "quick" but it is not correct --
        # some matrices with "transposed" values within a row can give True.
        if @!row-ptr ne $other.row-ptr { return False; }
        if @!col-index.sort ne $other.col-index.sort { return False; }
        my Numeric $diff = (@!values.sort <<->> $other.values.sort)>>.abs.max;
        return $diff ≤ $tol;
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
        return @!col-index.clone;
    }

    method columns-count() {
        return $!ncol;
    }

    method density() {
        return @!values.elems / ($!nrow * $!ncol);
    }

    method dimensions() {
        return ($!nrow, $!ncol);
    }

    method explicit-length() {
        return @!values.elems;
    }

    method explicit-positions() {
        return self.rules.keys;
    }

    method explicit-values() {
        return @!values.clone;
    }

    method row-pointers() {
        return @!row-ptr.clone;
    }

    method rows-count() {
        return $!nrow;
    }

    #=================================================================
    # Row-bind
    #=================================================================
    method row-bind(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        die 'The number of columns of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.ncol;

        die 'The the implicit value of the argument is expected to be same as the explicit value of the object.'
        unless $!implicit-value == $other.implicit-value;

        my @values = @!values.clone.append($other.values);
        my @col-index = @!col-index.clone.append($other.col-index);
        my @row-ptr = @!row-ptr.clone.Array.append($other.row-ptr.tail(*- 1).map({ $_ + @!row-ptr[*- 1] }));

        return Math::SparseMatrix::CSR.new(
                :@values,
                :@col-index,
                :@row-ptr,
                nrow => $!nrow + $other.nrow,
                :$!ncol,
                :$!implicit-value
                );
    }

    #=================================================================
    # Column-bind
    #=================================================================
    method column-bind(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        die 'The number of rows of the argument is expected to be equal to the number of rows of the object.'
        unless $!nrow == $other.nrow;

        # Not very effective, but quick to implement
        return self.transpose.row-bind($other.transpose).transpose;
    }

    #=================================================================
    # Slicing
    #=================================================================
    method head(Int $n = 1 --> Math::SparseMatrix::CSR:D) {
        my $end = min($n, $.nrow);
        my @values-slice;
        my @col-index-slice;
        my @row-ptr-slice = @.row-ptr[0 .. $end];

        for ^$end -> $i {
            my $start = @.row-ptr[$i];
            my $stop = @.row-ptr[$i + 1];
            @values-slice.append(@.values[$start ..^ $stop]);
            @col-index-slice.append(@.col-index[$start ..^ $stop]);
        }

        return Math::SparseMatrix::CSR.new(
                :values(@values-slice),
                :col-index(@col-index-slice),
                :row-ptr(@row-ptr-slice),
                :nrow($end),
                :ncol($.ncol),
                :$!implicit-value
                );
    }

    #=================================================================
    # Dense array
    #=================================================================
    #| (Dense) array of arrays representation.
    #| C<:$implicit-value> -- Implicit value to use.
    method Array(:i(:iv(:$implicit-value)) is copy = Whatever) {
        if $implicit-value.isa(Whatever) { $implicit-value = $!implicit-value }
        my @result;
        for ^$!nrow -> $i {
            my @row = ($implicit-value xx $!ncol);
            for @!row-ptr[$i] ..^ @!row-ptr[$i + 1] -> $j {
                @row[@!col-index[$j]] = @!values[$j];
            }
            @result.push(@row);
        }
        return @result;
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose() {
        my @IAT = 0 xx ($!ncol + 1);
        my @JAT;
        my @ANT;

        my $MH = $!ncol + 1;
        my $NH = $!nrow + 1;

        my $IAB = @!row-ptr[$NH - 1];

        for ^$IAB -> $i {
            my $J = @!col-index[$i] + 2;
            if $J < $MH {
                @IAT[$J] += 1;
            }
        }

        @IAT[0] = 0;
        @IAT[1] = 0;

        if $!ncol != 1 {
            for 2 ..^ $MH -> $i {
                @IAT[$i] += @IAT[$i - 1];
            }
        }

        for ^$!nrow -> $i {
            my $IAA = @!row-ptr[$i];
            my $IAB = @!row-ptr[$i + 1];
            next if $IAB < $IAA;
            for $IAA ..^ $IAB -> $jp {
                my $J = @!col-index[$jp] + 1;
                my $K = @IAT[$J];
                @JAT[$K] = $i;
                @ANT[$K] = @!values[$jp];
                @IAT[$J] = $K + 1;
            }
        }

        return Math::SparseMatrix::CSR.new(
                values => @ANT,
                col-index => @JAT,
                row-ptr => @IAT,
                nrow => $!ncol,
                ncol => $!nrow
                );
    }

    #=================================================================
    # Matrix-vector multiplication
    #=================================================================
    #| Dot product of a sparse matrix with dense vector
    multi method dot(@vector --> Array:D) {
        my @result = 0 xx @!row-ptr.elems - 1;
        for @!row-ptr.kv -> $i, $row-start {
            my $row-end = @!row-ptr[$i + 1] // @!values.elems;
            for $row-start ..^ $row-end -> $j {
                @result[$i] += @!values[$j] * @vector[@!col-index[$j]];
            }
        }
        return @result.Array;
    }

    #=================================================================
    # Matrix-matrix multiplication
    #=================================================================
    method dot-pattern(Math::SparseMatrix::CSR $other) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.nrow;

        my @IC = 0 xx ($!nrow + 1);
        my @JC;
        my @IX = 0 xx $other.ncol;
        my $IP = 0;

        for 0 ..^ $other.ncol -> $i {
            @IX[$i] = 0;
        }

        for ^$!nrow -> $i {
            @IC[$i] = $IP;
            my $IAA = @!row-ptr[$i];
            my $IAB = @!row-ptr[$i + 1] - 1;
            if $IAB >= $IAA {
                for $IAA .. $IAB -> $jp {
                    my $j = @!col-index[$jp];
                    my $IBA = $other.row-ptr[$j];
                    my $IBB = $other.row-ptr[$j + 1] - 1;
                    if $IBB >= $IBA {
                        for $IBA .. $IBB -> $kp {
                            my $k = $other.col-index[$kp];
                            if @IX[$k] != $i + 1 {
                                @JC[$IP++] = $k;
                                @IX[$k] = $i + 1;
                            }
                        }
                    }
                }
            }
        }

        @IC[$!nrow] = $IP;
        return Math::SparseMatrix::CSR.new(
                values => (1 xx $IP),
                col-index => @JC,
                row-ptr => @IC,
                :$!nrow,
                ncol => $other.ncol
                );
    }

    method dot-numeric(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $other.nrow;

        my $pattern = self.dot-pattern($other);
        my @IC = $pattern.row-ptr;
        my @JC = $pattern.col-index;
        #my @CN = self.values;
        my @IB = $other.row-ptr;
        my @JB = $other.col-index;
        my @BN = $other.values;
        my @X = 0 xx $other.ncol;
        my @result-values;

        for ^$.nrow -> $i {
            my $ICA = @IC[$i];
            my $ICB = @IC[$i + 1];

            next if $ICB ≤ $ICA;

            for $ICA ..^ $ICB -> $j {
                @X[@JC[$j]] = 0;
            }

            my $IAA = self.row-ptr[$i];
            my $IAB = self.row-ptr[$i + 1];
            for $IAA ..^ $IAB -> $jp {
                my $j = self.col-index[$jp];
                my $a = self.values[$jp];
                my $IBA = @IB[$j];
                my $IBB = @IB[$j + 1];

                next if $IBB ≤ $IBA;

                for $IBA ..^ $IBB -> $kp {
                    my $k = @JB[$kp];
                    @X[$k] += $a * @BN[$kp];
                }
            }

            for $ICA ..^ $ICB -> $j {
                @result-values.push: @X[@JC[$j]];
            }
        }

        Math::SparseMatrix::CSR.new(
                :values(@result-values),
                :col-index(@JC),
                :row-ptr(@IC),
                :nrow($.nrow),
                :ncol($other.ncol)
                );
    }

    #| Dot product of two sparse matrices
    multi method dot(Math::SparseMatrix::CSR:D $B, Bool:D :$sort-indexes = False --> Math::SparseMatrix::CSR:D) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $B.nrow;

        my @values;
        my @col-index;
        my @row-ptr = 0 xx (self.nrow + 1);
        #my @X = 0 xx $B.ncol;

        for ^self.nrow -> $i {
            my %accumulator;
            my $row-start = self.row-ptr[$i];
            my $row-end = self.row-ptr[$i + 1];

            for $row-start ..^ $row-end -> $jp {
                my $j = self.col-index[$jp];
                my $a = self.values[$jp].Num;
                my $col-start = $B.row-ptr[$j];
                my $col-end = $B.row-ptr[$j + 1];

                for $col-start ..^ $col-end -> $kp {
                    my $k = $B.col-index[$kp];
                    %accumulator{$k} += $a * $B.values[$kp];
                }
            }

            if $sort-indexes {
                # Having sorted keys is not essential.
                # It makes the execution of the algorithm 50% slower.
                for %accumulator.keys.sort -> $j {
                    @values.push(%accumulator{$j});
                    @col-index.push($j.Int);
                }
            } else {
                # Why is this so slow?
                # @values = @values.append(%accumulator.values);
                # @col-index = @col-index.append(%accumulator.keys>>.Int);
                for %accumulator.keys -> $j {
                    @values.push(%accumulator{$j});
                    @col-index.push($j.Int);
                }
            }
            @row-ptr[$i + 1] = @values.elems;
        }

        return Math::SparseMatrix::CSR.bless(
                :@values,
                :@col-index,
                :@row-ptr,
                nrow => self.nrow,
                ncol => $B.ncol,
                :$!implicit-value
                );
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition
    multi method add(Numeric:D $a --> Math::SparseMatrix::CSR:D) {
        my @values = @!values >>+>> $a;
        return Math::SparseMatrix::CSR.bless(
                :@values,
                :@!col-index,
                :@!row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $a
                );
    }

    #| Symbolic addition of two matrices
    method add-pattern(Math::SparseMatrix::CSR $other) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        my @IC = 0 xx ($!nrow + 1);
        my @JC;
        my @IX = 0 xx $!ncol;
        my $IP = 0;

        for ^$!nrow -> $i {
            @IC[$i] = $IP;
            my $IAA = @!row-ptr[$i];
            my $IAB = @!row-ptr[$i + 1] - 1;
            if $IAB >= $IAA {
                for $IAA .. $IAB -> $jp {
                    my $j = @!col-index[$jp];
                    @JC[$IP++] = $j;
                    @IX[$j] = $i + 1;
                }
            }

            my $IBA = $other.row-ptr[$i];
            my $IBB = $other.row-ptr[$i + 1] - 1;
            if $IBB >= $IBA {
                for $IBA .. $IBB -> $jp {
                    my $j = $other.col-index[$jp];
                    if @IX[$j] != $i + 1 {
                        @JC[$IP++] = $j;
                    }
                }
            }
        }

        @IC[$!nrow] = $IP;
        return Math::SparseMatrix::CSR.new(
                values => (1 xx $IP),
                col-index => @JC,
                row-ptr => @IC,
                :$!nrow,
                :$!ncol
                );
    }

    #| Numeric addition of two matrices
    multi method add(Math::SparseMatrix::CSR $other --> Math::SparseMatrix::CSR:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;
        if $!implicit-value == 0 && $other.implicit-value == 0 {
            return self!add0($other);
        } else {
            return self!add-iv($other);
        }
    }

    method !add0(Math::SparseMatrix::CSR $other --> Math::SparseMatrix::CSR:D) {
        my $pattern = self.add-pattern($other);

        my @CN = 0 xx $pattern.values.elems;
        # Should be filled-in with the implicit value?
        my @X = 0 xx $pattern.ncol;

        for ^$!nrow -> $i {
            my $IH = $i + 1;
            my $ICA = $pattern.row-ptr[$i];
            my $ICB = $pattern.row-ptr[$IH] - 1;

            next if $ICB < $ICA;

            for $ICA .. $ICB -> $ip {
                @X[$pattern.col-index[$ip]] = 0;
            }

            my $IAA = @!row-ptr[$i];
            my $IAB = @!row-ptr[$IH] - 1;

            if $IAB >= $IAA {
                for $IAA .. $IAB -> $ip {
                    @X[@!col-index[$ip]] = @!values[$ip];
                }
            }

            my $IBA = $other.row-ptr[$i];
            my $IBB = $other.row-ptr[$IH] - 1;

            if $IBB >= $IBA {
                for $IBA .. $IBB -> $ip {
                    my $J = $other.col-index[$ip];
                    @X[$J] += $other.values[$ip];
                }
            }

            for $ICA .. $ICB -> $ip {
                @CN[$ip] = @X[$pattern.col-index[$ip]];
            }
        }

        return Math::SparseMatrix::CSR.new(
                values => @CN,
                col-index => $pattern.col-index,
                row-ptr => $pattern.row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $other.implicit-value
                );
    }

    method !add-iv(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        my @result-values;
        my @result-col-index;
        my @result-row-ptr = 0;

        for ^$!nrow -> $i {
            my %row-a = self!row-map($i);
            my %row-b = $other!row-map($i);
            my @all-columns = [|%row-a.keys, |%row-b.keys].unique;
            my %row-result;

            # This can be optimized, without using so many lookups.
            for @all-columns -> $col {
                if (%row-a{$col}:exists) && (%row-b{$col}:exists) {
                    %row-result{$col} = %row-a{$col} + %row-b{$col};
                } elsif %row-b{$col}:exists {
                    %row-result{$col} = $!implicit-value + %row-b{$col};
                } elsif %row-a{$col}:exists {
                    %row-result{$col} = %row-a{$col} + $other.implicit-value;
                }
            }

            for %row-result.kv -> $col, $val {
                @result-values.push: $val;
                @result-col-index.push: $col;
            }

            @result-row-ptr.push: @result-values.elems;
        }

        return Math::SparseMatrix::CSR.new(
                values    => @result-values,
                col-index => @result-col-index,
                row-ptr   => @result-row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $other.implicit-value
                );
    }

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication
    multi method multiply(Numeric:D $a --> Math::SparseMatrix::CSR:D) {
        my @values = @!values >>*>> $a;
        return Math::SparseMatrix::CSR.bless(
                :@values,
                :@!col-index,
                :@!row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value * $a
                );
    }

    #| Matrix-matrix element-wise multiplication
    multi method multiply(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        die 'The dimensions of the argument must match the dimensions of the object.'
        unless $!nrow == $other.nrow && $!ncol == $other.ncol;

        if $!implicit-value == 0 || $other.implicit-value == 0 {
            return self!multiply0($other)
        } else {
            return self!multiply-iv($other)
        }
    }

    method !multiply0(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        my @result-values;
        my @result-col-index;
        my @result-row-ptr = 0;

        for ^$!nrow -> $i {
            my %row-a = self!row-map($i);
            my %row-b = $other!row-map($i);
            my %row-result;

            for %row-a.kv -> $col, $val-a {
                if %row-b{$col}:exists {
                    %row-result{$col} = $val-a * %row-b{$col};
                }
            }

            for %row-result.kv -> $col, $val {
                @result-values.push: $val;
                @result-col-index.push: $col;
            }

            @result-row-ptr.push: @result-values.elems;
        }

        return Math::SparseMatrix::CSR.new(
                values    => @result-values,
                col-index => @result-col-index,
                row-ptr   => @result-row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => 0 # This method assumes $!implicit-value * $other.implicit-value == 0
                );
    }

    method !multiply-iv(Math::SparseMatrix::CSR:D $other --> Math::SparseMatrix::CSR:D) {
        my @result-values;
        my @result-col-index;
        my @result-row-ptr = 0;

        for ^$!nrow -> $i {
            my %row-a = self!row-map($i);
            my %row-b = $other!row-map($i);
            my @all-columns = [|%row-a.keys, |%row-b.keys].unique;
            my %row-result;

            # This can be optimized, without using so many lookups.
            for @all-columns -> $col {
                if (%row-a{$col}:exists) && (%row-b{$col}:exists) {
                    %row-result{$col} = %row-a{$col} * %row-b{$col};
                } elsif %row-b{$col}:exists {
                    %row-result{$col} = $!implicit-value * %row-b{$col};
                } elsif %row-a{$col}:exists {
                    %row-result{$col} = %row-a{$col} * $other.implicit-value;
                }
            }

            for %row-result.kv -> $col, $val {
                @result-values.push: $val;
                @result-col-index.push: $col;
            }

            @result-row-ptr.push: @result-values.elems;
        }

        return Math::SparseMatrix::CSR.new(
                values    => @result-values,
                col-index => @result-col-index,
                row-ptr   => @result-row-ptr,
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
        @!values = 1 xx (@!values.elems);
        if $!implicit-value { $!implicit-value = 1 }
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
        @!values = @!values.map({
            if $_ < $v-min { $v-min }
            elsif $v-max < $_ { $v-max }
            else {$_}
        });
        # !! What about the implicit value?
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
        @!values = @!values>>.round($scale);
        if $!implicit-value { $!implicit-value .= round($scale) }
        return self;
    }

    #=================================================================
    # Apply elementwise
    #=================================================================
    #| Apply a function to the elements of sparse matrix.
    #| C<&func> -- Function to apply.
    #| C<:$skip-implicit-value> -- Should application to the implicit value be skipped or not?
    #| C<:$clone> -- Whether to operate in-place.
    method apply-elementwise(&func, Bool:D :$skip-implicit-value = False, Bool:D :$clone = True) {
        if $clone {
            return self.clone.apply-elementwise(&func, :$skip-implicit-value, :!clone);
        }
        @!values = @!values.map({ &func($_) });
        if !$skip-implicit-value { $!implicit-value .= &func }
        return self;
    }

    #=================================================================
    # Conjugate
    #=================================================================
    #| Conjugate the sparse matrix
    method conjugate(Bool:D :$clone = True -->Math::SparseMatrix::CSR:D) {
        if $clone {
            return self.clone.conjugate(:!clone);
        }
        @!values = @!values>>.conj;
        if $!implicit-value { $!implicit-value .= conj }
        return self;
    }

    #=================================================================
    # Pretty print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :$echo = True) {
        my $default = $show-implicit-value ?? $!implicit-value.Str !! '.';
        my @rows;
        my $max-len = $default.chars;
        # Minimum length for default / '.'

        for ^self.nrow -> $i {
            my @row = ($default xx self.ncol);
            for self.row-ptr[$i] ..^ self.row-ptr[$i + 1] -> $j {
                @row[self.col-index[$j]] = self.values[$j].Str;
                $max-len = @row[self.col-index[$j]].chars if @row[self.col-index[$j]].chars > $max-len;
            }
            @rows.push(@row);
        }

        @rows = @rows.map({ $_.map({ sprintf("%-*s", $max-len, $_) }) });

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
    #| HTML representation
    method to-html(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False) {
        my $default = $show-implicit-value ?? $!implicit-value.Str !! '.';
        my @rows;
        my $max-len = $default.chars;

        for ^self.nrow -> $i {
            my @row = ($default xx self.ncol);
            for self.row-ptr[$i] ..^ self.row-ptr[$i + 1] -> $j {
                @row[self.col-index[$j]] = self.values[$j].Str;
                $max-len = @row[self.col-index[$j]].chars if @row[self.col-index[$j]].chars > $max-len;
            }
            @rows.push(@row);
        }

        my $html = '<table border="1">';
        for @rows -> @row {
            $html ~= '<tr>' ~ @row.map({ "<td>" ~ $_ ~ "</td>" }).join ~ '</tr>';
        }
        $html ~= '</table>';

        return $html;
    }

    #| Wolfram Language (WL) representation
    method to-wl() {
        my $rules = self.rules.map({ "\{{$_.key.head+1},{$_.key.tail+1}\}->{$_.value}"}).join(',');
        return "SparseArray[\{$rules\}, \{{$!nrow}, {$!ncol}\}, {$!implicit-value}]"
    }


    #| To Hash
    multi method Hash(::?CLASS:D:-->Hash) {
        return
                {
                    specified-elements => self.explicit-length,
                    dimensions => ($!nrow, $!ncol),
                    default => $!implicit-value,
                    density => self.density,
                    elements => self.rules,
                    #                    column-indices => self.column-indices,
                    #                    row-pointers => self.row-pointers,
                    #                    explicit-values => self.explicit-values,
                };
    }

    #| To string
    multi method Str(::?CLASS:D:-->Str) {
        return self.gist;
    }

    #| To gist
    multi method gist(::?CLASS:D:-->Str) {
        return 'Math::SparseMatrix::CSR' ~ (specified-elements => self.explicit-length, dimensions => ($!nrow, $!ncol),
                                            density => self.density).List.raku;
    }
}