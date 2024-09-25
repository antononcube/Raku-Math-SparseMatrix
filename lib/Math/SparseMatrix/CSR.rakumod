class Math::SparseMatrix::CSR {
    has @.values;
    has @.col-index;
    has @.row-ptr;
    has UInt:D $.nrow is required;
    has UInt:D $.ncol is required;
    has Numeric:D $.implicit-value = 0;

    #=================================================================
    # Creators
    #=================================================================
    #    submethod BUILD(
    #            :@!values,
    #            :@!col-index,
    #            :@!row-ptr,
    #            UInt:D :$!nrow,
    #            UInt:D :$!ncol
    #                    ) {
    #    }
    #
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
        my %temp;
        for @rules -> $rule {
            my ($row, $col) = $rule.key;
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
                if $val != 0 {
                    @values.push: $val;
                    @col-index.push: $col;
                }
            }
            @row-ptr.push: @values.elems;
        }

        self.bless(:@values, :@col-index, :@row-ptr, :$nrow, :$ncol, :$implicit-value);
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

    #=================================================================
    # Rules and tuples
    #=================================================================
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

    method tuples() {
        return self.rules.map({ [|$_.key, $_.value].List }).Array;
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

    method density() {
        return @!values.elems / ($!nrow * $!ncol);
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
    method Array() {
        my @result;
        for ^$!nrow -> $i {
            my @row = (0 xx $!ncol);
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
    method transpose-by-rules(--> Math::SparseMatrix::CSR) {
        my @rules = |self.rules;
        return Math::SparseMatrix::CSR.new(
                rules => @rules.map({ $_.key.reverse => $_.value }),
                nrow => self.ncol,
                ncol => self.nrow
                );
    }

    method transpose() {
        my @IAT = 0 xx ($!ncol + 1);
        my @JAT;
        my @ANT;

        my $MH = $!ncol + 1;
        my $NH = $!nrow + 1;

        for 2 ..^ $MH -> $i {
            @IAT[$i] = 0;
        }

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
    multi method dot(Math::SparseMatrix::CSR:D $B --> Math::SparseMatrix::CSR:D) {
        die 'The number of rows of the argument is expected to be equal to the number of columns of the object.'
        unless $!ncol == $B.nrow;

        my @values;
        my @col-index;
        my @row-ptr = 0 xx (self.nrow + 1);
        my @X = 0 xx $B.ncol;

        for ^self.nrow -> $i {
            my %accumulator;
            my $row-start = self.row-ptr[$i];
            my $row-end = self.row-ptr[$i + 1];

            for $row-start ..^ $row-end -> $jp {
                my $j = self.col-index[$jp];
                my $a = self.values[$jp];
                my $col-start = $B.row-ptr[$j];
                my $col-end = $B.row-ptr[$j + 1];

                for $col-start ..^ $col-end -> $kp {
                    my $k = $B.col-index[$kp];
                    %accumulator{$k} += $a * $B.values[$kp];
                }
            }

            for %accumulator.keys.sort -> $j {
                @values.push(%accumulator{$j});
                @col-index.push($j.Int);
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
        my $pattern = self.add-pattern($other);

        my @CN = 0 xx $pattern.values.elems;
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

    #=================================================================
    # Mult
    #=================================================================
    multi method mult(Numeric:D $a --> Math::SparseMatrix::CSR:D) {
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

    #=================================================================
    # Pretty print
    #=================================================================
    method print() {
        my @rows;
        my $max-len = 1;
        # Minimum length for '.'

        for ^self.nrow -> $i {
            my @row = ('.' xx self.ncol);
            for self.row-ptr[$i] ..^ self.row-ptr[$i + 1] -> $j {
                @row[self.col-index[$j]] = self.values[$j].Str;
                $max-len = @row[self.col-index[$j]].chars if @row[self.col-index[$j]].chars > $max-len;
            }
            @rows.push(@row);
        }

        for @rows -> @row {
            say @row.map({ sprintf("%-*s", $max-len, $_) }).join(' ');
        }
    }

    #=================================================================
    # Representation
    #=================================================================
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