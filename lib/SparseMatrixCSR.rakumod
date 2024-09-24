class SparseMatrixCSR {
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
    # Rules
    #=================================================================
    method rules() {
        my @rules;
        for @.row-ptr.kv -> $i, $ptr {
            my $next_ptr = $i == @.row-ptr.end ?? @.values.elems !! @.row-ptr[$i + 1];
            for $ptr ..^ $next_ptr -> $j {
                if @.values[$j] != 0 {
                    @rules.push(Pair.new(($i, @.col-index[$j]), @.values[$j]))
                }
            }
        }
        return @rules;
    }

    #=================================================================
    # Slicing
    #=================================================================
    method head(Int $n = 1 --> SparseMatrixCSR:D) {
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

        return SparseMatrixCSR.new(
                :values(@values-slice),
                :col-index(@col-index-slice),
                :row-ptr(@row-ptr-slice),
                :nrow($end),
                :ncol($.ncol)
                );
    }

    method postcircumfix:<[ ]>(Range $range) {
        my $start = $range.min;
        my $end = min($range.max + 1, $.nrow);
        my @values-slice;
        my @col-index-slice;
        my @row-ptr-slice = @.row-ptr[$start .. $end];

        for $start ..^ $end -> $i {
            my $start-idx = @.row-ptr[$i];
            my $stop-idx = @.row-ptr[$i + 1];
            @values-slice.append(@.values[$start-idx ..^ $stop-idx]);
            @col-index-slice.append(@.col-index[$start-idx ..^ $stop-idx]);
        }

        @row-ptr-slice .= map(*- @.row-ptr[$start]);

        return SparseMatrixCSR.new(
                :values(@values-slice),
                :col-index(@col-index-slice),
                :row-ptr(@row-ptr-slice),
                :nrow($end - $start),
                :ncol($.ncol)
                );
    }

    #=================================================================
    # Transpose
    #=================================================================
    method transpose-first(--> SparseMatrixCSR) {
        my @rules = |self.rules;
        return SparseMatrixCSR.new(
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

        my $IAB = @!row-ptr[$NH - 1] - 1;

        for 0 .. $IAB -> $i {
            my $J = @!col-index[$i] + 2;
            if $J <= $MH {
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
            my $IAB = @!row-ptr[$i + 1] - 1;
            if $IAB < $IAA {
                next;
            }
            for $IAA .. $IAB -> $jp {
                my $J = @!col-index[$jp] + 1;
                my $K = @IAT[$J];
                @JAT[$K] = $i;
                @ANT[$K] = @!values[$jp];
                @IAT[$J] = $K + 1;
            }
        }

        return SparseMatrixCSR.new(
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
    multi method dot(SparseMatrixCSR:D $B --> SparseMatrixCSR:D) {
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

        return SparseMatrixCSR.bless(
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
    multi method add(Numeric:D $a --> SparseMatrixCSR:D) {
        my @values = @!values >>+>> $a;
        return SparseMatrixCSR.bless(
                :@values,
                :@!col-index,
                :@!row-ptr,
                :$!nrow,
                :$!ncol,
                :$!implicit-value + $a
                );
    }

    multi method add(SparseMatrixCSR $other --> SparseMatrixCSR:D) {
        my @CN = 0 xx @!values.elems;
        my @X = 0 xx $!ncol;

        for ^$!nrow -> $i {
            my $IH = $i + 1;
            my $ICA = @!row-ptr[$i];
            my $ICB = @!row-ptr[$IH] - 1;

            if $ICB < $ICA {
                next;
            }

            for $ICA .. $ICB -> $ip {
                @X[@!col-index[$ip]] = 0;
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
                @CN[$ip] = @X[@!col-index[$ip]];
            }
        }

        return SparseMatrixCSR.new(
                values => @CN,
                :@!col-index,
                :@!row-ptr,
                :$!nrow,
                :$!ncol,
                implicit-value => $!implicit-value + $other.implicit-value
                );
    }
    #=================================================================
    # Mult
    #=================================================================
    multi method mult(Numeric:D $a --> SparseMatrixCSR:D) {
        my @values = @!values >>*>> $a;
        return SparseMatrixCSR.bless(
                :@values,
                :@!col-index,
                :@!row-ptr,
                :$!nrow,
                :$!ncol,
                :$!implicit-value * $a
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
}