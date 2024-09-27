#unit module Math::SparseMatrix;

use Math::SparseMatrix::CSR;

#=====================================================================
# Math::SparseMatrix
#=====================================================================
class Math::SparseMatrix {
    has Math::SparseMatrix::CSR:D $.sparse-matrix
            is rw
            handles <columns-count explicit-length density dimensions rows-count>
            = Math::SparseMatrix::CSR.new(:0nrow, :0ncol);
    has %.row-names is rw = %();
    has %.column-names is rw = %();
    has %.dimension-names is rw = %();

    #=================================================================
    # Creators
    #=================================================================
    method !process-names($names, $n, $arg-name) {
        return do given $names {
            when ($_ ~~ Seq:D) {
                self!process-names($names.Array, $n, $arg-name)
            }
            when ($_ ~~ List:D | Array:D | Seq:D) && $_.unique.elems == $n {
                $_.kv.rotor(2)>>.reverse.Hash
            }
            when ($_ ~~ Map:D) && $_.elems == $n {
                $_
            }
            when Whatever {
                ((^$n) Z=> (^$n)).Hash
            }
            default {
                die "The argument $arg-name is expected to be a Positional or Map of length $n, or Whatever."
            }
        }
    }

    multi method new(Math::SparseMatrix::CSR:D :m(:matrix(:$sparse-matrix)),
                     :$row-names = Whatever,
                     :$column-names = Whatever,
                     :$dimension-names = Whatever) {
        self.bless(:$sparse-matrix,
                row-names => self!process-names($row-names, $sparse-matrix.nrow, 'row-names'),
                column-names => self!process-names($column-names, $sparse-matrix.ncol, 'column-names'),
                dimension-names => self!process-names($dimension-names, 2, 'dimension-names')
                );
    }

    multi method new(Math::SparseMatrix::CSR:D $sparse-matrix,
                     :$row-names = Whatever,
                     :$column-names = Whatever,
                     :$dimension-names = Whatever) {
        self.new(:$sparse-matrix, :$row-names, :$column-names, :$dimension-names);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix.new(
                sparse-matrix => $!sparse-matrix.clone,
                row-names => %!row-names.clone,
                column-names => %!column-names.clone,
                dimension-names => %!dimension-names.clone,
                );
    }
    #=================================================================
    # Access
    #=================================================================
    method elems(::?CLASS:D:) {
        return $!sparse-matrix.nrow;
    }

    multi method value-at(Int:D $row, Int:D $col) {
        return $!sparse-matrix.value-at($row, $col);
    }

    multi method value-at(Str:D $row, Str:D $col) {
        return $!sparse-matrix.value-at(%!row-names{$row}, %!column-names{$col});
    }

    multi method row-at(Int:D $row --> Math::SparseMatrix) {
        return $!sparse-matrix.row-at($row);
    }

    multi method row-at(Str:D $row --> Math::SparseMatrix) {
        return $!sparse-matrix.row-at(%!row-names{$row});
    }

    method row-slice(*@indexes) {
        die 'The indexes are expected to be non-negative integers or strings.'
        unless (@indexes.all ~~ Int:D | Str:D) && min(@indexes.grep(*~~ Int)) ≥ 0;

        my @indexes2 = @indexes.map({ %!row-names{$_} // $_ });
        return Math::SparseMatrix.new(
                sparse-matrix => $!sparse-matrix.row-slice(@indexes2),
                row-names => %!row-names.grep({ $_.key ∈ @indexes }).Hash,
                :%!column-names
                );
    }

    multi method column-at(Int:D $col --> Math::SparseMatrix) {
        my %columnInds = %!column-names.invert;
        return Math::SparseMatrix.new(
                sparse-matrix => $!sparse-matrix.column-at($col),
                :%!row-names,
                column-names => %(%columnInds{$col} => $col)
                );
    }

    multi method column-at(Str:D $col --> Math::SparseMatrix) {
        die 'Unknown column name.' if %!column-names{$col}:!exists;
        return Math::SparseMatrix.new(
                sparse-matrix => $!sparse-matrix.column-at(%!column-names{$col}),
                :%!row-names,
                column-names => %($col => %!column-names{$col})
                );
    }

    method AT-POS(*@index) {
        if @index.elems == 1 {
            return self.row-at(@index.head);
        }
        die "Only one index is expected.";
    }

    #=================================================================
    # Transpose
    #=================================================================

    #| Transpose the matrix
    method transpose(-->Math::SparseMatrix) {
        my $smat = $!sparse-matrix.transpose;
        return Math::SparseMatrix.new(sparse-matrix => $smat, row-names => %!column-names, column-names => %!row-names);
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition with another Math:SparseMatrix object,
    #| a Math::SparseMatrix::CSR object, or a scalar.
    method add($other, Bool:D $clone = True -->Math::SparseMatrix:D) {
        my $obj = $clone ?? self.clone !! self;
        if ($other ~~ Math::SparseMatrix:D)
                && %!row-names eqv $other.row-names
                && %!column-names eqv $other.column-names {
            $obj.sparse-matrix = $obj.sparse-matrix.add($other.sparse-matrix);
        }
        elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
            $obj.sparse-matrix = $obj.sparse-matrix.add($other);
        }
        else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Multiply
    #=================================================================
    #| Element-wise multiplication with another Math:SparseMatrix object,
    #| a Math::SparseMatrix::CSR object, or a scalar.
    method multiply($other, Bool:D $clone = True -->Math::SparseMatrix:D) {
        my $obj = $clone ?? self.clone !! self;
        if ($other ~~ Math::SparseMatrix:D)
                && %!row-names eqv $other.row-names
                && %!column-names eqv $other.column-names {
            $obj.sparse-matrix = $obj.sparse-matrix.multiply($other.sparse-matrix);
        }
        elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
            $obj.sparse-matrix = $obj.sparse-matrix.multiply($other);
        }
        else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Dot product
    #=================================================================
    method dot($other, Bool :$copy = True) {
        my $obj = $copy ?? Math::SparseMatrix.new() !! self;
        if $other ~~ Math::SparseMatrix:D {
            $obj.sparse-matrix = self.sparse-matrix.dot($other.sparse-matrix);
            #$obj.sparse-matrix.eliminate-zeros();
            $obj.column-names = $other.column-names;
            $obj.row-names = self.row-names;
        } elsif $other ~~ Math::SparseMatrix::CSR:D {
            $obj.sparse-matrix = self.sparse-matrix.dot($other);
            #$obj.sparse-matrix.eliminate-zeros();
            $obj.row-names = self.row-names;
        } elsif $other ~~ Seq:D {
            return self.dot($other.Array);
        } elsif $other ~~ Array:D | List:D {
            my @vec = self.sparse-matrix.dot($other);
            my $res = Math::SparseMatrix::CSR.new(dense-matrix => @vec.map({[$_,]}));
            #$res.eliminate-zeros();
            $obj.sparse-matrix = $res;
            $obj.row-names = self.row-names;
        } else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Print
    #=================================================================
    method print() {
        my @col-names = self.column-names.keys.sort;
        my @row-names = self.row-names.keys.sort;

        my $col-width = @col-names.map(*.chars).max // 0;
        my $row-width = @row-names.map(*.chars).max // 0;

        my @rows;
        my $max-len = 1;
        # Minimum length for '.'

        for ^self.sparse-matrix.nrow -> $i {
            my @row = ('.' xx self.sparse-matrix.ncol);
            for self.sparse-matrix.row-ptr[$i] ..^ self.sparse-matrix.row-ptr[$i + 1] -> $j {
                @row[self.sparse-matrix.col-index[$j]] = self.sparse-matrix.values[$j].Str;
                if @row[self.sparse-matrix.col-index[$j]].chars > $max-len {
                    $max-len = @row[self.sparse-matrix.col-index[$j]].chars
                }
            }
            @rows.push(@row);
        }

        $col-width = max($col-width, $max-len);

        my $header = (' ' x ($row-width + 3)) ~ @col-names.map({ sprintf("%-*s", $max-len, $_) }).join(' ');
        my $header-line1 = '=' x $header.chars;
        my $header-line2 = '-' x $header.chars;

        say $header-line1;
        say $header;
        say $header-line2;

        for ^@rows.elems -> $i {
            say [sprintf("%-*s", $row-width, @row-names[$i]), '|', |@rows[$i].map({ sprintf("%-*s", $max-len, $_) })]
                    .join(' ');
        }
    }


    #=================================================================
    # Representation
    #=================================================================
    #| Wolfram Language (WL) representation
    method wl() {
        my $rules = $!sparse-matrix.rules.map({ "\{{$_.key.head+1},{$_.key.tail+1}\}->{$_.value}"}).join(',');
        my $sp = "SparseArray[\{$rules\}, \{{$!sparse-matrix.nrow}, {$!sparse-matrix.ncol}\}, {$!sparse-matrix.implicit-value}]";
        my @row-names-list = %!row-names.pairs.sort({ $_.value })>>.key;
        my @column-names-list = %!row-names.pairs.sort({ $_.value })>>.key;
        my $rowNames = @row-names-list.raku.trans('[]'=>'{}');
        my $colNames = @column-names-list.raku.trans('[]'=>'{}');
        return "ToSSparseMatrix[$sp, \"RowNames\" -> $rowNames, \"ColumnNames\" -> $colNames]";
    }


    #| To Hash
    multi method Hash(::?CLASS:D:-->Hash) {
        return
                {
                    specified-elements => self.explicit-length,
                    dimensions => ($!sparse-matrix.nrow, $!sparse-matrix.ncol),
                    default => $!sparse-matrix.implicit-value,
                    :%!row-names,
                    :%!column-names,
                    :%!dimension-names,
                    density => self.density,
                    elements => self.sparse-matrix.rules,
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
        return 'Math::SparseMatrix' ~ (specified-elements => self.explicit-length, dimensions => ($!sparse-matrix.nrow, $!sparse-matrix.ncol),
                                            density => self.density).List.raku;
    }
}

#=====================================================================
#proto sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {*}
multi sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {
    return $mat.row-slice(@indexes);
}

multi sub postcircumfix:<[; ]>(Math::SparseMatrix::CSR:D $mat, @indexes) is export {
    return do given (@indexes[0], @indexes[1]) {
        when ($_.head ~~ Int || $_.tail ~~ Str) && ($_.tail ~~ Int || $_.tail ~~ Str) {
            $mat.value-at(@indexes[0], @indexes[1]);
        }
        when $_.head ~~ Range && ($_.tail ~~ Int || $_.tail ~~ Str) {
            $mat.row-slice($_.head).column-at($_.tail)
        }
        when $_.head.isa(Whatever) && ($_.tail ~~ Int || $_.tail ~~ Str) {
            $mat.column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Range {
            $mat.transpose.row-slice($_.tail).transpose
        }
        when $_.head ~~ Range && $_.tail ~~ Range {
            $mat.row-slice($_.head).transpose.row-slice($_.tail).transpose
        }
        default {
            die "Cannot process the given range"
        }
    }
}

#=====================================================================
#proto sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {*}
multi sub postcircumfix:<[ ]>(Math::SparseMatrix:D $mat, *@indexes) is export {
    return $mat.row-slice(@indexes);
}

multi sub postcircumfix:<[; ]>(Math::SparseMatrix:D $mat, @indexes) is export {
    return do given (@indexes[0], @indexes[1]) {
        when $_.head ~~ Int && $_.tail ~~ Int {
            $mat.value-at(@indexes[0], @indexes[1]);
        }
        when $_.head ~~ Range && ($_.tail ~~ Int || $_.tail ~~ Str) {
            $mat.row-slice($_.head).column-at($_.tail)
        }
        when $_.head.isa(Whatever) && ($_.tail ~~ Int || $_.tail ~~ Str) {
            $mat.column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Range {
            $mat.transpose.row-slice($_.tail).transpose
        }
        when $_.head ~~ Range && $_.tail ~~ Range {
            $mat.row-slice($_.head).transpose.row-slice($_.tail).transpose
        }
        default {
            die "Cannot process the given range"
        }
    }
}
