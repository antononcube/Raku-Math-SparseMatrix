#unit module Math::SparseMatrix;

use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::DOK;

#=====================================================================
# Math::SparseMatrix
#=====================================================================
class Math::SparseMatrix is Math::SparseMatrix::Abstract {
    has Math::SparseMatrix::Abstract:D $.core-matrix
            is rw
            handles <columns-count explicit-length density dimensions ncol nrow rows-count>
            = Math::SparseMatrix::CSR.new(:0nrow, :0ncol);
    has %.row-names is rw = %();
    has %.column-names is rw = %();
    has %.dimension-names is rw = %();

    #=================================================================
    # Creators
    #=================================================================
    method !process-names($names, $n, $arg-name) {
        return do given $names {
            when ($_ ~~ Seq:D | Range:D) {
                self!process-names($names.Array, $n, $arg-name)
            }
            when ($_ ~~ List:D | Array:D) && $_.unique.elems == $n {
                $_.kv.rotor(2)>>.reverse.flat.Hash
            }
            when ($_ ~~ Map:D) && $_.elems == $n {
                $_
            }
            when Whatever {
                ((^$n) Z=> (^$n)).Hash
            }
            default {
                die "The argument $arg-name is expected to be a Positional or a Map of length $n, or Whatever."
            }
        }
    }

    multi method new(Math::SparseMatrix::Abstract:D :m(:matrix(:$core-matrix)) is copy,
                     :$row-names is copy = Whatever,
                     :$column-names is copy = Whatever,
                     :$dimension-names is copy = Whatever) {
        if $core-matrix ~~ Math::SparseMatrix:D {
            if $row-names.isa(Whatever) { $row-names = $core-matrix.row-names }
            if $column-names.isa(Whatever) { $column-names = $core-matrix.column-names }
            if $dimension-names.isa(Whatever) { $dimension-names = $core-matrix.dimension-names }
            $core-matrix = $core-matrix.core-matrix;
        }
        self.bless(:$core-matrix,
                row-names => self!process-names($row-names, $core-matrix.nrow, 'row-names'),
                column-names => self!process-names($column-names, $core-matrix.ncol, 'column-names'),
                dimension-names => self!process-names($dimension-names, 2, 'dimension-names')
                );
    }

    multi method new(Math::SparseMatrix::Abstract:D $core-matrix,
                     :$row-names = Whatever,
                     :$column-names = Whatever,
                     :$dimension-names = Whatever) {
        self.new(:$core-matrix, :$row-names, :$column-names, :$dimension-names);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.clone,
                row-names => %!row-names.clone,
                column-names => %!column-names.clone,
                dimension-names => %!dimension-names.clone,
                );
    }
    #=================================================================
    # Access
    #=================================================================
    method elems(::?CLASS:D:) {
        return $!core-matrix.nrow;
    }

    multi method value-at(Int:D $row, Int:D $col) {
        return $!core-matrix.value-at($row, $col);
    }

    multi method value-at(Str:D $row, Str:D $col) {
        return $!core-matrix.value-at(%!row-names{$row}, %!column-names{$col});
    }

    multi method row-at(Int:D $row --> Math::SparseMatrix) {
        return $!core-matrix.row-at($row);
    }

    multi method row-at(Str:D $row --> Math::SparseMatrix) {
        return $!core-matrix.row-at(%!row-names{$row});
    }

    method row-slice(*@indexes) {
        die 'The indexes are expected to be non-negative integers or strings.'
        unless (@indexes.all ~~ Int:D | Str:D) && min(@indexes.grep(*~~ Int)) ≥ 0;

        my @indexes2 = @indexes.map({ %!row-names{$_} // $_ });
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.row-slice(@indexes2),
                row-names => %!row-names.grep({ $_.key ∈ @indexes }).Hash,
                :%!column-names
                );
    }

    multi method column-at(Int:D $col --> Math::SparseMatrix) {
        my %columnInds = %!column-names.invert;
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.column-at($col),
                :%!row-names,
                column-names => %(%columnInds{$col} => $col)
                );
    }

    multi method column-at(Str:D $col --> Math::SparseMatrix) {
        die 'Unknown column name.' if %!column-names{$col}:!exists;
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.column-at(%!column-names{$col}),
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
        my $smat = $!core-matrix.transpose;
        return Math::SparseMatrix.new(core-matrix => $smat, row-names => %!column-names, column-names => %!row-names);
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
            $obj.core-matrix = $obj.core-matrix.add($other.core-matrix);
        }
        elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
            $obj.core-matrix = $obj.core-matrix.add($other);
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
            $obj.core-matrix = $obj.core-matrix.multiply($other.core-matrix);
        }
        elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
            $obj.core-matrix = $obj.core-matrix.multiply($other);
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
            $obj.core-matrix = self.core-matrix.dot($other.core-matrix);
            #$obj.core-matrix.eliminate-zeros();
            $obj.column-names = $other.column-names;
            $obj.row-names = self.row-names;
        } elsif $other ~~ Math::SparseMatrix::Abstract:D {
            $obj.core-matrix = self.core-matrix.dot($other);
            #$obj.core-matrix.eliminate-zeros();
            $obj.row-names = self.row-names;
            $obj.column-names = self!process-names(Whatever, $obj.core-matrix.ncol, 'column-names');
        } elsif $other ~~ Seq:D {
            return self.dot($other.Array);
        } elsif $other ~~ Array:D | List:D {
            my @vec = self.core-matrix.dot($other);
            my $res = Math::SparseMatrix::CSR.new(dense-matrix => @vec.map({[$_,]}));
            #$res.eliminate-zeros();
            $obj.core-matrix = $res;
            $obj.row-names = self.row-names;
            $obj.column-names = self!process-names(Whatever, 1, 'column-names');
        } else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Print
    #=================================================================
    method print() {
        my $connector = '┼'; # '+'; # '┼';
        my $v-sep = '│'; #'|'; #'│';
        my $h-sep = '–'; # '-''–'
        my @col-names = self.column-names.keys.sort;
        my @row-names = self.row-names.keys.sort;

        my $col-width = @col-names.map(*.chars).max // 0;
        my $row-width = @row-names.map(*.chars).max // 0;

        my @rows;
        my $max-len = 1;
        # Minimum length for '.'

        @rows = $!core-matrix.print(:!echo);
        $max-len = @rows.map(*.Slip).map(*.chars).max;
        $col-width = max($col-width, $max-len);

        my $header = (' ' x ($row-width + 3)) ~ @col-names.map({ sprintf("%-*s", $max-len, $_) }).join(' ');
        my $header-line1 = $h-sep x $header.chars;
        my $header-line2 = ($h-sep x ($row-width + 1)) ~ $connector ~ ($h-sep x ($header.chars - $row-width - 2));

        say $header-line1;
        say $header;
        say $header-line2;

        for ^@rows.elems -> $i {
            say [sprintf("%-*s", $row-width, @row-names[$i]), $v-sep, |@rows[$i].map({ sprintf("%-*s", $max-len, $_) })].join(' ');
        }
    }


    #=================================================================
    # Representation
    #=================================================================
    method to-html() {
        my $html = '<table border="1">';
        $html ~= '<thead><tr><th></th>';
        for %.column-names.keys.sort -> $col {
            $html ~= "<th>{$col}</th>";
        }
        $html ~= '</tr></thead>';
        for %.row-names.keys.sort -> $row {
            $html ~= "<tr><th>{$row}</th>";
            for %.column-names.keys.sort -> $col {
                $html ~= "<td>{self.value-at(%!row-names{$row}, %!column-names{$col}) // self.implicit-value}</td>";
            }
            $html ~= '</tr>';
        }
        $html ~= '</table>';
        return $html;
    }

    #| Wolfram Language (WL) representation
    method to-wl() {
        my $sp = $!core-matrix.to-wl;
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
                    dimensions => ($!core-matrix.nrow, $!core-matrix.ncol),
                    default => $!core-matrix.implicit-value,
                    :%!row-names,
                    :%!column-names,
                    :%!dimension-names,
                    density => self.density,
                    elements => self.core-matrix.rules,
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
        return 'Math::SparseMatrix' ~ (specified-elements => self.explicit-length,
                                       dimensions => ($!core-matrix.nrow, $!core-matrix.ncol),
                                       density => self.density).List.raku;
    }
}

#=====================================================================
# Abstract
#=====================================================================
multi sub postcircumfix:<[ ]>(Math::SparseMatrix::Abstract:D $mat, *@indexes) is export {
    return $mat.row-slice(@indexes);
}

multi sub postcircumfix:<[; ]>(Math::SparseMatrix::Abstract:D $mat, @indexes) is export {
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
            die "Cannot process the given range."
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
            die "Cannot process the given range."
        }
    }
}
