unit module Math::SparseMatrix;

use Math::SparseMatrix::CSR;

#=====================================================================
# Math::SparseMatrix
#=====================================================================
my class Math::SparseMatrix {
    has Math::SparseMatrix::CSR:D $.sparse-matrix is required;
    has %.row-names is rw = %();
    has %.column-names is rw = %();
    has %.dim-names is rw = %();

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
        unless (@indexes.all ~~ Int:D | Str:D) && min(@indexes.grep(* ~~ Int)) ≥ 0;

        my @indexes2 = @indexes.map({ %!row-names{$_} // $_ });
        return $!sparse-matrix.row-slice(@indexes2);
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
                $max-len = @row[self.sparse-matrix.col-index[$j]].chars if @row[self.sparse-matrix.col-index[$j]]
                        .chars > $max-len;
            }
            @rows.push(@row);
        }

        $col-width = max($col-width, $max-len);

        my $header = (' ' x ($row-width + 3) ) ~ @col-names.map({ sprintf("%-*s", $max-len, $_) }).join(' ');
        my $header-line1 = '=' x $header.chars;
        my $header-line2 = '-' x $header.chars;

        say $header-line1;
        say $header;
        say $header-line2;

        for ^@rows.elems -> $i {
            say [sprintf("%-*s", $row-width, @row-names[$i]), '|', |@rows[$i].map({ sprintf("%-*s", $max-len, $_) })].join(' ');
        }
    }
}

#=====================================================================
#proto sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {*}
multi sub postcircumfix:<[ ]>(Math::SparseMatrix::CSR:D $mat, *@indexes) is export {
    return $mat.row-slice(@indexes);
}

multi sub postcircumfix:<[; ]>(Math::SparseMatrix::CSR:D $mat, @indexes) is export {
    return do given (@indexes[0], @indexes[1]) {
        when $_.head ~~ Int && $_.tail ~~ Int {
            $mat.value-at(@indexes[0], @indexes[1]);
        }
        when $_.head ~~ Range && $_.tail ~~ Int {
            $mat.row-slice($_.head).column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Int {
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
        when $_.head ~~ Range && $_.tail ~~ Int {
            $mat.row-slice($_.head).column-at($_.tail)
        }
        when $_.head.isa(Whatever) && $_.tail ~~ Int {
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
