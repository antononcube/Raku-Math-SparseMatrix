#unit module Math::SparseMatrix;

use Math::SparseMatrix::Abstract;
use Math::SparseMatrix::CSR;
use Math::SparseMatrix::Convertish;
use Math::SparseMatrix::DOK;

#=====================================================================
# Math::SparseMatrix
#=====================================================================
class Math::SparseMatrix
        does Math::SparseMatrix::Convertish
        is Math::SparseMatrix::Abstract {
    has Math::SparseMatrix::Abstract:D $.core-matrix
            is rw
            handles <columns-count explicit-length density dimensions implicit-value ncol nrow rows-count Array>
            = Math::SparseMatrix::CSR.new(:0nrow, :0ncol);
    has %.row-names-map;
    has %.column-names-map;
    has %.dimension-names-map;

    has @.row-names;
    has @.column-names;
    has @.dimension-names;

    method set-row-names(@names) {
        my @res = self!process-names(@names.Array, $!core-matrix.nrow, 'row-names');
        @!row-names = |@res.head;
        %!row-names-map = @res.tail;
        return self;
    }

    method set-column-names(@names) {
        my @res = self!process-names(@names.Array, $!core-matrix.ncol, 'column-names');
        @!column-names = |@res.head;
        %!column-names-map = @res.tail;
        return self;
    }

    method set-dimension-names(@names) {
        my @res = self!process-names(@names.Array, $!core-matrix.ncol, 'dimension-names');
        @!dimension-names = |@res.head;
        %!dimension-names-map = @res.tail;
        return self;
    }

    method get-row-name(Int:D $i) {
        return @!row-names[$i];
    }

    method get-column-name(Int:D $i) {
        return @!column-names[$i];
    }

    method get-dimension-name(Int:D $i) {
        return @!dimension-names[$i];
    }

    #=================================================================
    # Creators
    #=================================================================
    # Process various names specs into a list of two elements: (1) names, and (2) name-to-index hashmap.
    method !process-names($names, $n, $arg-name) {
        return do given $names {
            when ($_ ~~ Seq:D | Range:D) {
                self!process-names($names.Array, $n, $arg-name)
            }
            when ($_ ~~ List:D | Array:D) && $_.unique.elems == $n {
                ($names, $names.kv.rotor(2)>>.reverse.flat.Hash)
            }
            when ($_ ~~ Map:D) && $_.elems == $n {
                my @names = $_.sort(*.value)>>.key;
                (@names, $_)
            }
            when ($_ ~~ Map:D | Array:D | List:D) && $_.elems == 0 {
                # This is somewhat hacky.
                # It might happen with a step-by-step building of a Math::SparseMatrix object.
                (Empty, $_)
            }
            when Whatever {
                ((^$n)>>.Str, ((^$n) Z=> (^$n)).Hash)
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
        my @rowNames = self!process-names($row-names, $core-matrix.nrow, 'row-names');
        my @colNames = self!process-names($column-names, $core-matrix.ncol, 'column-names');
        my @dimNames = self!process-names($dimension-names, 2, 'dimension-names');

        self.bless(:$core-matrix,
                row-names => @rowNames.head.Array,
                row-names-map => @rowNames.tail,
                column-names => @colNames.head.Array,
                column-names-map => @colNames.tail,
                dimension-names => @dimNames.head.Array,
                dimension-names-map => @dimNames.tail,
                );
    }

    multi method new(Math::SparseMatrix::Abstract:D $core-matrix,
                     :$row-names = Whatever,
                     :$column-names = Whatever,
                     :$dimension-names = Whatever) {
        self.new(:$core-matrix, :$row-names, :$column-names, :$dimension-names);
    }

    multi method new(:@rules! where @rules.all ~~ Pair:D,
                     :$nrow is copy = @rules.map(*.key[0]).max + 1,
                     :$ncol is copy = @rules.map(*.key[1]).max + 1,
                     :$row-names is copy = Whatever,
                     :$column-names is copy = Whatever,
                     :$dimension-names = Whatever,
                     Numeric:D :$implicit-value = 0) {
        my $core-matrix = Math::SparseMatrix::CSR.new(:@rules, :$nrow, :$ncol, :$implicit-value);
        self.new(:$core-matrix, :$row-names, :$column-names, :$dimension-names);
     }

    multi method new(:@dense-matrix! where @dense-matrix ~~ (Array:D | List:D | Seq:D) && @dense-matrix.all ~~ (Array:D | List:D | Seq:D),
                     :$nrow is copy = @dense-matrix.elems,
                     :$ncol is copy = @dense-matrix>>.elems.max,
                     :$row-names is copy = Whatever,
                     :$column-names is copy = Whatever,
                     :$dimension-names = Whatever,
                     Numeric:D :$implicit-value = 0) {
        my $core-matrix = Math::SparseMatrix::CSR.new(:@dense-matrix, :$nrow, :$ncol, :$implicit-value);
        self.new(:$core-matrix, :$row-names, :$column-names, :$dimension-names);
    }

    multi method new(:@edge-dataset! where @edge-dataset.all ~~ Map:D,
                     Bool:D :$directed = False,
                     Numeric:D :$implicit-value = 0,
                     :$row-names is copy = Whatever,
                     :$column-names is copy = Whatever,
                     :$dimension-names is copy = Whatever) {
        if $row-names ~~ Seq:D { $row-names = $row-names.Array }
        if $column-names ~~ Seq:D { $column-names = $column-names.Array }

        die "Row names is expected to be a list of strings or Whatever."
        unless $row-names.isa(Whatever) || ($row-names ~~ Positional:D | Seq:D) && ($row-names.all ~~ Str:D);

        die "Column names is expected to be a list of strings or Whatever."
        unless $column-names.isa(Whatever) || ($column-names ~~ Positional:D | Seq:D) && ($column-names.all ~~ Str:D);

        my @rules;
        if $directed {
            # Row names
            my @rowNames = @edge-dataset.map(*<from>).unique.sort;
            if $row-names.isa(Whatever) {
                $row-names = @rowNames;
            } else {
                my $set = Set.new($row-names);
                @rowNames = @rowNames.grep({ $_ ∈ $set })
            }
            my %rowNames = $row-names.Array Z=> (^$row-names.elems);

            # Column names
            my @colNames = @edge-dataset.map(*<to>).unique.sort;
            if $column-names.isa(Whatever) {
                $column-names = @colNames;
            } else {
                my $set = Set.new($column-names);
                @colNames = @colNames.grep({ $_ ∈ $set })
            }
            my %colNames = $column-names.Array Z=> (^$column-names.elems);

            # Rules
            @rules = @edge-dataset.map({
                if (%rowNames{$_<from>}:exists) && (%colNames{$_<to>}:exists) {
                    (%rowNames{$_<from>}, %colNames{$_<to>}) => $_<weight>
                } else {
                    Empty
                }
            });
        } else {

            my @rowNames = [|@edge-dataset.map(*<from>), |@edge-dataset.map(*<to>)].unique.sort;

            given ($row-names, $column-names) {
                when (Whatever, Whatever) {
                    $row-names = @rowNames
                }
                when (Positional, Whatever) {
                }
                when (Whatever, Positional) {
                    $row-names = $column-names
                }
                when (Positional, Positional) {
                    die 'If the $directed is set to False, then $row-names and $column-names are expected to match if they are both lists.'
                    unless $row-names eqv $column-names
                }
                default {
                    # Should not happen
                    die 'If the $directed is set to False, then $row-names and $column-names are expected to be compatible.'
                }
            }
            $column-names = $row-names;
            my %allNames = $row-names.Array Z=> (^$row-names.elems);

            # Rules
            @rules = @edge-dataset.map({
                if (%allNames{$_<from>}:exists) && (%allNames{$_<to>}:exists) {
                    [
                        (%allNames{$_<from>}, %allNames{$_<to>}) => $_<weight>,
                        (%allNames{$_<to>}, %allNames{$_<from>}) => $_<weight>,
                    ]
                } else {
                    Empty
                }
            }).map(*.Slip);
        }

        # Core matrix
        my $core-matrix = Math::SparseMatrix::CSR.new(:@rules, nrow => $row-names.elems, ncol => $column-names.elems, :$implicit-value);

        # Object
        self.new(:$core-matrix, :$row-names, :$column-names, :$dimension-names);
    }

    #=================================================================
    # Clone
    #=================================================================
    method clone() {
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.clone,
                row-names => @!row-names.clone,
                column-names => @!column-names.clone,
                dimension-names => @!dimension-names.clone,
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
        die "Unknown row name: $row" unless %!row-names-map{$row}.defined;
        die "Unknown column name: $col" unless %!column-names-map{$col}.defined;
        return $!core-matrix.value-at(%!row-names-map{$row}, %!column-names-map{$col});
    }

    multi method value-at(Str:D $row, Int:D $col) {
        die "Unknown row name: $row" unless %!row-names-map{$row}.defined;
        return $!core-matrix.value-at(%!row-names-map{$row}, $col);
    }

    multi method value-at(Int:D $row, Str:D $col) {
        die "Unknown column name: $col" unless %!column-names-map{$col}.defined;
        return $!core-matrix.value-at($row, %!column-names-map{$col});
    }

    multi method row-at(Int:D $row --> Math::SparseMatrix) {
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.row-at($row),
                row-names => [self.get-row-name($row), ],
                :@!column-names
                );
    }

    multi method row-at(Str:D $row --> Math::SparseMatrix) {
        return self.row-at(%!row-names-map{$row});
    }

    method row-slice(*@indexes) {
        die 'The indexes are expected to be non-negative integers or strings.'
        unless (@indexes.all ~~ Int:D | Str:D) && min(@indexes.grep(*~~ Int)) ≥ 0;

        my @indexes2 = @indexes.map({ %!row-names-map{$_} // $_ });
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.row-slice(@indexes2),
                row-names => @indexes,
                :@!column-names,
                );
    }

    multi method column-at(Int:D $col --> Math::SparseMatrix) {
        my %columnInds = %!column-names-map.invert;
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.column-at($col),
                row-names => %!row-names-map,
                column-names => [%columnInds{$col},]
                );
    }

    multi method column-at(Str:D $col --> Math::SparseMatrix) {
        die 'Unknown column name.' if %!column-names-map{$col}:!exists;
        return Math::SparseMatrix.new(
                core-matrix => $!core-matrix.column-at(%!column-names-map{$col}),
                :%!row-names-map,
                column-names => [$col,]
                );
    }

    method AT-POS($_) {
        return self.row-at($_);
    }

    #=================================================================
    # Impose row names and column names
    #=================================================================
    #| Impose row names. (New Math::SparseMatrix object is created.)
    method impose-row-names(@names, Bool:D :$clone = True) {
        my $obj = $clone ?? self.clone !! self;

        die "The first argument is expected to be a list of strings."
        unless @names ~~ (Array:D | List:D | Seq:D) && @names.all ~~ Str:D;

        my @missing-rows = (@names (-) $obj.row-names).keys;
        my $n-missing-rows = @missing-rows.elems;

        if $n-missing-rows > 0 {
            my $compl-mat = Math::SparseMatrix.new(
                    rules => [ ],
                    nrow => $n-missing-rows,
                    ncol => $obj.ncol,
                    row-names => @missing-rows,
                    column-names => $obj.column-names
                    );

            $obj = $obj.row-bind($compl-mat);
        }

        return $obj[@names; *];
    }

    #| Impose column names. (New Math::SparseMatrix object is created.)
    method impose-column-names(@names, Bool:D :$clone = True) {
        # Not effective, but very quick to implement.
        # The method .transpose() always clones,
        # hence, this should be re-written to be same/similar to .impose-row-names().
        return self.transpose().impose-row-names(@names, :!clone).transpose();
    }

    #=================================================================
    # Row-bind
    #=================================================================
    method row-bind(Math::SparseMatrix:D $other) {
        die "The column names of the two Math::SparseMatrix objects are expected to be the same."
        unless self.column-names.sort eqv $other.column-names.sort;

        my $res-matrix;

        if self.column-names eqv $other.column-names {
            $res-matrix = $!core-matrix.row-bind($other.core-matrix);
        } else {
            my $other-matrix = $other[*;self.column-names].core-matrix;
            $res-matrix = $!core-matrix.row-bind($other-matrix);
        }

        my %newRowMap = %!row-names-map , $other.row-names-map;
        my @newRowNames = do if %newRowMap.elems == (self.nrow + $other.nrow) {
            self.row-names.clone.append($other.row-names)
        } else {
            self.row-names.clone.map({ $_ ~ ".1" }).Array.append($other.row-names.map({ $_ ~ ".2" }).Array)
        }

        return Math::SparseMatrix.new(core-matrix => $res-matrix, :@!column-names, row-names => @newRowNames);
    }

    #=================================================================
    # Column-bind
    #=================================================================
    method column-bind(Math::SparseMatrix:D $other) {
        die "The row names of the two Math::SparseMatrix objects are expected to be the same."
        unless self.row-names.sort eqv $other.row-names.sort;

        my $res-matrix;

        if self.row-names eqv $other.row-names {
            $res-matrix = $!core-matrix.column-bind($other.core-matrix);
        } else {
            my $other-matrix = $other[self.row-names].core-matrix;
            $res-matrix = $!core-matrix.column-bind($other-matrix);
        }

        my %newColumnMap = %!column-names-map , $other.column-names-map;
        my @newColumnNames = do if %newColumnMap.elems == (self.ncol + $other.ncol) {
            self.column-names.clone.append($other.column-names)
        } else {
            self.column-names.clone.map({ $_ ~ ".1" }).Array.append($other.column-names.map({ $_ ~ ".2" }).Array)
        }

        return Math::SparseMatrix.new(core-matrix => $res-matrix, :@!row-names, column-names => @newColumnNames);
    }

    #=================================================================
    # Rules and tuples
    #=================================================================

    method rules(Bool:D :n(:$names) = False) {
        my @res = do if $names {
            self.core-matrix.rules.map({ (self.row-names[$_.key[0]], self.column-names[$_.key[1]]) => $_.value }).Array
        } else {
            self.core-matrix.rules
        }

        return @res;
    }

    method tuples(Bool:D :n(:$names) = False, Bool:D :d(:$dataset) = False) {
        if !$names {
            return self.core-matrix.tuples(:$dataset)
        }

        my @res = do if $dataset {
            self.core-matrix.tuples(:dataset).map({ %(i => self.row-names[$_<i>], j => self.column-names[$_<j>], x => $_<x>) });
        } else {
            self.core-matrix.tuples.map({ (self.row-names[$_[0]], self.column-names[$_[1]], $_[2]) });
        }
        return @res;
    }

    #=================================================================
    # Transpose
    #=================================================================

    #| Transpose the matrix
    method transpose(-->Math::SparseMatrix) {
        my $smat = $!core-matrix.transpose;
        return Math::SparseMatrix.new(core-matrix => $smat, row-names => %!column-names-map, column-names => %!row-names-map);
    }

    #=================================================================
    # Add
    #=================================================================
    #| Element-wise addition with another Math:SparseMatrix object,
    #| a Math::SparseMatrix::CSR object, or a scalar.
    method add($other, Bool:D $clone = True -->Math::SparseMatrix:D) {
        my $obj = $clone ?? self.clone !! self;
        if ($other ~~ Math::SparseMatrix:D)
                && @!row-names eqv $other.row-names
                && @!column-names eqv $other.column-names {
            $obj.core-matrix = $obj.core-matrix.add($other.core-matrix);
        } elsif ($other ~~ Math::SparseMatrix:D)
                && $other.row-names eqv (^$other.nrow)».Str.Array
                && $other.column-names eqv (^$other.ncol)».Str.Array {
            $obj.core-matrix = $obj.core-matrix.add($other.core-matrix);
        } elsif ($other ~~ Math::SparseMatrix:D)
                && self.nrow eq $other.nrow
                && self.ncol eq $other.ncol {
            $obj = $obj.core-matrix.add($other.core-matrix);
            $obj = Math::SparseMatrix.new($obj);
        } elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
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
                && @!row-names eqv $other.row-names
                && @!column-names eqv $other.column-names {
            $obj.core-matrix = $obj.core-matrix.multiply($other.core-matrix);
        } elsif ($other ~~ Math::SparseMatrix:D)
                && $other.row-names eqv (^$other.nrow)».Str.Array
                && $other.column-names eqv (^$other.ncol)».Str.Array {
            $obj.core-matrix = $obj.core-matrix.multiply($other.core-matrix);
        } elsif ($other ~~ Math::SparseMatrix:D)
                && self.nrow eq $other.nrow
                && self.ncol eq $other.ncol {
            $obj = $obj.core-matrix.multiply($other.core-matrix);
            $obj = Math::SparseMatrix.new($obj);
        } elsif $other ~~ Numeric:D || $other ~~ Math::SparseMatrix::CSR {
            $obj.core-matrix = $obj.core-matrix.multiply($other);
        }
        else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Negate
    #=================================================================
    #| Negate the elements of a matrix
    method negate(Bool:D $clone = True -->Math::SparseMatrix:D) {
        # Maybe redundant, but useful for "symmetry" with the unary operation, e.g. -$m .
        return self.multiply(-1, :$clone);
    }

    #=================================================================
    # Dot product
    #=================================================================
    method dot($other, Bool :$copy = True) {
        my $obj = $copy ?? Math::SparseMatrix.new() !! self;
        if $other ~~ Math::SparseMatrix:D {
            # We have to make sure the row names and column names match!
            if self.column-names eqv $other.row-names {
                # Optimization
                $obj.core-matrix = self.core-matrix.dot($other.core-matrix);
            } elsif $other.row-names eqv (^$other.nrow)».Str.Array {
                $obj.core-matrix = self.core-matrix.dot($other.core-matrix);
            } elsif self.column-names eqv (^self.ncol)».Str.Array {
                $obj.core-matrix = self.core-matrix.dot($other.core-matrix);
            } else {
                $obj.core-matrix = self.core-matrix.dot($other.row-slice(self.column-names).core-matrix);
            }
            #$obj.core-matrix.eliminate-zeros();
            $obj.column-names = $other.column-names;
            $obj.column-names-map = $other.column-names-map;
            $obj.row-names = self.row-names;
            $obj.row-names-map = self.row-names-map;
        } elsif $other ~~ Math::SparseMatrix::Abstract:D {
            $obj.core-matrix = self.core-matrix.dot($other);
            #$obj.core-matrix.eliminate-zeros();
            $obj.row-names = self.row-names;
            $obj.row-names-map = self.row-names-map;
            my $colNames = self!process-names(Whatever, $obj.core-matrix.ncol, 'column-names');
            $obj.column-names = $colNames.head.Array;
            $obj.column-names-map = $colNames.tail;
        } elsif $other ~~ Seq:D {
            return self.dot($other.Array);
        } elsif $other ~~ Array:D | List:D {
            my @vec = self.core-matrix.dot($other);
            my $res = Math::SparseMatrix::CSR.new(dense-matrix => @vec.map({[$_,]}));
            #$res.eliminate-zeros();
            $obj.core-matrix = $res;
            $obj.row-names = self.row-names;
            $obj.row-names-map = self.row-names-map;
            my $colNames = self!process-names(Whatever, 1, 'column-names');
            $obj.column-names = $colNames.head.Array;
            $obj.column-names-map = $colNames.tail;
        } else {
            die "The first argument is expected to be a number, a Math::SparseMatrix object, or a Math::SparseMatrix::CSR object.";
        }
        return $obj;
    }

    #=================================================================
    # Row and column sums and maxes
    #=================================================================
    method !row-op(&op, Bool:D :p(:$pairs) = False) {
        # Probably more effective implementations can be provided by the core matrix classes.
        # But this is universal and quick to implement.

        my @sums = do for ^self.nrow -> $i {
            self.row-at($i).tuples.map(*[2]).&op
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

    method column-sums(Bool:D :p(:$pairs) = False) {
        # Quick to implement.
        return self.transpose.row-sums(:$pairs);
    }

    #| Row sums for a sparse matrix.
    method row-maxes(Bool:D :p(:$pairs) = False) {
        my &amax = -> @a { @a.Array.push(self.implicit-value).max };
        self!row-op(&amax, :$pairs);
    }

    method column-maxes(Bool:D :p(:$pairs) = False) {
        # Quick to implement.
        return self.transpose.row-maxes(:$pairs);
    }

    #=================================================================
    # Print
    #=================================================================
    method print(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, ) {
        my $connector = '┼'; # '+'; # '┼';
        my $v-sep = '│'; #'|'; #'│';
        my $h-sep = '–'; # '-''–'
        my @col-names = self.column-names;
        my @row-names = self.row-names;

        my $col-width = @col-names ?? @col-names.map(*.chars).max !! 0;
        my $row-width = @row-names ?? @row-names.map(*.chars).max !! 0;

        my @rows;
        # Minimum length for '.' is 1
        my $max-len = max($col-width, 1);

        @rows = $!core-matrix.print(:!echo, :$show-implicit-value);
        $max-len = max($max-len, @rows.map(*.Slip).map(*.chars).max);
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
    # Wrapper delegations
    #=================================================================
    #| Round the sparse matrix
    #| C<$scale> -- Scale to round to.
    #| C<:$clone> -- Whether to clone or not.
    method round(Numeric:D $scale = 1, Bool:D :$clone = True) {
        if $clone {
            return self.clone.round($scale, :!clone);
        }
        $!core-matrix.round($scale, :!clone);
        return self;
    }

    #| Conjugate the sparse matrix
    #| C<:$clone> -- Whether to clone or not.
    method conjugate(Bool:D :$clone = True) {
        if $clone {
            return self.clone.conjugate(:!clone);
        }
        $!core-matrix.conjugate(:!clone);
        return self;
    }

    #=================================================================
    # Representation
    #=================================================================
    #| HTML representation
    method to-html(Bool:D :iv(:implicit-value(:$show-implicit-value)) = False, Bool:D :v(:$vertical-column-names) = False) {
        my $default = $show-implicit-value ?? self.implicit-value.Str !! '.';
        my @indPairs = self.rules».key.map(*.Str);
        my $html = '<table border="1">';
        $html ~= '<thead><tr><th></th>';
        my $thLeft = $vertical-column-names ?? '<th style="writing-mode: vertical-rl; white-space: nowrap; vertical-align: bottom;">' !! '<th>';
        for @!column-names -> $col {
            $html ~= "{$thLeft}{$col}</th>";
        }
        $html ~= '</tr></thead>';
        for @!row-names -> $row {
            $html ~= "<tr><th>{$row}</th>";
            for @!column-names -> $col {
                my $value = do if (%!row-names-map{$row}, %!column-names-map{$col}).Str ∈ @indPairs {
                    self.value-at(%!row-names-map{$row}, %!column-names-map{$col});
                } else {
                    $default
                }
                $html ~= "<td>{$value}</td>";
            }
            $html ~= '</tr>';
        }
        $html ~= '</table>';
        return $html;
    }

    #| Wolfram Language (WL) representation
    method to-wl() {
        my $sp = $!core-matrix.to-wl;
        my @row-names-list = %!row-names-map.pairs.sort({ $_.value })>>.key;
        my @column-names-list = %!row-names-map.pairs.sort({ $_.value })>>.key;
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
                    :@!row-names,
                    :@!column-names,
                    :@!dimension-names,
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
#proto sub postcircumfix:<[ ]>(Math::SparseMatrix:D $mat, *@indexes) is export {*}
#multi sub postcircumfix:<[ ]>(Math::SparseMatrix:D $mat, *@indexes) is export {
#    note "HERE 1";
#    return $mat.row-slice(@indexes);
#}

multi sub postcircumfix:<[ ]>(Math::SparseMatrix:D $mat, $index) is export {
    given $index {
        when $_ ~~ Str | Int { return $mat.AT-POS($index) }
        when Iterable:D { $mat.row-slice(|$index) }
        default { die 'Cannot handle position spec.' }
    }
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
        when ($_.head ~~ Range || $_.head ~~ Iterable) && $_.tail.isa(Whatever) {
            $mat.row-slice($_.head);
        }
        when $_.head.isa(Whatever) && ($_.tail ~~ Range || $_.tail ~~ Iterable) {
            $mat.transpose.row-slice($_.tail).transpose
        }
        when ($_.head ~~ Range || $_.head ~~ Iterable) && ($_.tail ~~ Range || $_.tail ~~ Iterable) {
            $mat.row-slice($_.head).transpose.row-slice($_.tail).transpose
        }
        default {
            die "Cannot process the given range."
        }
    }
}

#=====================================================================
# Unary operator shortcuts
# These operators <@ %> are found "Math::Matrix" but they break typical code like:
#   |%opts
# multi sub prefix:<@>( Math::SparseMatrix:D $m --> Array) is export(:ALL) { $m.Array }

#
# multi sub prefix:<%>( Math::SparseMatrix:D $m --> Hash) is export(:ALL) { $m.rules.Hash }

#=====================================================================
# Arithmetic
multi sub prefix:<->( Math::SparseMatrix:D $m --> Math::SparseMatrix:D) is export(:ALL) { $m.multiply(-1) }

multi sub infix:<+>(Math::SparseMatrix:D $m1, Math::SparseMatrix:D $m2) is export(:ALL) { $m1.add($m2)}
multi sub infix:<+>(Math::SparseMatrix:D $m, Numeric:D $scalar) is export(:ALL) { $m.add($scalar) }
multi sub infix:<+>( Numeric:D $scalar, Math::SparseMatrix:D $m) is export(:ALL) { $m.add($scalar) }

multi sub infix:<->(Math::SparseMatrix:D $m1, Math::SparseMatrix:D $m2) is export(:ALL) { $m1.add($m2.multiply(-1))}
multi sub infix:<->(Math::SparseMatrix:D $m, Numeric:D $scalar) is export(:ALL) { $m.add(-$scalar) }
multi sub infix:<->( Numeric:D $scalar, Math::SparseMatrix:D $m) is export(:ALL) { $m.add(-$scalar) }

multi sub infix:<*>(Math::SparseMatrix:D $m1, Math::SparseMatrix:D $m2) is export(:ALL) { $m1.multiply($m2) }
multi sub infix:<*>(Math::SparseMatrix:D $m, Numeric:D $scalar) is export(:ALL) { $m.multiply($scalar)}
multi sub infix:<*>( Numeric:D $scalar, Math::SparseMatrix:D $m) is export(:ALL) { $m.multiply($scalar) }

# is equiv(&infix:<**>)
# Opt-8 for ⎡•⎦
multi sub infix:<•>( Math::SparseMatrix:D $m1, Math::SparseMatrix:D $m2) is tighter(&infix:<*>) is export(:ALL) { $m1.dot($m2) }
multi sub infix:<dot>( Math::SparseMatrix:D $m1, Math::SparseMatrix:D $m2) is equiv(&infix:<•>) is export(:ALL) { $m1.dot($m2) }

#=====================================================================
# Determinants and norms
# TBD...