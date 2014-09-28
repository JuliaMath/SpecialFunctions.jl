#!/bin/bash

for IFN in $(ls *.f); do
    OFN=${IFN%.*}.jl0

    # comments; convert to lowercase
    perl -p -e 's/^C(.*)/#\L\1/' $IFN > $OFN

    # join line breaks
    perl -0pi -e 's/\n     \S[ ]+/ /g' $OFN

    # subroutine, end
    perl -pi -e 's/^      SUBROUTINE/function/' $OFN
    perl -pi -e 's/^      END/end/' $OFN

    # function
    perl -pi -e 's/^      DOUBLE PRECISION FUNCTION (\w+)\(([\w, ]+)\)/function \1\(\2\)\n\@DOUBLE_PRECISION \1/' $OFN

    # handle these in the next pass
    perl -pi -e 's|^      DATA (.+)/ *$|\@DATA (\1,)|' $OFN
    perl -pi -e 's/^      DIMENSION /\@DIMENSION /' $OFN
    perl -pi -e 's/^      DOUBLE PRECISION /\@DOUBLE_PRECISION /' $OFN
    perl -pi -e 's/^      INTEGER /\@INTEGER /' $OFN

    # 'DO NNN X=Y ... NNN CONTINUE' --> 'DO X=Y ... NNN CONTINUE; end'
    # (run multiple times for nested loops)
    for i in 1 2 3; do
        perl -0pi -e 's/\n([ ]+)DO (\d+) ([A-Z][A-Z=0-9,]+)[ ]*\n(.+)(\n +\2 +CONTINUE)/\n\1DO \3\n\4\5\n\1end/sg' $OFN
    done

    # 'DO X=A,B[,C]' --> 'for X=A[:C]:B'
    perl -pi -e 's/ DO (\w+) *= *(\w+) *, *(\w+) *$/ for \1 = \2:\3/' $OFN
    perl -pi -e 's/ DO (\w+) *= *(\w+) *, *(\w+) *, *(\w+) *$/ for \1 = \2:\4:\3/' $OFN

    # goto, continue, return
    perl -pi -e 's/GO TO (\d+)$/\@goto line\1/' $OFN
    perl -pi -e 's/^ +(\d+) +CONTINUE$/      \@label line\1/' $OFN
    perl -pi -e 's/^   (\d\d) (\S.+)$/      \@label line\1\n      \2/' $OFN
    perl -pi -e 's/RETURN$/return/' $OFN

    # if statements
    perl -pi -e 's/^([ ]+)IF ?\(([^=]+)\)(.+)$/\1if \2\n\1    \3\n\1end/' $OFN

    # call statements
    perl -pi -e 's/^([ ]+)CALL /\1/' $OFN

    # doubles
    perl -pi -e 's/(\d+\.\d+)D0/\1/g' $OFN
    perl -pi -e 's/(\d+\.\d*)D([+-]\d+)/\1e\2/g' $OFN

    # powers
    perl -pi -e '/^[^#]/ && s/(\w *)\*\*( *\w)/\1\^\2/g' $OFN

    # comparison & logical operators
    perl -pi -e 's/\.GT\./ > /g' $OFN
    perl -pi -e 's/\.LT\./ < /g' $OFN
    perl -pi -e 's/\.GE\./ >= /g' $OFN
    perl -pi -e 's/\.LE\./ <= /g' $OFN
    perl -pi -e 's/\.EQ\./ == /g' $OFN
    perl -pi -e 's/\.NE\./ != /g' $OFN
    perl -pi -e 's/\.OR\./ || /g' $OFN
    perl -pi -e 's/\.AND\./ && /g' $OFN

    # all done; convert comments back to uppercase
    perl -pi -e 's/^#(.+)/#\U\1/' $OFN
done
