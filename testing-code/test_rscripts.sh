#!/bin/bash

# tests of the Rscript XXXX scripts.

###
# Simple-minded ones:

TEMPFILE=$(mktemp tempXXXXX.RData)

# TEST 1
echo "Test 1:"
Rscript ../scripts/make-genmat.R -c <(echo '{ "bases":["X","O"], "mutpats":[[["X","O"],["O","X"]]], "mutrates":[1] }') -w 2 -o $TEMPFILE 2>/dev/null
TEST1=$( diff <(Rscript -e "library(contextual); load(\"$TEMPFILE\"); as.matrix(genmatrix)" 2>/dev/null | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') <(( cat <<END
      XX OX XO OO
   XX  0  1  1  0
   OX  1  0  0  1
   XO  1  0  0  1
   OO  0  1  1  0
END
) | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') )

if [ $TEST1 ]; then
    echo "  Failed test 1."
else
    echo "  Passed."
fi

# Test 2
echo "Test 2:"
Rscript ../scripts/make-genmat.R -c <(echo '{ "bases":["X","O"], "mutpats":[[["X","O"]],[["O","X"]]], "mutrates":[1,2] }') -w 2 -o $TEMPFILE 2>/dev/null
TEST2=$( diff <(Rscript -e "library(contextual); load(\"$TEMPFILE\"); as.matrix(genmatrix)" 2>/dev/null | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') <(( cat <<END
      XX OX XO OO
   XX  0  1  1  0
   OX  2  0  0  1
   XO  2  0  0  1
   OO  0  2  2  0
END
) | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') )

if [ "$TEST2" ]; then
    echo "Failed test 2."
else
    echo "Passed."
fi

# Test 3
Rscript ../scripts/make-genmat.R -c <(echo '{ "bases":["X","O"], "mutpats":[[["X","O"]],[["O","X"]],[["OX","XX"]],[["XO","XX"]]], "mutrates":[1,2,3,5] }') -w 3 -o $TEMPFILE 2>/dev/null
TEST2=$( diff <(Rscript -e "library(contextual); load(\"$TEMPFILE\"); as.matrix(genmatrix)" 2>/dev/null | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') <(( cat <<END
    XXX OXX XOX OOX XXO OXO XOO OOO
XXX   0   1   1   0   1   0   0   0
OXX   5   0   0   1   0   1   0   0
XOX  10   0   0   1   0   0   1   0
OOX   0   5   2   0   0   0   0   1
XXO   7   0   0   0   0   1   1   0
OXO   0   7   0   0   5   0   0   1
XOO   0   0   2   0   7   0   0   1
OOO   0   0   0   2   0   2   2   0
END
) | sed -e 's/^ *//' -e 's/ *$//' -e 's/  */ /g') )

if [ $TEST3 ]; then
    echo "Failed test 3."
else
    echo "Passed."
fi

# Test 4
Rscript -e "library(contextual); x <- read.counts('unit-test.counts',leftwin=0); stopifnot(all(as.numeric(x@counts)==c(5,5,5,2,6,6,3,4,3,6,6,4,6,6,5,6))); stopifnot(all(x@colpatterns\$sp2==c('A','T','A','T'))); stopifnot(all(x@colpatterns\$sp3==c('A','A','T','T'))); cat('read in count data: passed\n') "

rm $TEMPFILE
