#!/bin/bash
# 31/10/2023: shuffle the order of the rows in the input job array file
# so you dont get similar jobs being ran successively on the same worker
# since these can step on each others toes.

# NOTE GIBLIN! YOU NEED TO ALTER THE FOLLOWING:
# - input
# - nrows
# - awk line; $2, $3, $4, $5 for X-maps; $2, $3, $4 for auto.

# THE INPUT SHOULD HAVE THE HEADER REMOVED!
input=config_cross_S.txt 
output=${input}_rand
nrows=322

# shuffle the order of the rows:
awk 'BEGIN {srand()} {print rand(), $0}' $input | sort -n | cut -d ' ' -f2- > randomfile.txt

# now make it so the first column (job number) is set back to sequential
awk '{print $2, $3, $4, $5}' < randomfile.txt > randomfile2.txt # get rid of first number 
for i in `seq 1 $nrows`; do echo $i >> tmp2.txt; done
paste tmp2.txt randomfile2.txt > $output

rm -f randomfile.txt randomfile2.txt tmp2.txt 


#awk '{print $2, $3, $4, $5}' < randomfile.txt > randomfile2.txt
#awk '{print $0,NR}' < randomfile2.txt > randomfile2.txt
#sed = randomfile2.txt | sed 'N;s/\(.*\)\n\(.*\)/\2 \1/'

