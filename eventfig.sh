#!/bin/bash

savepath=result
foldername=clustering_result
outfile=png     #png or eps
text=N          #Y or N
catalog=tremor  #tremor or all or original

mkdir -p $savepath/{eventfig,eventfig/eps}

cat eventfig.txt | parallel -j 3 -a - python3 src/eventfig.py --savepath $savepath --foldername $foldername --outfile $outfile --text $text --catalog $catalog --ymdhms