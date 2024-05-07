#!/bin/bash

year=2010               #year of interest
savepath=result          #save dir

current=$year-01-01
#end=$year-12-01
end=$(($year + 1))-01-01
    #current; first day / end; the day after the last day

data=raws.in

if [ -e date.txt ]; then
    rm date.txt
fi

while [ ! $current = $end ] ; do
    echo `date -d "$current" "+%y%m%d"` >> date.txt
    current=`date -d "$current 1day" "+%Y-%m-%d"`
done

mkdir -p $savepath/$year/{csv_files,figs,json}

# if parallel is available
cat date.txt | parallel -j 3 -a - python3 src/main.py --format $data --savepath $savepath --quiet --date 

# if not
#while read line
#do
#  python3 src/main.py --format hinet.in --date $line
#done < date.txt