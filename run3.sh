#!/bin/bash

year=2010            #year or "all"
savepath=result

foldername=clustering_result_$year

mkdir -p $savepath/$foldername/{figs,clustering}

python3 src/clustering.py --year $year --savepath $savepath --foldername $foldername