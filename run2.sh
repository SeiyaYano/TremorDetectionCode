#!/bin/bash

year=2010               #year of interest
savepath=result          #save dir

python3 src/duplication_removal.py --year $year --savepath $savepath
python3 src/event_features.py --year $year --savepath $savepath