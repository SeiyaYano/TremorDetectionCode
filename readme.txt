Corresponding Paper
    Yano, S., & Ide, S. (2024)
    Event-feature-based clustering reveals continuous distribution of tectonic tremors of 0.3–100 s: Application to western Japan
    Geophysical Research Letters, 51, e2024GL108874
    https://doi.org/10.1029/2024GL108874

Structure
TremorDetectionCode
    JMA2001
    result (this directory will be made through the process)
        year1
            csv_files
            figs
            json
        year2
        year3
            ⋮
        clustering_result
            clustering
            figs
            events.csv
            nearest.csv
            parameters.txt
            tremor_catalog.csv
        eventfig
            figure1
            figure2
                ⋮
    src
    date.txt    (this file will be made through the process)
    errors1.txt (this file will be made through the process)
    errors2.txt (this file will be made through the process)
    eventfig.sh
    eventfig.txt
    raws.in
    readme
    run.sh
    run2.sh
    run3.sh

Usage
Seismic Event Detection ("raws.in" and "run.sh")
1. List the paths of all SAC files you use in the "raws.in" file. SAC files should be prepared for each channel with a length of one day.
2. Set the year of interest in "run.sh"
3. Specify the analysis period in "run.sh" by setting "current" as the first day and "end" as the day after the last day
4. Run "run.sh"

Event Feature Calculations ("run2.sh")
1. Set the year of interest in "run2.sh"
2. Run "run2.sh"
3. When prompted with "time difference of duplication is automatically determined as:**** [sec]" and "Do you use this value as threshold? [y/n]", check the figure "savepath/year/figs/duplication.png" and confirm that the vertical dashed line separates duplicate and non-duplicate events. If it seems to work well, enter "y". Otherwise, enter the threshold of the time difference [sec].

Tremor Catalog Compilation ("run3.sh")
*Please note that this step may not work properly if there are not enough events because we employ density-based clustering algorithm in this step.
1. Set the year of interest in "run3.sh". If you want to analyze whole period, you can specify "all" instead.
2. Set the latitudes and longitudes in "src/const" by changing "lat_range" and "lon_range".
3. When "label the cluster *:" appears, refer to the figures in "savepath/result/clustering_result/clustering" and label the clusters appropriately, such as "tremor", "fast_shallow", "anthropogenic", etc. Note that the cluster "-1" should be labeled as "noise".
4. When "threshold for T:" appears, review the figures "savepath/result/clustering_result/clustering/mZB.png" and "(same)/mZB_T.png", then enter the temporal threshold for identifying isolated events. Repeat this process for "threshold for R:" with "(same)/mZB_R.png", and repeat this process for "threshold for eta:" with "(same)/mZB_eta.png".
5. When "Do you finish compiling quality A catalog? [y/n]" appears, examine the figures "savepath/result/clustering_result/clustering/QualityA.png" and "(same)/isolated.png". If the removal of isolated events seems successful, enter "y". Otherwise, you can adjust the thresholds for T, R, and eta again.
6. When "maximum depth [km]:" and "minimum depth [km]:" appear, check the figures in "savepath/result/clustering_result/clustering" and enter the minimum and maximum depth of events to be reclassified.

You can change the set of event features used in the clustering and reclassificaion steps.
Default is ["energy dur","dep","f ratio","f high sn","f low sn"] for the clustering and ["energy dur","f low sn","f high sn","f ratio"] for the reclassification. (See Line 127-131 and Line 345 in "scr/clustering.py")
Since these sets are determined for western Japan, there might be more preferable sets for other tectonic settings.
