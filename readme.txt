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
1. In the "raws.in" file, list the paths of all SAC files you Usage
2. Set the year of interest in "run.sh"
3. Specify the analysis period in "run.sh" by setting "current" as the first day and "end" as the day after the last day
4. Run "run.sh"

Event Feature Calculations ("run2.sh")
1. Set the year of interest in "run2.sh"
2. Run "run2.sh"
3. When prompted with "time difference of duplication is automatically determined as:**** [sec]" and "Do you use this value as threshold? [y/n]", check the figure "savepath/year/figs/duplication.png" and confirm that the vertical dashed line separates duplicate and non-duplicate events. If it seems to work well, enter "y". Otherwise, enter the threshold of the time difference [s].

Tremor Catalog Compilation ("run3.sh")
1. Set the year of interest in "run2.sh". If you want to analyze whole period, you can specify "all" instead of year.
2. When prompted with "label the cluster *:", refer to the figures in "savepath/result/clustering_result/clustering" and label the clusters appropriately, such as "tremor", "fast_shallow", "anthropogenic", etc. Note that the cluster "-1" is supposed to be labeled as "noise".
3. When "threshold for T:" appears, review the figures "savepath/result/clustering_result/clustering/mZB.png" and "(same)/mZB_T.png", then enter the temporal threshold for identifying isolated events. Repeat this process for "threshold for R:" with "(same)/mZB_R.png" and "threshold for eta:" with "(same)/mZB_eta.png".
4. When "Do you finish compiling quality A catalog? [y/n]" appears, examine the figures "savepath/result/clustering_result/clustering/QualityA.png" and "(same)/isolated.png". If the removal of isolated events seems successful, enter "y". Otherwise, you can adjust the thresholds for T, R, and eta again.
5. When "maximum depth [km]:" and "minimum depth [km]:" appear, after examininig the figures in "savepath/result/clustering_result/clustering", determine and enter the minimum and maximum depth of events to be reclassified.