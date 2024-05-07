import hdbscan
from sklearn.preprocessing import StandardScaler


def HDBSCAN(df_clustering,min_sample):

    sc = StandardScaler()
    clustering_sc = sc.fit_transform(df_clustering)
    clusterer = hdbscan.HDBSCAN(gen_min_span_tree=True, min_cluster_size=min_sample)
    clusterer.fit(clustering_sc)
    cluster_num = clusterer.labels_
    score = clusterer.relative_validity_
    
    return {
        "clu num" : cluster_num,
        "score" : score
    }