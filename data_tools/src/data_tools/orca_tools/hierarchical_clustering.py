def hierarchical_clustering(points, cutoff_distance, method='single', metric='euclidean'):
    from scipy.cluster.hierarchy import linkage, fcluster
    """
    Cluster points using agglomerative hierarchical clustering.

    Parameters:
    points : numpy.ndarray
        An array of shape (n, d) representing n points in d-dimensional space.
    cutoff_distance : float
        The maximum distance between two points for them to be considered in the same cluster.
    method : str, optional
        The linkage method to use (default is 'single').
    metric : str, optional
        The distance metric to use (default is 'euclidean').

    Returns:
    numpy.ndarray
        An array of shape (n,) where the element at index i represents the cluster label of points[i].
    """
    linkage_matrix = linkage(points, method=method, metric=metric)
    cluster_labels = fcluster(linkage_matrix, t=cutoff_distance, criterion='distance')
    return cluster_labels

