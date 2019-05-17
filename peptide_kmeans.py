import numpy as np
from sklearn.cluster import KMeans
from sklearn.neighbors import DistanceMetric
from sklearn.decomposition import PCA
import sys

import pickle

def printHelp():
    print('Usage: peptide_kmeans.py [target_file] [distance_filename (NO EXTENSION)]')
    exit()
    
ALPHABET = ['A','R','N','D','C','Q','E','G','H','I', 'L','K','M','F','P','S','T','W','Y','V','*']

CLUSTERS_RANGE = range(2, 23)
        
def pep_to_int_list( pep):
    '''Takes a single string of amino acids and translates to a list of ints'''
    return(list(map(ALPHABET.index, pep.replace('\n', '')))) 

def main():
    if(len(sys.argv) != 3):
        printHelp()
    fname = sys.argv[1]
    pickle_filename = sys.argv[2]
    with open(fname, 'r') as f:
        lines = f.readlines()
    peptide_strings = []
    peptides = []
    for line in lines:
        peptide_strings.append(line.replace('\n',''))
    for pep_str in peptide_strings:
        peptides.append(pep_to_int_list(pep_str))
    peptides = np.array(peptides)
    print('Peptide data: {} peptides of length {}'.format(
        np.size(peptides,0),
        np.size(peptides,1))
    )

    dist = DistanceMetric.get_metric('hamming')
    #for large datasets, this will require a lot of memory!
    peptides_distance_matrix = dist.pairwise(peptides)
    pca = PCA(n_components = 12)
    pca.fit(peptides_distance_matrix)
    projected_vectors = pca.transform(peptides_distance_matrix)

    cluster_centers = {}
    inertias = {}
    point_labels = {}
    for num_clusters in CLUSTERS_RANGE:
        print('Fitting for {} clusters...'.format(num_clusters))
        clusters = KMeans(n_clusters=num_clusters, max_iter=500)
        clusters.fit(projected_vectors)
        cluster_centers[num_clusters] = clusters.cluster_centers_
        inertias[num_clusters] = clusters.inertia_
        point_labels[num_clusters] = clusters.labels_

    data_dict = {'cluster_centers':cluster_centers,
                 'inertias':inertias,
                 'point_labels':point_labels}
    pickle.dump(data_dict, open('trained_models.pickle', 'wb+'))
    #e.g. /home/pca_projected_distances.pickle
    pickle.dump(projected_vectors, open('{}.pickle'.format(pickle_filename),'wb+'))


if __name__ == '__main__':
    main()
