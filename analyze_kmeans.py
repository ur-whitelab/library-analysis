from matplotlib import pyplot as plt
from sklearn.decomposition import PCA
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
import pickle
import sys

if(len(sys.argv) != 3):
    print('Usage: analyze_kmeans.py [data_pickle_file] [distances_pickle_file]. Must be formatted correctly as per peptide_kmeans.py')

ALPHABET = ['A','R','N','D','C','Q','E','G','H','I', 'L','K','M','F','P','S','T','W','Y','V', '*']

def pep_to_int_list( pep):
    '''Takes a single string of amino acids and translates to a list of ints'''
    return(list(map(ALPHABET.index, pep.replace('\n', '')))) 

data_file = sys.argv[1]
distances_file = sys.argv[2]
data_dict = pickle.load(open(data_file, 'rb'))

with open('best_elbow_value.txt', 'r') as f:
    lines = f.readlines()
best_kval = int(lines[0].replace('\n',''))

sequence_file = 'HL2m5_lib_Shh_120_AA_SEQS.txt'
with open(sequence_file, 'r') as f:
    lines = f.readlines()
peptide_strings = [line.replace('\n','') for line in lines]
peptide_ints = [pep_to_int_list(item) for item in peptide_strings]
cluster_labels = data_dict['point_labels'][best_kval]
assert(len(cluster_labels) == len(peptide_ints))
print(peptide_strings)


cluster_pos_counts = {}

print('Calculating cluster modes...')

cluster_mode_peptides = []

for i in range(len(peptide_ints)):
    label = cluster_labels[i]
    if label not in cluster_pos_counts.keys():
        cluster_pos_counts[label] = np.zeros((len(peptide_ints[i]), len(ALPHABET)))
    for idx, value in enumerate(peptide_ints[i]):
        cluster_pos_counts[label][idx][value] += 1

for i in range(best_kval):
    pep = ''.join([ALPHABET[np.argmax(aa_dist)] for aa_dist in cluster_pos_counts[i]])
    cluster_mode_peptides.append(pep)
    print('MODE of cluster {}: {}'.format(i, pep))



print('length of point labels: {}'.format(len(data_dict['point_labels'][3])))

#load the pickled data, e.g. /home/pca_projected_distances.pickle
peptide_distances = pickle.load(open('{}'.format(distances_file),'rb'))

#project down to 2d for plotting
pca = PCA(n_components = 2)
pca.fit(peptide_distances)
two_d_distances = pca.transform(peptide_distances)

cluster_mean_positions = np.zeros((best_kval, 2))
cluster_counts = np.zeros(best_kval)
for i, dist in enumerate(two_d_distances):
    #get mean of all cluster point locations
    #print(dist, cluster_labels[i])
    cluster_mean_positions[cluster_labels[i]] += dist
    cluster_counts[cluster_labels[i]] += 1

for i in range(best_kval):
    cluster_mean_positions[i] /= cluster_counts[i]

print('Cluster mean positions: {}'.format(cluster_mean_positions))


plt.figure()
cmap=plt.cm.get_cmap('Dark2')
plt.scatter(two_d_distances[:,0], two_d_distances[:,1], c=cluster_labels, cmap=cmap, marker='.', alpha=0.3)# linewidths=0.05, 
plt.scatter(cluster_mean_positions[:,0], cluster_mean_positions[:,1], c=list(range(best_kval)), cmap=cmap, marker='*', s=[200 for i in range(best_kval)])
custom_lines = [Line2D([0], [0], color=cmap(float(i)/float(best_kval-1)), lw=4) for i in range(best_kval)]
plt.legend(custom_lines, cluster_mode_peptides)#, bbox_to_anchor=(1.1,0.5))
plt.savefig('{}_clusters_kmeans_pca.svg'.format(best_kval))
plt.savefig('{}_clusters_kmeans_pca.png'.format(best_kval))

print(data_dict.keys())
