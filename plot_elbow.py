from matplotlib import pyplot as plt
import pickle
import sys

if(len(sys.argv) != 2):
    print('Usage: elbow_plot.py [target_pickle_file]. Must be formatted correctly as per peptide_kmeans.py')

fname = sys.argv[1]

data_dict = pickle.load(open(fname, 'rb'))

k_vals = list(data_dict['inertias'].keys())
inertias = list(data_dict['inertias'].values())

plt.figure()
plt.plot(k_vals, inertias)
plt.xlabel('Number of Clusters')
plt.ylabel('Inertia (Mean Squared Distance)')

plt.savefig('elbow_plot.svg')
plt.savefig('elbow_plot.png')

