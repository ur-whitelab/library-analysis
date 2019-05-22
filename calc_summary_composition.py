import matplotlib.pyplot as plt
import numpy as np
import operator
import sys

ALPHABET = ['A','R','N','D','C','Q','E','G','H','I', 'L','K','M','F','P','S','T','W','Y','V', '*']
TEMPLATE_SEQ = "TLSW*EAMDMCTDTG"

def pep_to_int_list( pep):
    '''Takes a single string of amino acids and translates to a list of ints'''
    return(list(map(ALPHABET.index, pep.replace('\n', '')))) 

if(len(sys.argv) != 2):
    print('Usage: calc_summary_composition.py [target_AA_strings_file]. Must be formatted correctly.')

fname = sys.argv[1]

with open(fname, 'r') as f:
    lines = f.readlines()

for i in range(len(lines)):
    lines[i] = lines[i].replace('\n','')

print(lines)

#do summary

summary_dict = {}
peptide_strings = [line.replace('\n','') for line in lines]
peptide_ints = [pep_to_int_list(item) for item in peptide_strings]

for line in lines:
    if line not in summary_dict.keys():
        summary_dict[line] = 1
    else:
        summary_dict[line] += 1

sorted_summary = sorted(summary_dict.items(), key = operator.itemgetter(1), reverse=True)

print(sorted_summary[:100])

with open('{}_SUMMARY.txt'.format(fname.split('.')[0]), 'w+') as f:
    for item in sorted_summary:
        f.write('{}, {}\n'.format(item[0], item[1]))

#do composition

global_pos_counts = np.zeros((len(peptide_ints[0]), len(ALPHABET)))

for i in range(len(peptide_ints)):
    for idx, value in enumerate(peptide_ints[i]):
        global_pos_counts[idx][value] += 1

np.savetxt('{}_COMPOSITION.txt'.format(fname.split('.')[0]), global_pos_counts)

#make stacked bar chart
x_indices = np.arange(len(TEMPLATE_SEQ))
plots = []
width = 0.35
bottom = np.zeros(len(TEMPLATE_SEQ))

for i in range(len(ALPHABET)):
    plt.bar(x_indices, global_pos_counts[:,i], width, bottom=bottom, label=ALPHABET[i])
    bottom += global_pos_counts[:,i]

plt.xticks(x_indices, TEMPLATE_SEQ)
plt.legend(bbox_to_anchor=(1.15,1.15))

plt.savefig('{}_composition_bar_chart.svg'.format(fname.split('.')[0]))
