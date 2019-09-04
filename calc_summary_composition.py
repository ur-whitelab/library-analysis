from Bio.Alphabet import IUPAC
from Bio.Seq import Seq
from matplotlib.cm import get_cmap
import matplotlib.pyplot as plt
import numpy as np
import operator
from sys import argv

ALPHABET = ['A','R','N','D','C','Q','E','G','H','I', 'L','K','M','F','P','S','T','W','Y','V', '*']

def pep_to_int_list(pep):
    '''Takes a single string of amino acids and translates to a list of ints'''
    return(list(map(ALPHABET.index, pep.replace('\n', '')))) 

if(len(argv) != 3):
    print('Usage: calc_summary_composition.py [target_AA_seqs_file] [barcodes_file]. Must be formatted correctly.')
    exit(1)

seqs_fname = argv[1]
barcode_fname = argv[2]

barcode = seqs_fname.split('_')[-3]
TEMPLATE_SEQ = ''
with open(barcode_fname, 'r') as f:
    lines = f.read().splitlines()
for line in lines:
    this_barcode = line.split()[0]
    if this_barcode == barcode:
        TEMPLATE_SEQ = str(Seq(line.split()[2],
                                    IUPAC.unambiguous_dna).translate())
    elif barcode == 'UNSORTED':
        TEMPLATE_SEQ = '*'*(len(str(Seq(lines[0].split()[2],
                                    IUPAC.unambiguous_dna).translate())) - 1)

with open(seqs_fname, 'r') as f:
    lines = f.read().splitlines()
    if len(lines) == 0:
        print('FOR BARCODE {}, NO SEQUENCES WERE FOUND!'.format(barcode))
        exit(0)

#do summary

summary_dict = {}
peptide_ints = [pep_to_int_list(line) for line in lines]

for line in lines:
    if line not in summary_dict.keys():
        summary_dict[line] = 1
    else:
        summary_dict[line] += 1

sorted_summary = sorted(summary_dict.items(), key = operator.itemgetter(1), reverse=True)

with open('{}_SUMMARY.txt'.format(seqs_fname.split('.')[0]), 'w+') as f:
    for item in sorted_summary:
        f.write('{}, {}\n'.format(item[0], item[1]))

#do composition

global_pos_counts = np.zeros((len(peptide_ints[0]), len(ALPHABET)))

for i in range(len(peptide_ints)):
    for idx, value in enumerate(peptide_ints[i]):
        global_pos_counts[idx][value] += 1

np.savetxt('{}_RAW_COMPOSITION.txt'.format(seqs_fname.split('.')[0]), global_pos_counts)

#make stacked bar chart
x_indices = np.arange(global_pos_counts.shape[0])
plots = []
width = 0.35
bottom = np.zeros(global_pos_counts.shape[0])
#do it this way to avoid matplotlib colormap version conflicts
color_list=['r','g','b','c','m','y','k','w']

for i in range(len(ALPHABET)-1, -1, -1):
    plt.bar(x_indices, global_pos_counts[:,i], width, bottom=bottom, label=ALPHABET[i], edgecolor='black', color=color_list[i%len(color_list)])
    bottom += global_pos_counts[:,i]

plt.xticks(x_indices, TEMPLATE_SEQ)
plt.legend(bbox_to_anchor=(1.15,1.15))

plt.savefig('{}_composition_bar_chart.svg'.format(seqs_fname.split('.')[0]))
plt.savefig('{}_composition_bar_chart.png'.format(seqs_fname.split('.')[0]))
