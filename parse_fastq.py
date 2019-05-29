import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator as fastqIterator
from Bio.Seq import Seq 
from Bio.Alphabet import IUPAC
from scipy.stats import mode

def printHelp():
    print('Usage: parse_fastq.py [target_FASTQ_file] [optional: barcode_designation_file]')
    exit()
'''This script parses and translates a fastq file into a number of datasets, 
   determined by the length of the barcode_designation_file. If not provided,
   created a single dataset instead. Searches for peptides with the template
   sequence TLSW*EAMDMCTDTG, using this subsequence as the reference location
   for alignment. Saves the Amino Acid sequences (AA), codons of the sequence,
   labels of the sequences, and scores to separate .txt files for each dataset.
   The barcode designation file should contain nothing but one codon per line 
   indicating the "barcodes" that should be sorted and processed separately.'''



def main():
    if(len(sys.argv) != 2 and len(sys.argv) != 3):
        printHelp()

    fastq_filename = sys.argv[1]
    barcode_filename = None
    if len(sys.argv) == 3:
        barcode_filename = sys.argv[2]
    if barcode_filename is not None:
        barcode_lines = open(barcode_filename, 'r').readlines()
        barcodes = []
        for line in barcode_lines:
            barcodes.append(line.replace('\n', ''))
        print('Barcodes: {}'.format(barcodes))
    names = []
    codon_seqs = []
    aa_seqs = []
    qualities = []
    seq_lengths = []
    with open(fastq_filename,'r') as f:
        for name, sequence, quality in fastqIterator(f):
            if 'N' not in sequence:
                seq_str = str(Seq(sequence, IUPAC.unambiguous_dna).translate())
                codon_idx = None
                if 'TLSW*' in seq_str:
                    # first region not mutated
                    codon_idx = seq_str.index('TLSW*')*3
                    aa_idx = seq_str.index('TLSW*')
                    # 129 is length back to barcode codon from start of 'TLSW*'
                    barcode_idx = codon_idx - 129
                elif '*EAMDMC' in seq_str:
                    #second region not mutated
                    codon_idx = (seq_str.index('*EAMDMC') - len('TLSW') )*3
                    aa_idx = seq_str.index('*EAMDMC') - len('TLSW')
                    barcode_idx = codon_idx - 129
                elif 'CTDTG' in seq_str:
                    #third region not mutated. shouldn't ever get here but just in case
                    codon_idx = (seq_str.index('CTDTG') - len('TLSW*EAMDM') )*3
                    aa_idx = seq_str.index('CTDTG') - len('TLSW*EAMDM')
                    barcode_idx = codon_idx - 129
                if (codon_idx is not None
                    and len(str(seq_str)[aa_idx:aa_idx + 15]) == 15
                    and barcode_idx >= 0) :
                    barcode = sequence[barcode_idx:barcode_idx + 3]
                    names.append(name)
                    codon_seqs.append(sequence[codon_idx:codon_idx + 45])
                    qualities.append(quality[codon_idx:codon_idx + 45])
                    aa_seqs.append(str(seq_str)[aa_idx:aa_idx + 15])

    fastq_filename = fastq_filename.split('.')[0]
    with open(fastq_filename +'_LABELS.txt', 'w+') as f:
        for name in names:
            f.write('{}\n'.format(name))
    with open(fastq_filename +'_CODONS.txt', 'w+') as f:
        for codon in codon_seqs:
            f.write('{}\n'.format(codon))
    with open(fastq_filename +'_SCORES.txt', 'w+') as f:
        for qual in qualities:
            f.write('{}\n'.format(qual))
    with open(fastq_filename +'_AA_SEQS.txt', 'w+') as f:
        for aa_seq in aa_seqs:
            f.write('{}\n'.format(aa_seq))
    
    
        

if __name__ == '__main__':
    main()
