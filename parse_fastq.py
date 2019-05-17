import sys
from Bio.SeqIO.QualityIO import FastqGeneralIterator as fastqIterator
from Bio.Seq import Seq 
from Bio.Alphabet import IUPAC
from scipy.stats import mode

def printHelp():
    print('Usage: parse_fastq.py [target_FASTQ_file]')
    exit()


def main():
    if(len(sys.argv) != 2):
        printHelp()

    fname = sys.argv[1]
    names = []
    codon_seqs = []
    aa_seqs = []
    qualities = []
    seq_lengths = []
    with open(fname,'r') as f:
        for name, sequence, quality in fastqIterator(f):
            if 'N' not in sequence:
                seq_str = str(Seq(sequence, IUPAC.unambiguous_dna).translate())
                codon_idx = None
                if 'TLSW*' in seq_str:
                    #first region not mutated
                    codon_idx = seq_str.index('TLSW*')*3
                    aa_idx = seq_str.index('TLSW*')
                elif '*EAMDMC' in seq_str:
                    #second region not mutated
                    codon_idx = (seq_str.index('*EAMDMC') - len('TLSW') )*3
                    aa_idx = seq_str.index('*EAMDMC') - len('TLSW')
                elif 'CTDTG' in seq_str:
                    #third region not mutated. shouldn't ever get here but just in case
                    codon_idx = (seq_str.index('CTDTG') - len('TLSW*EAMDM') )*3
                    aa_idx = seq_str.index('CTDTG') - len('TLSW*EAMDM')
                if codon_idx is not None and len(str(seq_str)[aa_idx:aa_idx + 15]) == 15:
                    names.append(name)
                    codon_seqs.append(sequence[codon_idx:codon_idx + 45])
                    qualities.append(quality[codon_idx:codon_idx + 45])
                    aa_seqs.append(str(seq_str)[aa_idx:aa_idx + 15])
                    #print('processing {}...'.format(str(seq_str)[aa_idx:aa_idx + 15]))


    fname = fname.split('.')[0]
    with open(fname +'_LABELS.txt', 'w+') as f:
        for name in names:
            f.write('{}\n'.format(name))
    with open(fname +'_CODONS.txt', 'w+') as f:
        for codon in codon_seqs:
            f.write('{}\n'.format(codon))
    with open(fname +'_SCORES.txt', 'w+') as f:
        for qual in qualities:
            f.write('{}\n'.format(qual))
    with open(fname +'_AA_SEQS.txt', 'w+') as f:
        for aa_seq in aa_seqs:
            f.write('{}\n'.format(aa_seq))
    
    
        

if __name__ == '__main__':
    main()
