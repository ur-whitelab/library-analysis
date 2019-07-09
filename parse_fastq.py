import re
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

    names = {}
    codon_seqs = {}
    aa_seqs = {}
    qualities = {}

    for barcode in barcodes:
        names[barcode] = []
        codon_seqs[barcode] = []
        aa_seqs[barcode] = []
        qualities[barcode] = []

    hits = 0
    misses = 0
    count = 0
    barcodes_missing = 0
    bad_barcodes = 0
    sequence_end_missing = 0
    template_not_found = 0
    bad_reads = 0
    with open(fastq_filename,'r') as f:
        for name, sequence, quality in fastqIterator(f):
            count += 1
            if 'N' not in sequence:
                seq_str = str(Seq(sequence, IUPAC.unambiguous_dna).translate())
                codon_idx = None
                temp_str = 'TLSW\*'
                for i, item in enumerate(temp_str):
                    re_str = temp_str[:i]+'.'+temp_str[i+1:]
                    try:
                        start_idx = re.search(re_str, seq_str).start()
                    except AttributeError:
                        start_idx = None
                    if start_idx is not None:
                        codon_idx = start_idx * 3
                        aa_idx = start_idx
                        barcode_idx = codon_idx - 132
                temp_str = '\*EAMDMC'
                if codon_idx is None:
                    for i, item in enumerate(temp_str):
                        re_str = temp_str[:i]+'.'+temp_str[i+1:]
                        try:
                            start_idx = re.search(re_str, seq_str).start()
                        except AttributeError:
                            start_idx = None
                        if start_idx is not None:
                            codon_idx = (start_idx - len('TLSW')) * 3
                            aa_idx = start_idx - len('TLSW')
                            barcode_idx = codon_idx - 132
                if codon_idx is None:
                    temp_str = 'CTDTG'
                    for i, item in enumerate(temp_str):
                        re_str = temp_str[:i]+'.'+temp_str[i+1:]
                        try:
                            start_idx = re.search(re_str, seq_str).start()
                        except AttributeError:
                            start_idx = None
                        if start_idx is not None:
                            codon_idx = (start_idx - len('TLSW*EAMDMC')) * 3
                            aa_idx = start_idx - len('TLSW*EAMDMC')
                            barcode_idx = codon_idx - 132
                if (codon_idx is not None
                    and len(str(seq_str)[aa_idx:aa_idx + 15]) == 15
                    and barcode_idx >= 0) :
                    barcode = sequence[barcode_idx:barcode_idx + 3]
                    if barcode in barcodes:
                        names[barcode].append(name)
                        codon_seqs[barcode].append(sequence[codon_idx:codon_idx + 45])
                        qualities[barcode].append(quality[codon_idx:codon_idx + 45])
                        aa_seqs[barcode].append(str(seq_str)[aa_idx:aa_idx + 15])
                        #print('barcode for {} is {}'.
                        #      format(aa_seqs[barcode][-1], barcode, sequence[6:9]))
                        hits += 1
                    else:
                        misses += 1
                        bad_barcodes += 1
                elif barcode_idx < 0:
                    misses += 1
                    barcodes_missing += 1
                elif len(str(seq_str)[aa_idx:aa_idx + 15]) != 15:
                    misses += 1
                    sequence_end_missing += 1
                elif codon_idx is None:
                    misses += 1
                    template_not_found += 1
            else:
                misses += 1
                bad_reads += 1
    print(
        'Hits: {}\nMisses: {}\n\
        Missing Barcodes: {}\n\
        Incorrect Barcodes: {}\n\
        Incomplete Template Sequence: {}\n\
        Template Not Found: {}\n\
        Low-Quality Reads: {}\nTotal: {}'.format(hits,
                                                 misses,
                                                 barcodes_missing,
                                                 bad_barcodes,
                                                 sequence_end_missing,
                                                 template_not_found,
                                                 bad_reads,
                                                 count)
    )

    for barcode in barcodes:
        filename = fastq_filename.split('.')[0] + '_{}'.format(barcode)
        with open(filename +'_LABELS.txt', 'w+') as f:
            for name in names[barcode]:
                f.write('{}\n'.format(name))
        with open(filename +'_CODONS.txt', 'w+') as f:
            for codon in codon_seqs[barcode]:
                f.write('{}\n'.format(codon))
        with open(filename +'_SCORES.txt', 'w+') as f:
            for qual in qualities[barcode]:
                f.write('{}\n'.format(qual))
        with open(filename +'_AA_SEQS.txt', 'w+') as f:
            for aa_seq in aa_seqs[barcode]:
                f.write('{}\n'.format(aa_seq))

if __name__ == '__main__':
    main()
