import sys
import warnings
from Bio.SeqIO.QualityIO import FastqGeneralIterator as fastqIterator
from Bio import Align
from Bio import BiopythonWarning
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC

def printHelp():
    print('Usage: parse_fastq.py [target_FASTQ_file] [barcode_designation_file]')
    exit()
'''This script parses and translates a fastq file into a number of datasets, 
   determined by the length of the barcode_designation_file. If not provided,
   created a single dataset instead. Searches for peptides with the template
   sequence TLSW*EAMDMCTDTG, using this subsequence as the reference location
   for alignment. Saves the Amino Acid sequences (AA), codons of the sequence,
   labels of the sequences, and scores to separate .txt files for each dataset.
   The barcode designation file should contain nothing but one codon per line 
   indicating the "barcodes" that should be sorted and processed separately.'''

class Error(Exception):
    '''Base class for exceptions'''
    pass

class BarcodeFormatError(Error):
    '''Raised when the barcode file is incorrectly formatted'''
    pass

class BarcodeBadCharacterError(Error):
    '''Raised when the barcode file has disallowed characters.'''
    pass

class BarcodeSizeError(Error):
    '''Raised when the barcode file is too small.'''
    pass

def find_barcode_and_template(sequence, aligner, barcodes, barcode_templates,
                              sequence_templates, cutoff=0.8, retry=True):
    '''Given a BioPython aligner instance, a set of barcodes to look for, and a set of
       target sequences to check, and an alignment "score" cutoff, try to find a good
       barcode and target sequence alignment. If no barcode alignment scores above cutoff,
       will try to align with the reversed sequence instead. If that doesn't work either,
       returns None. Else, returns the barcode and aligned targe sequence.'''
    scores = [0 for barcode in barcode_templates]
    barcode, template = None, None
    for i, barcode_seq in enumerate(barcode_templates):
        barcode_alignments = aligner.align(sequence, barcode_seq)
        try:
            scores[i] = (sorted(barcode_alignments)[-1].score)/float(len(sequence_templates[i]))
        except IndexError:
            pass
    # make sure we have at least one score above some cutoff
    if any(scores) and max(scores) > cutoff: 
        maximal_idx = scores.index(max(scores))
        barcode = barcodes[maximal_idx]
        template = sequence_templates[maximal_idx]
    elif retry:
        # recurse once with reversed sequence
        barcode, template = find_barcode_and_template(sequence[::-1],
                                                      aligner,
                                                      barcodes,
                                                      barcode_templates,
                                                      sequence_templates,
                                                      retry=False)
    return(barcode, template)

ALLOWED_CHARS = set('PARENT ACTG -\n')

def main():
    if(len(sys.argv) != 3):
        printHelp()

    warnings.simplefilter('ignore', BiopythonWarning)

    fastq_filename = sys.argv[1]
    barcode_filename = sys.argv[2]
    if barcode_filename is not None:
        barcode_lines = open(barcode_filename, 'r').read().splitlines()
        try:
            if not len(barcode_lines) >= 2:
                raise BarcodeSizeError
            if ' ' not in barcode_lines[0] or \
            any([len(line.split()) != 3 for line in barcode_lines]):
                raise BarcodeFormatError
            #make a set of all chars in the barcode file
            lineset = {l for line in barcode_lines for l in line}
            if not lineset <= ALLOWED_CHARS:
                raise BarcodeBadCharacterError
        except BarcodeFormatError:
            print('ERROR: barcode file incorrectly formatted.'
                  '\nSee sample_barcode_file.txt for example of correct formatting.')
            exit(1)
        except BarcodeBadCharacterError:
            print('ERROR: barcode file has illegal '
                  'characters: {}\n'
                  ' Barcode file should'
                  ' only contain "A", "C", "T", "G", "-", "PARENT", spaces, and newlines.'
                  '\nSee sample_barcode_file.txt for '
                  'example of correct formatting.'.format(lineset - ALLOWED_CHARS))
            exit(1)
        except BarcodeSizeError:
            print('ERROR: barcode file must contain a target region and at least'
                  ' one barcode.'
                  '\nSee sample_barcode_file.txt for example of correct formatting.')
            exit(1)
        BARCODES = []
        BARCODE_TEMPLATES = []
        SEQUENCE_TEMPLATES = []
        for line in barcode_lines:
            BARCODES.append(line.split()[0])
            BARCODE_TEMPLATES.append(line.split()[1])
            SEQUENCE_TEMPLATES.append(line.split()[2])
        print('Target templates: '
              '{} ({})'.format(SEQUENCE_TEMPLATES,
                               [str(Seq(template,
                                    IUPAC.unambiguous_dna).translate())
                                for template in SEQUENCE_TEMPLATES]))
        print('Barcodes: {}'.format(BARCODES))
        print('Barcode templates: {}'.format(BARCODE_TEMPLATES))

    names = {}
    codon_seqs = {}
    aa_seqs = {}
    qualities = {}

    for barcode in BARCODES:
        names[barcode] = []
        codon_seqs[barcode] = []
        aa_seqs[barcode] = []
        qualities[barcode] = []
    names['UNSORTED'] = []
    codon_seqs['UNSORTED'] = []
    aa_seqs['UNSORTED'] = []
    qualities['UNSORTED'] = []
    
    hits = 0
    misses = 0
    count = 0
    barcodes_missing = 0
    sequence_end_missing = 0
    bad_reads = 0

    aligner = Align.PairwiseAligner()
    aligner.mode = 'local'
    #reward matches, punish mismatches
    aligner.match = 2
    aligner.mismatch = -1
    #heavily penalize gaps
    aligner.open_gap_score = -100.

    with open(fastq_filename,'r') as f:
        # this parses a fastq file into names, sequences of DNA codons, 
        # and the quality score for the sequence transcription
        print('Processing sequences...')
        for name, sequence, quality in fastqIterator(f):
            count += 1
            if 'N' not in sequence:
                this_barcode, this_template = find_barcode_and_template(sequence=sequence,
                                                                        aligner=aligner,
                                                                        barcodes=BARCODES,
                                                                        barcode_templates=BARCODE_TEMPLATES,
                                                                        sequence_templates=SEQUENCE_TEMPLATES)
                if this_barcode is None or this_template is None:
                    barcodes_missing += 1
                    misses += 1
                    #continue
                    this_barcode='UNSORTED'
                    this_template = SEQUENCE_TEMPLATES[-1]#TODO: make this assigned manually?
                seq_str = str(Seq(sequence, IUPAC.unambiguous_dna).translate())
                fwd_alignments = aligner.align(sequence, this_template)
                rev_alignments = aligner.align(sequence[::-1], this_template)
                split_fwd_alignment = str(sorted(fwd_alignments)[-1]).split('\n')
                split_rev_alignment = str(sorted(rev_alignments)[-1]).split('\n')
                fwd_score = sorted(fwd_alignments)[-1].score/float(len(this_template))
                rev_score = sorted(rev_alignments)[-1].score/float(len(this_template))
                split_alignment = split_fwd_alignment if fwd_score > rev_score else split_rev_alignment
                codon_idx = split_alignment[2].index(this_template[0])
                aligned_seq_str = split_alignment[0][codon_idx:codon_idx+len(this_template)]
                #screen for completely-missing first part of the sequence
                if '.' not in aligned_seq_str and '-' not in aligned_seq_str:
                    #print(Seq(aligned_seq_str, IUPAC.unambiguous_dna).translate())
                    translated_aligned_seq = (Seq(aligned_seq_str, IUPAC.unambiguous_dna).translate())
                    names[this_barcode].append(name)
                    codon_seqs[this_barcode].append(sequence[codon_idx:codon_idx +\
                                                             len(this_template)])
                    qualities[this_barcode].append(quality[codon_idx:codon_idx +\
                                                           len(this_template)])
                    aa_seqs[this_barcode].append(translated_aligned_seq)
                    if this_barcode != 'UNSORTED':
                        hits += 1
                elif this_barcode != 'UNSORTED':
                    misses += 1 
                    sequence_end_missing += 1
            else:
                misses += 1
                bad_reads += 1
    print(
        'Hits: {}\nMisses: {}\n\
        Missing Barcodes: {}\n\
        Incomplete Template Sequence: {}\n\
        Low-Quality Reads: {}\nTotal: {}'.format(hits,
                                                 misses,
                                                 barcodes_missing,
                                                 sequence_end_missing,
                                                 bad_reads,
                                                 count)
    )

    BARCODES.append('UNSORTED')
    for barcode in BARCODES:
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
