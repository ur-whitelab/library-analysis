# library-analysis
Collaboration project with Rudi Fasan's group to implement a peptide library analysis tool in python and do QSAR modeling.

## Installation

1. Make sure you have all the requisite python packages installed:
* biopython
* matplotlib -- make sure you have the `'tab20'` color set
* numpy


## Running the Software

1. Change to the `library-analysis` directory.
1. Create an appropriate barcode and target sequence file, as detailed below.
1. Make sure your `run.sh` file is executable: `chmod u+x run.sh`.
1. Invoke the program: `./run.sh <fastq-filename> <barcode-filename>`. Here, `<fastq-filename>` is the name of your fastq file, and `<barcode-filename>` is the name of the barcode file, e.g.:
>`./run.sh test.fastq sample_barcode_file.txt`

This will process the chosen fastq file and create a directory of the same name. In that directory, the output will be saved. Output consists of:
* Composition bar chart (as SVG and PNG images)
* Composition text file
* Summary file, listing the counts of each amino acid sequence in the aligned target regions
* Codons file, containing codon sequences of peptides from the fastq file
* Labels file, containing the labels from the fastq file
* Scores file, containing the scores from the fastq file

### Creating the Barcode File

The barcode file follows a *strict* file format. It must have exactly three columns, with no leading empty lines. The first column contains the barcode(s) to be used as a _single word_, the second contains the codon sequence to use to search for the barcode, and the final column contains the "target region" codon sequence that corresponds to this barcode. The barcode names appearing in this file are used to generate names for output files.

The sample barcode file looks like this:
>`ATG GCTCGTATGTTGTGTGGAATTG ACCCTGTCCTGGTAGGAAGCCATGGACATGTGCACCGATACCGGA`
>
>`AAG GCTCGTAAGTTGTGTGGAATTG ACCCTGTCCTGGTAGGAAGCCATGGACATGTGCACCGATACCGGA`

Each line contains `<barcode-name> <barcode-sequence> <target-sequence>`, in that order. 
Note that the columns are separated by one space only. Valid characters that may appear in the barcode file are: "A", "C", "T", "G", "-", "PARENT", " ", and "\n" (newline). Some example valid barcode names are: `ATG-TCG`, `PARENT`, and `TTC`.

### Todo
* Make this a proper python package to ease installation.