# iRep.py and bPTR.py

Described in:

“Measurement of bacterial replication rates in microbial communities” - Christopher T. Brown, Matthew R. Olm, Brian C. Thomas, Jillian F. Banfield (Nature Biotechnology 2016).

http://dx.doi.org/10.1038/nbt.3704

Requires python3 and the following packages:
lmfit
numpy
scipy
pandas
seaborn
matplotlib

## Example usage:

### iRep
./bin/iRep -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.iRep

### bPTR
./bin/bPTR -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.bPTR.tsv -plot test.bPTR.pdf -m coverage

### GC Skew
./bin/gc_skew -f sample_data/l_gasseri.fna

#### Example output is provided in sample_output. 

## Installation:

`pip install iRep`

## Usage notes:

Running bPTR and iRep requires that you have genome sequences for each of the organisms for which you want to measure replication rates. Each program takes the file paths to each of the genomes, in separate FASTA files, as input. bPTR requires complete (closed) genome sequences, and iRep requires high-quality draft genome (≥75% complete, ≤175 fragments/mbp sequencing, and ≤2% contamination). Both methods are most accurate when the genomes have been assembled from metagenomes from the samples being studied, or if it is known that an organism is present in the system with a highly similar genome sequence.

The second set of inputs are mapping files in SAM format, which you can generate by mapping your DNA sequencing reads against all of your assembled genomes using Bowtie 2. In this case, provide the path to each SAM file generated for each sample. IMPORTANT: this requires ordered SAM files that can be generated using the Bowtie2 --reorder flag.

bPTR.py and iRep.py output both a table (TSV) with the results and a PDF with plots showing genome coverage and how the results were determined.
