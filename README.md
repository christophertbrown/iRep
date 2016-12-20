# iRep.py and bPTR.py

Described in:

“Measurement of bacterial replication rates in microbial communities” - Christopher T. Brown, Matthew R. Olm, Brian C. Thomas, Jillian F. Banfield

http://dx.doi.org/10.1101/057992

Requires python3 and the following packages:
lmfit
numpy
scipy
pandas
seaborn
matplotlib

##Example usage:

### iRep
./iRep/iRep.py -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.iRep

### bPTR
./iRep/bPTR.py -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.bPTR.tsv -plot test.bPTR.pdf -m coverage

### GC Skew
./iRep/gc_skew.py -f sample_data/l_gasseri.fna

#### Example output is provided in sample_output. 

## Installation:

`pip install iRep`
