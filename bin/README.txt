# iRep.py and bPTR.py

These scripts are provided for peer-review of the manuscript:

“In situ replication rates for uncultivated bacteria in natural microbial communities” - Christopher T. Brown, Matthew R. Olm, Brian C. Thomas, Jillian F. Banfield

All materials will be made publicly available upon acceptance.

Scripts require python2.7 and the following packages:
lmfit
numpy
scipy
pandas
seaborn
matplotlib
cPickle

Example usage:

# iRep
./bin/iRep.py -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.iRep

# bPTR
./bin/bPTR.py -f sample_data/l_gasseri.fna -s sample_data/l_gasseri*sam -o test.bPTR.tsv -plot test.bPTR.pdf -m coverage

Example output is provided in sample_output. 