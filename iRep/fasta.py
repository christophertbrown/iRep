#!/usr/bin/env python3

"""
script for parsing a fasta file
"""

import sys
import os

def iterate_fasta(fasta, length=0, string = False):
	sequence = []
	if type(fasta) is str and string is False:
		fasta = open(fasta)
	elif type(fasta) is str and string is True:
		fasta = fasta.split('\n')
	for line in fasta:
		if line == '\n':
			continue
		sequence, formatted = parse_fasta(line, sequence, length)
		if formatted != []:
			yield formatted
	yield format_print(sequence, length)

def parse_fasta(line, sequence, length=0):
	line = line.strip()
	formatted = []
	if line.startswith('>'):
		if sequence != []:
			formatted = format_print(sequence, length)
		sequence = [line, []]
	else:
		sequence[1].append(line.replace(' ', ''))
	return sequence, formatted

def format_print(sequence, length=0):
	if sequence == []:
		return [[], []]
	if length == 0:
		formatted = [sequence[0], ''.join(sequence[1])]
	else:
		sequence[1] = ''.join(sequence[1])
		formatted =  [sequence[0], '\n'.join(sequence[1][i:i+length] for i in range(0, len(sequence[1]), length))]
	return formatted

if __name__ == "__main__":
	if len(sys.argv) == 1:
		for sequence in iterate_fasta(sys.stdin):
			print('\n'.join(sequence))
	elif len(sys.argv) == 2:
		length = int(sys.argv[1])
		for sequence in iterate_fasta(sys.stdin, length):
			print('\n'.join(sequence))
	elif len(sys.argv) != 3:
		print('please specify the fasta file and the number of characters \
'			'to print on each line, or 0 to print all characters on one line')
		exit()
	else:
		fasta, length = sys.argv[1], int(sys.argv[2])
		for sequence in iterate_fasta(fasta, length):
			print('\n'.join(sequence))
