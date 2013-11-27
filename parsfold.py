"""
parsfold.py

Use experimental PARS data to generate a set of potential foldings ignoring
energy considerations. Only returns fully matching sets of pairs in Vienna
format (using dots and brackets) - discounts those cases that are invalid.

Author:  Thomas E. Gorochowski 
Updated: 26/11/2013
"""

import csv
from Bio.Seq import Seq
from Bio.SeqUtils import *

def load_pars_data (pars_file):
	"""
	Load PARS data FASTA-like format:
	  >gene_name
	  ATGAGCGCGATA...
	  ...|||..|-||...

	Key:
	  |  structured
	  .  unstructed
	  -  unknown (not handled at moment)

	Returns dictionary with gene_name the key and each element a dict with keys
	  'seq' is the DNA sequence.
	  'pars' is the PARS data.
	"""
	all_data = {}
	f = open(pars_file, 'r')
	line = f.readline() 
	while line != '':
		if len(line) > 0 and line[0] == '>':
			k = line[1::].rstrip()
			gene_data = {}
			gene_data['seq'] = f.readline().rstrip()
			gene_data['pars'] = f.readline().rstrip()
			all_data[k] = gene_data
		line = f.readline()
	return all_data

def find_regions (pars_data):
	"""
	Find consecutive stretches of structure in the PARS data. These regions
	form the central part of the matching.
	"""
	region_data = []
	in_region = False
	region_start_idx = 0
	region_end_idx = 0
	for i in range(len(pars_data)):
		if pars_data[i] == '|' and in_region:
			# Continue existing region
			region_end_idx = region_end_idx+1
		elif pars_data[i] == '|' and not in_region:
			# Start new region
			region_start_idx = i
			region_end_idx = i+1
			in_region = True
		elif pars_data[i] == '.' and in_region:
			# End existing region
			region_data.append([region_start_idx, region_end_idx])
			in_region = False
	# Check if in region at end (if so add last region)
	if in_region:
		region_data.append([region_start_idx, region_end_idx])
	return region_data

def find_seq_regions (region_data, seq_data):
	"""
	Find all the sequences that corrispond to each of the structured regions.
	Assumes 5'->3' sequences. 
	"""
	seq_regions = []
	for r in region_data:
		seq_regions.append(seq_data[r[0]:r[1]])
	return seq_regions

def match_recursion (regions, seq_regions, list_to_add_to, sub_solution, 
	                 matched_regions, region_to_match, cur_region):
	"""
	Recursive function that generates the entire possible matching tree. If
	a non-complete matching of all regions is found then the matches for that
	branch are removed.
	"""
	r1 = str(Seq(seq_regions[region_to_match]).reverse_complement())
	r1_matched = False
	# Loop through the remaining regions
	for r2 in range(cur_region,len(regions)):
		# Check if r2 is past the furthest matched region
		if len(matched_regions) > 0 and region_to_match < max(matched_regions) and r2 > max(matched_regions):
			# If so, end branch
			return
		# See if we have found a match
		if r2 not in matched_regions and seq_regions[r2] == r1:
			# Search forward
			match_recursion(regions, seq_regions, list_to_add_to, 
				            list(sub_solution), list(matched_regions), 
				            region_to_match, r2+1)
			# Add match, continue
			sub_solution.append([region_to_match, r2])
			matched_regions.append(region_to_match)
			matched_regions.append(r2)
			r1_matched = True
			break
	# Didn't find match so dead end
	if r1_matched == False:
		return
	# Move to next region if possible, else end
	next_region = -1
	for i in range(region_to_match+1, len(regions)):
		if i not in matched_regions:
			next_region = i
			break
	if next_region != -1:
		# Start on the next free region
		match_recursion(regions, seq_regions, list_to_add_to, 
		                list(sub_solution), list(matched_regions), 
		                next_region, next_region+1)
	else:
		# At end so have a solution
		list_to_add_to.append(sub_solution)
	return list_to_add_to

def generate_all_matches (seq_data, pars_data):
	"""
	Calculates all the potential matches and outputs list of possibilities
	in Vienna format.
	"""
	# Find all regions in the pars data
	regions = find_regions(pars_data)
	# Using these extract their sequences
	seq_regions = find_seq_regions(regions, seq_data)
	# Start at beginning and... match! Use a recursive method to generate full
	# matching tree.
	return match_recursion(regions, seq_regions, [], [], [], 0, 1)

def test_hard_coded ():
	data = {}
	gene_data1 = {}
	gene_data1['seq']  = 'ATGAAGGTCGCGGACCTCTCGATACGCATAGGCTAAGCGGTATCATA'
	gene_data1['pars'] = '.....||||...||||.....||||..............||||....'
	gene_data2 = {}
	gene_data2['seq']  = 'ATGAAGGGGGCGAAAATCTCGTTTTGCATAGGCTAAGCGCCCCCATA'
	gene_data2['pars'] = '.....||||...||||.....||||..............||||....'
	data['test1'] = gene_data1
	data['test2'] = gene_data2
	matches1 = generate_all_matches(data['test1']['seq'], data['test1']['pars'])
	print matches1
	matches2 = generate_all_matches(data['test2']['seq'], data['test2']['pars'])
	print matches2

def test_from_file ():
	data = load_pars_data('test.pars')
	matches = generate_all_matches(data['test']['seq'], data['test']['pars'])
	print matches
