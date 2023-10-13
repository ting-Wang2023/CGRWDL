import pandas as pd
import sys
import os
from collections import defaultdict
import matplotlib.pyplot as plt
import numpy as np
from functools import reduce
from Bio import SeqIO
import argparse

try:
	import helper
except ModuleNotFoundError:
	import helper
# defining cgr graph
# CGR_CENTER = (0.5, 0.5)
CGR_X_MAX = 1
CGR_Y_MAX = 1
CGR_X_MIN = 0
CGR_Y_MIN = 0
CGR_A = (CGR_X_MIN, CGR_Y_MIN)
CGR_T = (CGR_X_MAX, CGR_Y_MIN)
CGR_G = (CGR_X_MAX, CGR_Y_MAX)
CGR_C = (CGR_X_MIN, CGR_Y_MAX)
CGR_CENTER = ((CGR_X_MAX - CGR_Y_MIN) / 2, (CGR_Y_MAX - CGR_Y_MIN) / 2)

# Add color code for each element


def empty_dict():
	"""
	None type return vessel for defaultdict
	:return:
	"""
	return None


CGR_DICT = defaultdict(
	empty_dict,
	[
		('A', CGR_A),  # Adenine
		('T', CGR_T),  # Thymine
		('G', CGR_G),  # Guanine
		('C', CGR_C),  # Cytosine
		('U', CGR_T),  # Uracil demethylated form of thymine
		('a', CGR_A),  # Adenine
		('t', CGR_T),  # Thymine
		('g', CGR_G),  # Guanine
		('c', CGR_C),  # Cytosine
		('u', CGR_T)  # Uracil/Thymine
		]
)


def fasta_reader(fasta):
	"""Return a generator with sequence description and sequence
	:param fasta: str filename
	"""
	# TODO: modify it to be capable of reading genebank etc
	flist = SeqIO.parse(fasta, "fasta")
	for i in flist:

		yield i.description, i.seq


def mk_cgr(seq):

	"""Generate cgr
	:param seq: list of nucleotide
	:return cgr: [['nt', (x, y)]] List[List[Tuple(float, float)]]
	"""
	cgr = []
	cgr_marker = CGR_CENTER[:
		]    # The center of square which serves as first marker
	for s in seq:
		cgr_corner = CGR_DICT[s]
		if cgr_corner:
			cgr_marker = (
				(cgr_corner[0] + cgr_marker[0]) / 2,
				(cgr_corner[1] + cgr_marker[1]) / 2
			)
			cgr.append([s, cgr_marker])
            
		else:
			sys.stderr.write("Bad Nucleotide: " + s + " \n")
        
	return cgr
    

def cgr_count(seq,cgr,pks,k,z,pk,A):
    n=0
    for i in range(len(seq) - k + 1):
            kmer = seq[i:i+k]
            n=n+1
            if kmer not in pk.keys():
                    pk[kmer] = []
            pk[kmer].append(cgr[n+k-2][1])
            
    for key,value in pk.items():
        p=np.sum([value[i][0] for i in range(len(value))])/len(value)
        q=np.sum([value[i][1] for i in range(len(value))])/len(value)
        pks[key]=[p,q]        
    pks=dict(sorted(pks.items(), key=lambda kv: (kv[0])))
    l=-1
    for value in pks.values():
            l=l+1
            A[z,l]=value[0]  #x coordinate
            A[z,4**k+l]=value[1] #y coordinate
    return A







