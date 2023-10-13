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
CGR_A = (CGR_X_MIN, CGR_Y_MIN)#A,I,L,M,F,P,W,V is 0
CGR_I = (CGR_X_MIN, CGR_Y_MIN)
CGR_L = (CGR_X_MIN, CGR_Y_MIN)
CGR_M = (CGR_X_MIN, CGR_Y_MIN)
CGR_F = (CGR_X_MIN, CGR_Y_MIN)
CGR_P = (CGR_X_MIN, CGR_Y_MIN)
CGR_W = (CGR_X_MIN, CGR_Y_MIN)
CGR_V = (CGR_X_MIN, CGR_Y_MIN)
CGR_D = (CGR_X_MIN, CGR_Y_MAX)#D,E is 1
CGR_E = (CGR_X_MIN, CGR_Y_MAX)
CGR_N = (CGR_X_MAX, CGR_Y_MAX)#N,C,Q,G,S,T,Y is 2
CGR_C = (CGR_X_MAX, CGR_Y_MAX)
CGR_Q = (CGR_X_MAX, CGR_Y_MAX)
CGR_G = (CGR_X_MAX, CGR_Y_MAX)
CGR_S = (CGR_X_MAX, CGR_Y_MAX)
CGR_T = (CGR_X_MAX, CGR_Y_MAX)
CGR_Y = (CGR_X_MAX, CGR_Y_MAX)
CGR_R = (CGR_X_MAX, CGR_Y_MIN)#R,H,K is 3
CGR_H = (CGR_X_MAX, CGR_Y_MIN)
CGR_K = (CGR_X_MAX, CGR_Y_MIN)


CGR_CENTER = ((CGR_X_MAX - CGR_Y_MIN) / 2, (CGR_Y_MAX - CGR_Y_MIN) / 2)

def empty_dict():
	"""
	None type return vessel for defaultdict
	:return:
	"""
	return None


CGR_DICT = defaultdict(
	empty_dict,
	[
		('A', CGR_A),  # Alanine
		('I', CGR_I),  # Isoleucine
		('L', CGR_L),  # Leucine
		('M', CGR_M),  # Methionine
        ('F', CGR_F),  # Phenylalanine
		('P', CGR_P),  # Proline
		('W', CGR_W),  # Tryptophan
		('V', CGR_V),  # Valine
        ('D', CGR_D),  # Aspartic acid
		('E', CGR_E),  # Glutamic acid
		('N', CGR_N),  # Asparagine
		('C', CGR_C),  # Cysteine
        ('Q', CGR_Q),  # Glutanine
		('G', CGR_G),  # Glicine
		('S', CGR_S),  #Serine
		('T', CGR_T),  # Threonine
        ('Y', CGR_Y),  # Tyrosine
		('R', CGR_R),  # Arginine
		('H', CGR_H),  # Histidine
		('K', CGR_K),  # Lysine
        ('a', CGR_A),  # Alanine
		('i', CGR_I),  # Isoleucine
		('l', CGR_L),  # Leucine
		('m', CGR_M),  # Methionine
        ('f', CGR_F),  # Phenylalanine
		('p', CGR_P),  # Proline
		('w', CGR_W),  # Tryptophan
		('v', CGR_V),  # Valine
        ('d', CGR_D),  # Aspartic acid
		('e', CGR_E),  # Glutamic acid
		('n', CGR_N),  # Asparagine
		('c', CGR_C),  # Cysteine
        ('q', CGR_Q),  # Glutanine
		('g', CGR_G),  # Glicine
		('s', CGR_S),  #Serine
		('t', CGR_T),  # Threonine
        ('y', CGR_Y),  # Tyrosine
		('r', CGR_R),  # Arginine
		('h', CGR_H),  # Histidine
		('k', CGR_K),  # Lysine
		
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
def cgr_jisu(seq,cgr,pks,k,z,pk,A): 
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
            A[z,20**k+l]=value[1] #y coordinate
    return A

