
from functools import reduce
import numpy as np
from Bio import SeqIO
import CGR_Protein as CP
import DL_Protein as DLP
import CGR_DNA as CD
import DL_DNA as DLD
import os
import argparse


def write_megaa(path, name, dist):
	with open(path, "w") as f:
		f.write('#mega\n!TITLE;\n')
		for namei in name:
			f.write('#' + namei + '\n')
		l = dist.shape[0]
		for i in range(l):
			for j in range(i):
				f.write(str('{:.10f}'.format(dist[i,j])) + '\t')
			f.write('\n')
	f.close()

def mtx_distance(A, z) -> float:
    array_distance = np.mat(np.zeros((z, z)))

    for i in range(z):
        if i == z:
            break
        else:
            for j in range(i + 1, z):
                differ = A[i, :] - A[j, :]
                # print(numer)
                # print(denom)
                array_distance[i, j] = np.sum(abs(differ))  

    return array_distance.T
def get_args():

	"""Get args"""
	parser = argparse.ArgumentParser(
		description="CGRWDL",
		formatter_class=argparse.ArgumentDefaultsHelpFormatter
	)
	parser.add_argument('seqtype', type=str, choices=['dna', 'protein'],help=' the sequence stype')

	parser.add_argument('--savefolder', default='./output', type=str,help='position of save distance file')

	parser.add_argument('--k', type=int, help='kmer')

	parser.add_argument('--seqformat', default='fasta', type=str,help='the data stype')

	parser.add_argument('--files',default='./input/coding_acid.fasta')

	args = parser.parse_args()
	return args

def CGRWDL(args):
	if args.seqtype=="dna":
		k = args.k
		r = reduce(lambda x, y: [i + j for i in x for j in y], [['A', 'T', 'C', 'G']] * k)
		CGR_pks = {}
		pk = {}
		DL_pks = {}
		for i in range(0, 4 ** k):
			CGR_pks[r[i]] = [0, 0]
	

		ll = [args.files]
		
		seqs = [fa.seq for fa in SeqIO.parse(args.files, "fasta")]

		l = len(seqs)

		for i in ll:
			fasta_seq = CD.fasta_reader(i)
			# cgr result
			A = np.mat(np.zeros((l, 2 * 4 ** k)))
			# DLTree result
			B = np.mat(np.zeros((l, 4 ** k)))
			z = -1
			names = []
			for name, seq in fasta_seq:
				names.append(name)
				cgr = CD.mk_cgr(seq)
				z = z + 1
				A = CD.cgr_count(seq, cgr, CGR_pks, k, z, pk, A)
				B = DLD.DlTree(seq, DL_pks, k, r, B, z)
			E = np.append(B, B, axis=1)
		cgr_dl = np.multiply(A, E)
		z = A.shape[0]  

		array_distance = mtx_distance(cgr_dl, z)

		savefile = str(k) +'_'+ args.seqtype + '_' + splitn(args.files)
		# savefile = str(k)  + '_' + splitn(args.files)
		path = os.path.join(args.savefolder, savefile)
		write_megaa(path, names, array_distance)
		del array_distance

	elif  args.seqtype=="protein":
		k = args.k
		r = reduce(lambda x,y: [i+j for i in x for j in y], [['A','I','L','M','F','P','W','V','D','E','N','C','Q','G','S','T','Y','R','H','K']] * k)
		CGR_pks = {}
		pk = {}
		DL_pks = {}

		for i in range(0, 20 ** k):
			CGR_pks[r[i]] = [0, 0]
		

		ll = [args.files]

		seqs = [fa.seq for fa in SeqIO.parse(args.files, "fasta")]
		l = len(seqs)

		for i in ll:
			fasta_seq = CP.fasta_reader(i)
			# cgr result
			A = np.mat(np.zeros((l, 2 * 20 ** k)))
			# DLTree result
			B = np.mat(np.zeros((l, 20 ** k)))
			z = -1
			names = []
			for name, seq in fasta_seq:
				names.append(name)
				cgr = CP.mk_cgr(seq)
				z = z + 1
				A = CP.cgr_count(seq, cgr, CGR_pks, k, z, pk, A)
				B = DLP.DlTree(seq, DL_pks, k, r, B, z)
			E = np.append(B, B, axis=1)
		cgr_dl = np.multiply(A, E)
		z = A.shape[0]  

		array_distance = mtx_distance(cgr_dl, z)

		savefile = str(k) + '_'+ args.seqtype + '_' + splitn(args.files)
		#savefile = str(k)  + '_' + splitn(args.files)
		path = os.path.join(args.savefolder, savefile)
		write_megaa(path, names, array_distance)
		del array_distance
	else:
		exit('data upload error, check your input !')
def splitn(s):
    return os.path.basename(s).replace('.fasta','.meg')


if __name__=="__main__":

	args = get_args()
	if not os.path.exists(args.savefolder):
		os.makedirs(args.savefolder)
	print(' -sequence type is {}'.format(args.seqtype))
	print(' -kmer is {}'.format(args.k))
	CGRWDL(args)



