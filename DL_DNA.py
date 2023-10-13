from Bio import SeqIO
from numpy import mat
import numpy as np
import copy
from functools import reduce
from collections import Counter

def get_all_kmer(string,k):
    length = len(string)
    return[string[i:i+k]for i in range(length - k+1)]

#seqs = [fa.seq for fa in SeqIO.parse("E:/xtu/shuju/HPV/326/test.fasta","fasta")]
#print(seqs)


def DlTree(s,pks,k,r,A,z):
    #A = np.mat(np.zeros((len(seqs), 4 ** k)))
    #r = reduce(lambda x, y: [i + j for i in x for j in y], [['A', 'T', 'C', 'G']] * k)

    for i in range(0, 4 ** k):
        pks[r[i]] = 0
    ck = Counter(get_all_kmer(s, k))
    k_1 = Counter(get_all_kmer(s, 1))  # k=1
    k_2 = Counter(get_all_kmer(s, k - 1))  # k=k-1
    # Counter(get_all_kmer(s,k-1))

    for key in ck.keys():
        p = ck[key] / (len(s) - k + 1)
        q_1 = k_1[key[0]] / len(s)
        q_2 = k_2[key[1:k]] / (len(s) - k + 2)
        q_3 = k_2[key[0:k - 1]] / (len(s) - k + 2)
        q_4 = k_1[key[k - 1]] / len(s)
        # ck=Counter(get_all_kmer(s,k)) 
        f = (p - ((q_1 * q_2 + q_3 * q_4) / 2)) / ((q_1 * q_2 + q_3 * q_4) / 2)
        pks[key] = f

    pks = dict(sorted(pks.items(), key=lambda kv: (kv[0])))
    # print(pks)

    l = -1
    for value in pks.values():
        l = l + 1
        A[z, l] = value
    return A



