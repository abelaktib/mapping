from copy import deepcopy
import numpy as np 
import copy
import argparse
import os
import sys
SEQ = "CACATCG"

def transf(SEQ):
    n_str = SEQ + "$"
    matrix = []
    for i in range(len(n_str)):
        ephem = []
        ephem.append(n_str[len(n_str)-1])
        ephem.append(n_str[0:len(n_str)-1])
        ephem = "".join(ephem)
        matrix.append(ephem)
        n_str = ephem
    matrix = sorted(matrix)
    return matrix

def extract(mat):
    seque = []
    for i in mat:
        seque.append(i[len(mat)-1])
    seque = "".join(seque)
    return seque

def positioning(seq):
    alpha = {}
    count = []
    for x in seq:
        if x not in alpha:
            alpha[x] = 0
        count.append(alpha[x])
        alpha[x]+= 1
    t = 0 
    for k in sorted(alpha.keys()):
        tmp = alpha[k]
        alpha[k] = t
        t+= tmp
    return alpha,count


def counting(seq):
    memoire = []
    for res in seq:
        if res not in memoire.keys():
            memoire[res] = seq.count(res)
    return memoire


def reverse(seq_a):
    (alpha,count) = positioning(bw)
    seq = "$"
    x = 0
    p = 0
    while x != "$":
        x = seq_a[p]
        seq = x +seq
        p = count[p] + alpha[x]

    return(seq[1:])

def vpos(seq_a):
    (alpha,count) = positioning(seq_a)
    taille = len(seq_a)
    pos = [-1]*taille
    p=0
    while seq_a[p] != "$":
        pos[p] = taille - 2
        p = alpha[seq_a[p]] + count[p]
        taille -= 1
    return pos


def genere_fm (bw ,alpha):
    dico = dict(alpha)
    for k  in dico:
        dico[k] =0
    fmindex = []
    for x in bw:
        fmindex.append(dict(dico))
        dico[x] += 1
    fmindex.append(dico)
    return fmindex


def trouv_motif(motif,fmindex,alpha):
    b = alpha[motif[-1]]
    f = b + fmindex[-1][motif[-1]]
    c = len (motif) - 2
    while c>= 0 and f > b :
        x = motif[c]
        tmp = b
        b = alpha [x] + fmindex [b][x]
        f = b + fmindex[f][x]-fmindex[tmp][x]
        c-=1
    return(b,f)

def trouv_pos(b,f,pos):
    read_pos = []
    for x in range (b,f):
        tmp = pos[x] +1
        if tmp == len(pos):
            tmp = 0
        read_pos.append(tmp)
    return read_pos




#############################
def get_arguments():
    """Retrieves the arguments of the program.
      Returns: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description="Alignement")
    parser.add_argument('-r', dest='ref', help="Fasta file")
    parser.add_argument('-f', dest='file', help="reads")
    return parser.parse_args()



######################################################
def readfasta(ref):
    with open(ref, 'r') as file :
        l = file.readlines()
        ref = ""
        for i in l[1:]:
           ref += i.strip()
    return ref




def readreads(files):
    with open(files, 'r') as file :
        lines = file.readlines()
        ref = []
        for line in lines:
            if line.startswith(">"):
                continue
            ref.append(line.strip())
    return ref
        

def main(ref,file):
    reference = readfasta(ref)
    reads = readreads(file)
    mat = transf(reference)
    bw = extract(mat) 
    (alpha,count) = positioning(bw)
    pos = vpos(bw)
    fmindex = genere_fm(bw,alpha)
    with open ("result.txt","w") as r:
        for read in reads :
            b,f = trouv_motif(read,fmindex,alpha)
            xpos = str(trouv_pos(b,f,pos))
            long = len(read)
            xpos = xpos.replace("[","").replace("]","")
            r.write(f" {xpos}  {long}\n ")

        
    

if __name__ == '__main__':
    result_args = get_arguments()
    main(result_args.ref,result_args.file)