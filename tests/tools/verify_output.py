#! /usr/bin/env python2.7

"""A program to verify output of prophasm.

Author: Karel Brinda <kbrinda@hsph.harvard.edu>

Licence: MIT
"""

import collections
import re
import sys

in1_fn = sys.argv[1]
in2_fn = sys.argv[2]
out1_fn = sys.argv[3]
out2_fn = sys.argv[4]
inter_fn = sys.argv[5]
k = int(sys.argv[6])

reg_splitting = re.compile("[^ACGT]")

comp_dict = {
    "A": "T",
    "C": "G",
    "G": "C",
    "T": "A",
}

def load_fasta(fn):
    d=collections.defaultdict(list)
    with open(fn) as f:
        seqname="_"
        for line in f:
            line=line.strip()
            if len(line)==0:
                continue
            if line[0]==">":
                seqname, _, _ = line.partition(" ")
            else:
                d[seqname].append(line)
    dd={}
    for k in d:
        dd[k]="".join(d[k])

    return dd


def reverse_complement_str(dna):
    reverse_complement = "".join([comp_dict[x] for x in dna[::-1]])
    return reverse_complement


def get_canonical_kmers_from_fasta(fasta_fn, k):
    kmers = set()

    reg_splitting = re.compile("[^ACGT]")
    set_of_kmers = set()
    fasta_sequences = load_fasta(fasta_fn)
    for name in fasta_sequences:
        sequence = fasta_sequences[name].upper()
        sequences_ok = reg_splitting.split(sequence)
        for seq in sequences_ok:
            for i in range(len(seq) - k + 1):
                kmer = seq[i:i + k]
                kmer_rc = reverse_complement_str(kmer)
                set_of_kmers.add(min(kmer, kmer_rc))
    return set_of_kmers


print("Loading {}".format(in1_fn))
in1 = get_canonical_kmers_from_fasta(in1_fn, k)
print("Loading {}".format(in2_fn))
in2 = get_canonical_kmers_from_fasta(in2_fn, k)
print("Loading {}".format(out1_fn))
out1 = get_canonical_kmers_from_fasta(out1_fn, k)
print("Loading {}".format(out2_fn))
out2 = get_canonical_kmers_from_fasta(out2_fn, k)
print("Loading {}".format(inter_fn))
inter = get_canonical_kmers_from_fasta(inter_fn, k)

print("Is ok", in1 & in2 == inter)

s1 = in1 | in2
s2 = inter | out1 | out2
print("Is ok {} (sizes: {}, {})".format(s1 == s2, len(s1), len(s2)))
print("sym. difference: ", s1 ^ s2)
print()
print("out1 - in1")
print(out1 - in1)
print()
print("out2 - in2")
print(out2 - in2)
