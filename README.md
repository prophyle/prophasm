# ProphAsm

[![Build Status](https://travis-ci.org/prophyle/prophasm.svg?branch=master)](https://travis-ci.org/prophyle/prophasm)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Prerequisities](#prerequisities)
* [Getting started](#getting-started)
* [Examples](#examples)
* [Command line parameters](#command-line-parameters)
* [Algorithm](#algorithm)
* [Similar programs](#similar-programs)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Author](#author)

<!-- vim-markdown-toc -->

## Introduction

ProphAsm is a tool for computing simplitigs from k-mer sets and for k-mer set manipulation. Simplitigs are genomic sequences computed as disjoint paths in a bidirectional vertex-centric de Bruijn graph. Compared to unitigs, simplitigs provide an improvement in the total number of sequences and their cumulative length, while both representations contain exactly the same k-mers

Upon execution, ProphAsm first loads all specified datasets (see the `-i` param) and computes their k-mer sets (see the `-k` param). If the `-x` param is provided, ProphAsm then computes their intersection, subtracts the intersection from the individual k-mer sets and computes unitigs for the intersection. If output files are specified (see the `-o` param).


## Cite

If you want to cite ProphAsm, please use the following reference:

> Brinda K, Baym M, and Kucherov G. **Simplitigs as an efficient and scalable representation of de Bruijn graphs**
. bioRxiv 2020.01.12.903443, 2020. https://doi.org/10.1101/2020.01.12.903443

If you want to cite the concept of simplitigs, please include also the parallel manuscript (the same concept discovered independently and simultaneously):

> Rahman A and Medvedev P. **Representation of k-mer sets using spectrum-preserving string sets**
. bioRxiv 2020.01.07.896928, 2020. https://doi.org/10.1101/2020.01.07.896928


## Prerequisities

* GCC 4.8+ or equivalent
* ZLib


## Getting started

```
git clone https://github.com/prophyle/prophasm
cd prophasm && make -j
./prophasm -k 15 -i tests/test1.fa -i tests/test2.fa -o _out1.fa -o _out2.fa -x _intersect.fa -s _stats.tsv
```

## Examples

```bash
prophasm -k 31 -i input.fa -o simplitigs.fa  # compute simplitigs for a single dataset
prophasm -k 31 -i inset1.fa -i inset2.fa -o outset1.fa outset2.fa  # compute simplitigs for two datasets
prophasm -k 31 -i inset1.fa -i inset2.fa -x intersect.fa -o outset1.fa outset2.fa  # compute simplitigs for two datasets and subtract their intersection

```


## Command line parameters

<!---
USAGE-BEGIN
-->
```
Program:  prophasm (a greedy assembler for k-mer set compression)
Version:  0.1.1
Contact:  Karel Brinda <karel.brinda@hms.harvard.edu>

Usage:    prophasm [options]

Examples: prophasm -k 15 -i f1.fa -i f2.fa -x fx.fa
             - compute intersection of f1 and f2
          prophasm -k 15 -i f1.fa -i f2.fa -x fx.fa -o g1.fa -o g2.fa
             - compute intersection of f1 and f2, and subtract it from them
          prophasm -k 15 -i f1.fa -o g1.fa
             - re-assemble f1 to g1

Command-line parameters:
 -k INT   K-mer size.
 -i FILE  Input FASTA file (can be used multiple times).
 -o FILE  Output FASTA file (if used, must be used as many times as -i).
 -x FILE  Compute intersection, subtract it, save it.
 -s FILE  Output file with k-mer statistics.
 -S       Silent mode.

Note that '-' can be used for standard input/output.

```
<!---
USAGE-END
-->


## Algorithm

```python
def extend_simplitig_forward (K, simplitig):
	extending = True
	while extending:
		extending = False
		q = simplitig[-k+1:]
		for x in [‘A’, ‘C’, ‘G’, ‘T’]:
			kmer = q + x
			if kmer in K:
				extending = True
				simplitig = simplitig + x
				S.remove (kmer)
				S.remove (reverse_complement (kmer))
				break
	return S, s

def get_maximal_simplitig (K, initial_kmer):
	simplitig = initial_kmer
	K.remove (initial_kmer)
	K.remove (reverse_completement (initial_kmer))
	K, simplitig = extend_simplitig_forward (K, simplitig)
	simplitig = reverse_completent (simplitig)
	K, simplitig = extend_simplitig_forward (K, simplitig)
	return K, simplitig

def compute_simplitigs (kmers):
	K = set()
	for kmer in kmers:
		K.add (kmer)
		K.add (reverse_completement(kmer))
	simplitigs = set()
	while |K|>0:
		initial_kmer = K.random()
		K, simplitig = get_maximal_simplitig (K, initial_kmer)
		simplitigs.add (simplitig)
	return simplitigs
```

<!--
<img alt="Greedy assembly" src="figures/greedy_assembly.png" height="150px" width="540px" /><img alt="Subtraction of k-mer sets" src="figures/subtraction.png" height="180px" width="355px" />
-->

## Similar programs

* [BCalm 2](https://github.com/GATB/bcalm) - Computation of unitigs.
* [Unikmer](https://github.com/shenwei356/unikmer) - K-mer set operations.


## Issues

Please use [Github issues](https://github.com/prophyle/prophasm/issues).


## Changelog

See [Releases](https://github.com/prophyle/prophasm/releases).


## Licence

[MIT](https://github.com/prophyle/prophasm/blob/master/LICENSE)


## Author

Karel Brinda \<karel.brinda@hms.harvard.edu\>
