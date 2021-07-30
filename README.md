# ProphAsm

[![Build Status](https://travis-ci.org/prophyle/prophasm.svg?branch=master)](https://travis-ci.org/prophyle/prophasm)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Cite](#cite)
* [Prerequisities](#prerequisities)
* [Getting started](#getting-started)
* [Links](#links)
* [Issues](#issues)
* [Changelog](#changelog)
* [Licence](#licence)
* [Author](#author)

<!-- vim-markdown-toc -->

## Introduction

ProphAsm is a tool for computing *simplitigs* from *k-mer sets*. Simplitigs are
strings obtained as disjoint paths in a bidirectional vertex-centric
de Bruijn graph. Compared to unitigs, simplitigs provide an improvement in the
number of sequences and their cumulative length, while both representations
carry the same k-mers. For more details, see the
[paper](https://doi.org/10.1101/2020.01.12.903443).

Various types of sequencing datasets can be used as the input for
ProphAsm, including genomes, pan-genomes, metagenomes or sequencing reads.
Besides computing simplitigs, ProphAsm can also compute intersection
and set differences of k-mer
sets (while set unions are easy to compute simply by merging the source files).

Upon execution, ProphAsm first loads all specified datasets (see the `-i`
param) and the corresponding k-mer sets (see the `-k` param). If the `-x` param
is provided, ProphAsm then computes their intersection, subtracts the
intersection from the individual k-mer sets and computes simplitigs for the
intersection. If output files are specified (see the `-o` param), it computes
also set differences.


## Cite

To cite _ProphAsm_ as a tool, please use the following reference:

> Brinda K, Baym M, and Kucherov G. **Simplitigs as an efficient and scalable representation of de Bruijn graphs**
. bioRxiv 2020.01.12.903443, 2020. https://doi.org/10.1101/2020.01.12.903443

For _the concept of simplitigs_, please use the reference above and include
also the following paper, introducing independently the same concept under the
name spectrum-preserving string sets:

> Rahman A and Medvedev P. **Representation of k-mer sets using
  spectrum-preserving string sets**. Research in Computational Molecular
  Biology - 24th Annual International Conference, {RECOMB} 2020, Padua, Italy,
  May 10-13, 2020, Proceedings. Lecture Notes in Computer Science 12074, pp.
  152-168, Springer, 2020. https://doi.org/10.1007/978-3-030-45257-5_10


## Prerequisities

* GCC 4.8+ or equivalent
* ZLib


## Getting started

Download and compile ProphAsm:

```
git clone https://github.com/prophyle/prophasm
cd prophasm && make -j
```

Compute simplitigs:

```
./prophasm -k 15 -i tests/test1.fa -o simplitigs.fa
```

Set operations:
```
./prophasm -k 15 -i tests/test1.fa -i tests/test2.fa -o _out1.fa -o _out2.fa -x _intersect.fa -s _stats.tsv
   ```


## Command line parameters

<!---
USAGE-BEGIN
-->
```
Program:  prophasm (computation of simplitigs and k-mer set operations)
Version:  0.1.2
Contact:  Karel Brinda <karel.brinda@hms.harvard.edu>

Usage:    prophasm [options]

Examples: prophasm -k 31 -i ref.fa -o simplitigs.fa
           - compute simplitigs of ref.fa
          prophasm -k 31 -i ref1.fa -i ref2.fa -x inter.fa
           - intersect the k-mers sets of ref1 and ref2
          prophasm -k 31 -i ref1.fa -i ref2.fa -x inter.fa -o dif1.fa -o dif2.fa
           - intersect ref1 and ref2, and compute the set differences

Command-line parameters:
 -k INT   k-mer length (from [1, 32])
 -i FILE  input FASTA file (can be used multiple times)
 -o FILE  output FASTA file (if used, must be used as many times as -i)
 -x FILE  compute intersection, subtract it, save it
 -s FILE  output file with k-mer statistics
 -S       silent mode

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

## Links

* [Sneak peek at the -tigs!](https://kamimrcht.github.io/webpage/tigs.html) - An overview of different *tigs in computational biology.
* [UST](https://github.com/medvedevgroup/UST/) - Another tool for computing simplitigs. Unlike ProphAsm, UST requires pre-computed unitigs as the input, therefore the method is overall more resource-demanding.
* [BCalm 2](https://github.com/GATB/bcalm) - The best available tool for computing unitigs.
* [Unikmer](https://github.com/shenwei356/unikmer) - Another tool for k-mer set operations.


## Issues

Please use [Github issues](https://github.com/prophyle/prophasm/issues).


## Changelog

See [Releases](https://github.com/prophyle/prophasm/releases).


## Licence

[MIT](https://github.com/prophyle/prophasm/blob/master/LICENSE)


## Author

Karel Brinda \<karel.brinda@hms.harvard.edu\>
