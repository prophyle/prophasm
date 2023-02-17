# ProphAsm

[![Build Status](https://travis-ci.org/prophyle/prophasm.svg?branch=master)](https://travis-ci.org/prophyle/prophasm)

<!-- vim-markdown-toc GFM -->

* [Introduction](#introduction)
* [Cite](#cite)
* [Prerequisities](#prerequisities)
* [Getting started](#getting-started)
* [How to use](#how-to-use)
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

To cite the concept of simplitigs and _ProphAsm_ as a tool, please use the following reference:

> Břinda K, Baym M, Kucherov G.
> **Simplitigs as an efficient and scalable representation of de Bruijn graphs.**
*Genome Biology* **22**(96), 2021; doi: https://doi.org/10.1186/s13059-021-02297-z

```
@article{brinda2021-simplitigs,
  title   = { Simplitigs as an efficient and scalable representation of de {Bruijn} graphs },
  author  = { Karel B{\v r}inda and Michael Baym and Gregory Kucherov },
  journal = { Genome Biology },
  volume  = { 22 },
  number  = { 96 },
  year    = { 2021 },
  doi     = { 10.1186/s13059-021-02297-z }
}
```


For the concept of simplitigs, you might also consider citing the following paper from another group,
introducing independently the same concept under the name
_spectrum-preserving string sets_ (SPSS):

> Rahman A and Medvedev P. **Representation of k-mer sets using
  spectrum-preserving string sets.** *Journal of Computational Biology* **28**(4), pp. 381-394, 2021. https://doi.org/10.1089/cmb.2020.0431


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


## How to use

Set operations:
```
./prophasm -k 15 -i tests/test1.fa -i tests/test2.fa -o _out1.fa -o _out2.fa -x _intersect.fa -s _stats.tsv
   ```


## Command-line arguments

<!---
USAGE-BEGIN
-->
```
Program:  prophasm (computation of simplitigs and k-mer set operations)
Version:  0.1.2
Contact:  Karel Brinda <karel.brinda@inria.fr>

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
				K.remove (kmer)
				K.remove (reverse_complement (kmer))
				break
	return K, simplitig

def get_maximal_simplitig (K, initial_kmer):
	simplitig = initial_kmer
	K.remove (initial_kmer)
	K.remove (reverse_complement (initial_kmer))
	K, simplitig = extend_simplitig_forward (K, simplitig)
	simplitig = reverse_complement (simplitig)
	K, simplitig = extend_simplitig_forward (K, simplitig)
	return K, simplitig

def compute_simplitigs (kmers):
	K = set()
	for kmer in kmers:
		K.add (kmer)
		K.add (reverse_complement(kmer))
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
