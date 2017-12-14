SHELL=/bin/bash -e -o pipefail

.SECONDARY:

IND=../../prophyle_index
F2K=../fa_to_kmers.py

FQ=../simulation_bacteria.1000.fq
FA=index.fa

d=$(shell pwd)
$(info )
$(info Directory: $(d))
$(info )
