.PHONY: all help clean

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

.SUFFIXES:

all: _intersect.fa
	../tools/verify_output.py test1.fa test2.fa _out1.fa _out2.fa _intersect.fa 22

_intersect.fa:
	../../prophasm -i test1.fa -i test2.fa -o _out1.fa -o _out2.fa -x _intersect.fa -k 22

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

clean: ## Clean
	rm -f _*.fa

