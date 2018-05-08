.PHONY: all help clean test prophasm

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all: prophasm

prophasm:
	$(MAKE) -C src

help: ## Print help message
	    @echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

test:
	$(MAKE) -C tests

clean: ## Clean
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	rm -f prophasm
