.PHONY: all help clean test prophasm readme

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all: prophasm

prophasm:
	$(MAKE) -C src

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

test:
	$(MAKE) -C tests

readme:
	f=$$(mktemp);\
	  echo $$f;\
	  sed '/USAGE-BEGIN/q' README.md >> $$f; \
	  printf -- '-->\n```' >> $$f; \
	  ./prophasm -h 2>&1 | perl -pe 's/^(.*)$$/\1/g' >> $$f; \
	  printf '```\n<!---\n' >> $$f; \
	  sed -n '/USAGE-END/,$$ p' README.md >> $$f;\
	  cat $$f \
	  | perl -pe 's/^[\s]+$$/\n/g' \
	  | perl -pe 's/[\s]+$$/\n/g' \
	  > README.md
	markdown_py README.md > README.html


clean: ## Clean
	$(MAKE) -C src clean
	$(MAKE) -C tests clean
	rm -f prophasm README.html
