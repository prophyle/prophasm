TESTS = $(wildcard */.)

.PHONY: all help clean $(TESTS) parallel

SHELL=/usr/bin/env bash -eo pipefail

.SECONDARY:

all: ## Run tests
all: $(TESTS)

help: ## Print help message
	@echo "$$(grep -hE '^\S+:.*##' $(MAKEFILE_LIST) | sed -e 's/:.*##\s*/:/' -e 's/^\(.\+\):\(.*\)/\\x1b[36m\1\\x1b[m:\2/' | column -c2 -t -s : | sort)"

V=1

parallel: ## Run tests in parallel
	$(MAKE) -j

$(TESTS): compile
	@echo
	@echo "========================================="
	@echo "$@"
	@echo "========================================="
	@echo
	$(MAKE) -C $@

compile: ## Compile prophyle_index
	$(MAKE) -C ../src

clean: ## Clean
	@for a in $(TESTS); do \
		$(MAKE) -C $$a clean; \
	done;
