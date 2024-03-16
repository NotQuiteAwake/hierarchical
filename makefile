MAKEFLAGS := --jobs=$(shell nproc)
# MAKEFLAGS += --output-sync=target

.PHONY: default all hchl test docs clean

default: hchl

all: hchl test docs

hchl:
	make -f makefiles/makefile.hchl

test:
	make -f makefiles/makefile.test

docs:
	make -f makefiles/makefile.docs

clean:
	make -f makefiles/makefile.hchl clean
	make -f makefiles/makefile.test clean
	make -f makefiles/makefile.docs clean
	rm *.plist
