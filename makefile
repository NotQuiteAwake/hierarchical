MAKEFLAGS := --jobs=$(shell nproc)
# MAKEFLAGS += --output-sync=target

all: hchl hchl_test

hchl:
	make -f makefiles/makefile.hchl

hchl_test:
	make -f makefiles/makefile.test

clean:
	make -f makefiles/makefile.hchl clean
	make -f makefiles/makefile.test clean
	rm *.plist
