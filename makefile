MAKEFLAGS := --jobs=$(shell nproc)
# MAKEFLAGS += --output-sync=target

MKF = $(wildcard makefiles/makefile.*)

.PHONY: clean hchl default

default: hchl

hchl:
	make -f makefiles/makefile.hchl

clean:
	make -f makefiles/makefile.hchl clean
	for mk in $(MKF); do make -f $$mk clean; done
	rm *.plist

%:
	make -f makefiles/makefile.$@



