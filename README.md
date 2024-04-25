# Hierarchical (Overview) {#mainpage}

_Hierarchical_ contains my attempt to implement serial versions of hierarchical
N-body simulation algorithms in C++. This includes the Barnes-Hut(BH) method and
the Fast Multipole Method (FMM). A brute-force algorithm is also included for
comparison.

## Directory structure

For explanation on the role of each file, see **File Index** in the
Doxygen-generated documentation.

- `src`: Where the main body of the project is stored (C++).
- `tests`: C++ test files, depends on `doctest`.
- `build`: Created on calling `make`; Where binaries and object files are
  stored.
- `pysrc`: Where all Python source code for data processing is stored.
- `data`: Where pre-written experiments dump their test results. This folder is
  expected to be present for them, or else they will just abort.
- `docs`: Documentation, `Doxyfile`.
- `notes`: Contains report and some other notes in markdown format.
- `scripts`: Very short bash scripts for small tasks.
- `makefiles`: Stores all `make` targets.

## Build Instructions

### `hchl`: Main C++ binary

Assuming one is on a Unix-like system, simply navigate to the project root, then
type:

~~~{.bash}
make
~~~

This requires that one has `GNU Make` and `clang` installed. Alternatively to
use `gcc`, replace occurrences of `clang++` with `g++` for all files in the
`makefile` folder.

This will compile all files in `src` and `main.cpp`, which has the `main`
function, the entry point of the program. Uncomment some of the lines in `main`
then build to see the effect of individual tests.

After compilation finishes, the binary can then be found at `build/hchl`. When
benchmarking, I run `./run` to set CPU affinities.

### `test`: C++ tests with `doctest`

I have written some tests for parts of the `C++` source. In order to compile
this, make sure the `doctest` header is installed. On my system (Arch Linux), I
can then write in the test files

~~~{.cpp}
#include "doctest/doctest.h"
~~~

This may have to be replaced with your own library location, which one specify
by adding `-Ipath/to/doctest/lib` in the `CXX_FLAGS` variable in
`makefile.test`.

To then compile the test run

~~~{.bash}
make test
~~~

The binary will be found at `build/test`.

Note that this will be compiled with `-O0` flag and is incompatible with
previous object files from the `hchl` target compiled with `-O3`. Therefore
`make clean` (see below) must be run first before changing C++ target.

### `report`: PDF Generation from markdown

To compile my report, written in `pandoc markdown` at `notes/report.md`, call

~~~{.bash}
make report
~~~

The compiled report will be found at `notes/report.pdf`. All other markdown
files found in `notes/` will also be compiled. Requires `pandoc` and the
`pandoc-crossref` filter.

### `doc`: `Doxygen` documentation generation

I have written `Doxygen`-compatible docstrings in both C++ and Python code. To
make the documentation with `doxygen`, run

~~~{.bash}
make doc
~~~

This calls `doxygen` with file options set in `docs/Doxyfile`. For `html`
version open `docs/html/index.html` with a browser. For `pdf`, open
`docs/latex/refman.pdf`. Requires `doxygen`. PDF compilation requires
appropriate `latex` packages.

### `clean`: Directory clean-up

Cleans all files resulting from compilation and `make` commands. Please note it
will also **remove .plist files in the project root**, since my language servers
are generating a lot of garbage in the format.

### `Python`

To run Python code, simply edit the `main` function of `main.py` and execute the
file. Comment and uncomment lines to see the working of data processing
functions, provided that data has been generated previously from `C++`, stored
to `data/`, and that the script is executed at the project root.

Some common python libraries such as `matplotlib` and `scipy` will be required.
Additionally, parts of the code (animation) depends on `ffmpeg` and a Unix-like
shell environment.

## Coding style

The coding style of the project largely follows that outlined in Section 6.6 of
_Guide to Scientific Computing in C++, Second Edition_ by Pitt-Francis &
Whiteley.
