#!/bin/bash

# some outputs didn't have number of repeats (10) on top.
# this scripts manually add it to them to make Python's life easier.

for dump in *.dump; do
    sed -i '1i\10' "$dump"
done
