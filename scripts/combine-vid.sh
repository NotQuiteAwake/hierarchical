#!/bin/bash

set -e

FOLDER=${1:-"./"}

BRUTE=$FOLDER/brute_ani.mp4
BH=$FOLDER/bh_ani.mp4
FMM=$FOLDER/fmm_ani.mp4
COMB=$FOLDER/combined.mp4

ffmpeg $2 -i $BRUTE -i $BH -i $FMM -filter_complex hstack=inputs=3 $COMB
