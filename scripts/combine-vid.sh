#!/bin/bash

ffmpeg -i brute_ani.mp4 -i bh_ani.mp4 -i fmm_ani.mp4 -filter_complex hstack=inputs=3 combined.mp4
