#!/bin/bash -ue
tracer tool macs 3.0.1
macs3 callpeak -t P1s1.sorted.bam -f BAMPE -p 0.05 --outdir ./P1s1_peaks
