#!/bin/bash -ue
tracer start
tracer log "Started Tracer run"
tracer tool fastqc 0.12.1
fastqc s1_1.fq s1_2.fq -o .
