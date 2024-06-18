#!/bin/bash -ue
tracer start
tracer log "Tracer run initialized..."
tracer tool index v1
salmon index --threads 4 -t ref.fa -i salmon_index
