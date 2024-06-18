#!/bin/bash -ue
printf "%s %s\n" gut_1.fq gut.fq | while read old_name new_name; do
    [ -f "${new_name}" ] || ln -s $old_name $new_name
done

fastqc \
     \
    --threads 1 \
    gut.fq

cat <<-END_VERSIONS > versions.yml
"QC":
    fastqc: $( fastqc --version | sed '/FastQC v/!d; s/.*v//' )
END_VERSIONS
