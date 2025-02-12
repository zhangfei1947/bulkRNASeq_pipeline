#!/bin/bash

echo -e "Sample\tReads_brfore\tReads_after\tQ20_R1_before\tQ20_R2_before\tQ20_R1_after\tQ20_R2_after\tduplication_rate\tinsert_size" > $1
for f in {input}; do
    sample=$(basename $f .json)
    total_reads_b=$(jq '.summary.before_filtering.total_reads' $f)
    total_reads_a=$(jq '.summary.after_filtering.total_reads' $f)
    q20b_r1=$(jq '.summary.before_filtering.q20_rate' $f)
    q20b_r2=$(jq '.summary.before_filtering.q20_rate' $f)
    q20a_r1=$(jq '.summary.after_filtering.q20_rate' $f)
    q20a_r2=$(jq '.summary.after_filtering.q20_rate' $f)
    dup_rate=$(jq '.summary.duplication.rate' $f)
    insert_size=$(jq '.summary.insert_size.peak' $f)
    echo -e "$sample\t$total_reads_b\t$total_reads_a\t$q20b_r1\t$q20b_r2\t$q20a_r1\t$q20a_r2\t$dup_rate\t$insert_size" >> $1
done