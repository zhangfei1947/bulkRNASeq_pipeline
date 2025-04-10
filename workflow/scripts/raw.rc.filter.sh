#keep if at least one sample got >10 reads
#awk 'BEGIN NR==2{print $0; next} NR>2{for(i=7;i<=NF;i++) if($i>10){print $0; next}}' {input} > {output}

cut -f1,7- ${snakemake_input[0]} | awk '
BEGIN{FS="\t"}
NR==1{next}
NR==2{
    printf "Geneid"
    for(i=2;i<=NF;i++){
        split($i,a,"/")
        printf "\t%s", a[2]
    }
    printf "\n"
    next
}
NR>2{
    for(i=2;i<=NF;i++) if($i>10){print $0; next}
}' > ${snakemake_output[0]}

#keep if avg(reads) > 1
#awk 'BEGIN NR==2{print $0; next} NR>2{sum=0; for(i=7;i<=NF;i++) sum+=$i; if(sum>(NF-6)) print $0}'