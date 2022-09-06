smoove call \
        --outdir "/output/results-smoove/" \
        -p $1 \
        --duphold \
        -x --name $2 \
        --fasta "/genomedir/$3" \
        --genotype \
        /bamdir/*.dupremoved.sorted.bam
