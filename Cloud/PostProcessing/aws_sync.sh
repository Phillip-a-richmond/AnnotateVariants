aws s3 sync \
        ./ \
        s3://thousandgenometriovariants.store.ubc-hpc.cloud \
        --include '*' --exclude '' \
        --dryrun

aws s3 sync \
        ./ \
        s3://thousandgenometriovariants.store.ubc-hpc.cloud \
        --include '*' --exclude '' 
