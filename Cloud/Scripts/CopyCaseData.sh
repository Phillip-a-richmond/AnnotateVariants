#!/bin/bash
Case=Case1
S3Bucket=genomeanalysistraining.store.ubc-hpc.cloud
mkdir -p $Case
aws s3 cp --recursive \
	s3://${S3Bucket}/CaseFiles/${Case} \
	/shared/$Case \
	--dryrun

# Once validated the above, you can run this
aws s3 cp --recursive \
        s3://${S3Bucket}/CaseFiles/${Case} \
        /shared/$Case 

