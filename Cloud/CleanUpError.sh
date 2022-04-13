#!/bin/bash

Analysis=Excord_Batch1
mkdir ${Analysis}_OutputError
mv *error ${Analysis}_OutputError
mv *out ${Analysis}_OutputError
tar -zcvf ${Analysis}_OutputError.tar.gz ${Analysis}_OutputError

aws s3 cp ${Analysis}_OutputError.tar.gz s3://thousandgenometriovariants.store.ubc-hpc.cloud  

