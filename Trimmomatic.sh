#!/bin/bash

for FWD in ../*R1_Combined.fq
do
    REV=${FWD/R1/R2}
    OutFWDpe=`basename $FWD`
    OutFWDse=`basename ${FWD/fq/UNPAIRED.fq}`
    OutREVpe=`basename $REV`
    OutREVse=`basename ${REV/fq/UNPAIRED.fq}`

    TrimmomaticPE -threads 8 -phred33 \
        $FWD \
        $REV \
        $OutFWDpe \
        $OutFWDse \
        $OutREVpe \
        $OutREVse \
        ILLUMINACLIP:TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:35
done
