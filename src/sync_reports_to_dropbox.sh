#!/bin/bash

#Must be run from local computer (with rclone installed) after mounting postpartum drive
#Do not ssh into harris first

#From local:
sshfs harris:/DATA/share/postpartumbc ~/mnt/postpartumbc

#Html reports
cd ~/mnt/postpartumbc/reports
for REPORT in `ls *.html`; do echo $REPORT; rclone sync -P $REPORT remote:ppbc/reports; done

#Metadata
cd ~/mnt/postpartumbc/data/metadata
for f in `ls | grep -E ".tsv|.csv|.txt|.pdf|.xlsx"`; do echo $f; rclone sync -P $f remote:ppbc/reports; done

#Excel differential expression lists
cd ~/mnt/postpartumbc/results/diffex
for EXCEL in `ls *.xlsx`; do echo $EXCEL; rclone sync -P $EXCEL remote:ppbc/results/diffex/excel; done
rclone sync -P figs remote:ppbc/results/diffex/figs

#Results directories
cd ~/mnt/postpartumbc/results/
rclone sync -P cibersortX remote:ppbc/results/cibersortX
rclone sync -P clustering remote:ppbc/results/clustering
rclone sync -P survival remote:ppbc/results/survival

cd ~/mnt/postpartumbc/results/flexgsea/deseq
rclone sync -P *.xlsx remote:ppbc/results/flexgsea

