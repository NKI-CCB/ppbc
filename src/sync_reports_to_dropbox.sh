#!/bin/bash

#Must be run from local computer (with rclone installed) after mounting postpartum drive
#Do not ssh into harris first

#From local:
#ssh -f k.moore@rhpc.nki.nl -L 2222:harris:22 -N && sshfs -p 2222 k.moore@localhost:/DATA/share/postpartumbc ~/mnt/postpartumbc/

#Html reports
cd ~/mnt/postpartumbc/reports
for REPORT in `ls *nb.html`; do echo $REPORT; rclone sync -P $REPORT remote:ppbc/reports; done

#Excel differential expression lists
cd ~/mnt/postpartumbc/results/diffex
for EXCEL in `ls *.xlsx`; do echo $EXCEL; rclone sync -P $EXCEL remote:ppbc/results/diffex/excel; done
rclone sync -P figs remote:ppbc/results/diffex/figs

cd ~/mnt/postpartumbc/results/
rclone sync -P cibersortX remote:ppbc/results/cibersortX
rclone sync -P clustering remote:ppbc/results/clustering