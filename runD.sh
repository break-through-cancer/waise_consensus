cd /pn/ct
rm -rf outD
nextflow run /wc/repo/main.nf --entry cohort --csv /pn/ct/sheetD.csv --publishDir /pn/ct/outD -w /pn/ct/workB -resume -with-docker -process.containerOptions='--platform linux/amd64' < /dev/null > logD.txt 2>&1
echo done > doneD.txt
