cd /pn/ct
rm -rf outD outA3
nextflow run /wc/repo/main.nf --entry cohort --csv /pn/ct/sheetD.csv --publishDir /pn/ct/outD -w /pn/ct/workB -resume -with-docker -process.containerOptions='--platform linux/amd64' < /dev/null > logD.txt 2>&1
nextflow run /wc/repo/main.nf --entry cohort --csv /pn/ct/sheetA.csv --matched_normal_subtract false --normal_alt_only false --publishDir /pn/ct/outA3 -w /pn/ct/workB -resume -with-docker -process.containerOptions='--platform linux/amd64' < /dev/null > logA3.txt 2>&1
echo done > doneDA.txt
