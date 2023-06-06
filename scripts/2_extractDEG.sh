#!bin/bash

# regulated gene
echo -e 'project\tRAD51C(-1)\tRAD51C(0)\tRAD51C(1)' > RAD51C_regulated.temp 
awk -v lines=3 '/RAD51C regulated/ {for(i=lines;i;--i) getline; split(FILENAME, file, "."); split(file[1], project, "_");printf "%s_%s\t%s\n",project[2],project[3], $0;}' mRNA_*_meth.o* | column -t >> RAD51C_regulated.temp

echo -e 'project\tBRCA1(-1)\tBRCA1(0)\tBRCA1(1)' > BRCA1_regulated.temp 
awk -v lines=3 '/BRCA1 regulated/ {for(i=lines;i;--i) getline; split(FILENAME, file, "."); split(file[1], project, "_");printf "%s_%s\t%s\n",project[2],project[3], $0;}' mRNA_*_meth.o* | column -t >> BRCA1_regulated.temp

# silencing from DGE
echo -e 'project\tgene\tlogFC\tlogCPM\tF\tPValue\tFDR' > silence.temp 
awk -v lines=3 '/BRCA1 p value/ {for(i=lines;i;--i) getline; split(FILENAME, file, "."); split(file[1], project, "_");printf "%s_%s\t%s\n",project[2],project[3], $0;}' mRNA_*_meth.o* | column -t >> data/processed/mRNA/silence_DEG.tsv

# merge files
awk -F, 'FNR==1{if (NR==1) printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\n","project","RAD51C(-1)","RAD51C(0)","RAD51C(1)","BRCA1(-1)", "BRCA1(0)", "BRCA1(1)";next} FNR==NR{a[$1]=$0;next} {print $0 "\t" (($1 in a)? a[$1]:"NA\tNA\tNA"); delete a[$1]}END{for (i in a) split(a[i], parts, "\t"); print parts[1]"\tNA\tNA\tNA\t"parts[2]"\t"parts[3]"\t"parts[3]}' RAD51C_regulated.temp BRCA1_regulated.temp | column -t > data/processed/mRNA/DGE.tsv

# clean temp
rm *.temp
