######################################################
###                                                ###
### **Script**  02_RCOR1_in_PD_and_DLB.sh          ###
###                                                ###
### **Project** Genome-wide association study of   ### 
###             REM-sleep behaviour disorder       ###
###             identifies new risk loci           ###
###                                                ###
### **Last update** June, 2025                     ###
###                                                ###
######################################################


# Investigation of common and rare variants in synucleinopathy GWASs
### UCSC Genome Browser was used to gather information on the RCOR1 locus positions


## ---------------- Analysis Start --------------- ##

## ===== Common variants =====

# Pull RCOR1 region +/- 500kb
# GWAS_studies.txt just has the names of each external GWAS summary statistics
cat GWAS_studies.txt | while read line
do
    awk '{if ($1 == 14 && $2 > 102108461 && $2 < 103108461) print}' ${line}.txt | sort -k2,2 > stats.txt # Pull 500kb +/- region
    awk 'BEGIN {OFS="\t"} {print $1,$2,$3,$5,$4,".",".",".","."}' stats.txt > temp.vcf # VEP prep
    cat header.txt temp.vcf > input.vcf # VEP prep
    bash ~/runs/emsom/scripts/VEP_simple_hg38.sh # VEP annotation
    grep -v '#' output.vcf | awk 'BEGIN {OFS="|"} {print $3,$8}' > annos.txt # Format VEP output 
    awk -F"|" 'BEGIN {OFS="\t"} {print $1,$19,$5}' annos.txt > temp.txt # Keep SNP, rsID, and gene
    awk -F"\t" 'BEGIN {OFS="\t"} {if ($2 == "") $2 = "."; print}' temp.txt > temp2.txt # Fill blanks in rsID column
    awk 'FNR==NR {a[$1]=$2; b[$1]=$3; next} {print $0, a[$3], b[$3]}' temp2.txt stats.txt > RCOR1_500kb_${line}.txt
done

## ===== Rare variants ======

# SKAT-O burden analysis 
# Loops through gene list to apply burden testing to all 
# Utilized WGS data from AMP-PD

## Make bfiles for all genes in PD and DLB datasets
while read -r gene; do
    mkdir $gene # Make a directory for each gene

    awk -v gene="$gene" '{if ($4 == gene) print $1, $2-300, $3+300, $4}' RCOR1_genes.bed > temp.txt # pull coordinates for current gene
    
    if [ -s temp.txt ]; then
        read -r chr start end gene < temp.txt
            plink --bfile PD_for_SKATO --chr "$chr" --from-bp "$start" --to-bp "$end" --make-bed --out PD_"$gene" # pull gene from rare PD bfiles
            plink --bfile DLB_for_SKATO --chr "$chr" --from-bp "$start" --to-bp "$end" --make-bed --out DLB_"$gene" # pull gene from rare DLB bfiles
    else
        echo "No data found for gene $gene"
    fi

    mv PD_$gene.* $gene
    mv DLB_$gene.* $gene
done < gene_list.txt

##### Run ANNOVAR #####
while read -r gene; do
    plink --bfile $gene/PD_$gene --recode vcf --out $gene/PD_$gene
    plink --bfile $gene/DLB_$gene --recode vcf --out $gene/DLB_$gene   
done < gene_list.txt

## Convert to .vcf to .avinput and annotate with ANNOVAR 
while read -r gene; do
    perl ~/runs/emsom/softwares/annovar/convert2annovar.pl --format vcf4 "$gene/PD_$gene.vcf" --allsample --withfreq --outfile "temp_PD.vcf" 2>&1
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl "temp_PD.vcf" ~/runs/emsom/softwares/annovar/humandb/ --buildver hg38 --out "$gene/PD_anno_$gene" --remove --protocol refGene,avsnp150,clinvar_20221231,ljb26_all --operation g,f,f,f --nastring . 2>&1

    perl ~/runs/emsom/softwares/annovar/convert2annovar.pl --format vcf4 "$gene/DLB_$gene.vcf" --allsample --withfreq --outfile "temp_DLB.vcf" 2>&1
    perl ~/runs/emsom/softwares/annovar/table_annovar.pl "temp_DLB.vcf" ~/runs/emsom/softwares/annovar/humandb/ --buildver hg38 --out "$gene/DLB_anno_$gene" --remove --protocol refGene,avsnp150,clinvar_20221231,ljb26_all --operation g,f,f,f --nastring . 2>&1
done < gene_list.txt

while read -r gene; do
    sed '1d' $gene/PD_anno_$gene.hg38_multianno.txt > temp.txt
    awk '{print $2}' $gene/PD_$gene.bim > temp2.txt
    paste temp2.txt temp.txt > $gene/final_PD_anno_$gene.hg38_multianno.txt

    sed '1d' $gene/DLB_anno_$gene.hg38_multianno.txt > temp.txt
    awk '{print $2}' $gene/DLB_$gene.bim > temp2.txt
    paste temp2.txt temp.txt > $gene/final_DLB_anno_$gene.hg38_multianno.txt
done < gene_list.txt


##### Prep variants categories #####

# ----------- PD ----------- #
### Keep columns for SNP/rsID, Func.refGene, ExonicFunc.refGene, CLNSIG, and CADD_phred

while read -r gene; do
    awk '{print $1,$7,$10,$18,$39}' $gene/final_PD_anno_$gene.hg38_multianno.txt > temp.txt

    awk '{print "rare",$1}' temp.txt > rare.txt # Test for all rare vars
    
    awk '{if ($3=="nonsynonymous") print "NS",$1}' temp.txt > NS.txt # Test for all nonsynonymous vars

    grep splic temp.txt | awk '{print "LoF",$1}' > LoF.txt
    grep frameshift temp.txt | awk '{print "LoF",$1}' >> LoF.txt 
    grep stopgain temp.txt | awk '{print "LoF",$1}' >> LoF.txt

    awk '{print "functional",$2}' LoF.txt > temp_functional.txt 
    awk '{print "functional",$2}' NS.txt >> temp_functional.txt 
    sort temp_functional.txt | uniq > functional.txt

    awk '{if ($5>20) print "CADD",$1}' temp.txt > CADD.txt # Test for 'CADD > 20' vars

    cat rare.txt NS.txt LoF.txt functional.txt CADD.txt > $gene/PD_$gene.SETID
    sort $gene/PD_$gene.SETID | uniq > temp.txt
    mv temp.txt $gene/PD_$gene.SETID
done < gene_list.txt

# ----------- DLB ----------- #
### Keep columns for SNP/rsID, Func.refGene, ExonicFunc.refGene, CLNSIG, and CADD_phred
while read -r gene; do
    awk '{print $1,$7,$10,$18,$39,$40}' $gene/final_DLB_anno_$gene.hg38_multianno.txt > temp.txt

    awk '{print "rare",$1}' temp.txt > rare.txt # test for all rare vars
    
    awk '{if ($3=="nonsynonymous") print "NS",$1}' temp.txt > NS.txt # Test for all nonsynonymous vars

    grep splic temp.txt | awk '{print "LoF",$1}' > LoF.txt
    grep frameshift temp.txt | awk '{print "LoF",$1}' >> LoF.txt 
    grep stopgain temp.txt | awk '{print "LoF",$1}' >> LoF.txt

    awk '{print "functional",$2}' LoF.txt > temp_functional.txt 
    awk '{print "functional",$2}' NS.txt >> temp_functional.txt 
    sort temp_functional.txt | uniq > functional.txt

    awk '{if ($5>20) print "CADD",$1}' temp.txt > CADD.txt # Test for 'CADD > 20' vars

    cat rare.txt NS.txt LoF.txt functional.txt CADD.txt > $gene/DLB_$gene.SETID
    sort $gene/DLB_$gene.SETID | uniq > temp.txt
    mv temp.txt $gene/DLB_$gene.SETID
done < gene_list.txt


##### Prep files and run SKATO #####

# ----------- PD ----------- #
while read -r gene; do
cp $gene/PD_$gene.bed temp.bed
cp $gene/PD_$gene.bim temp.bim
cp $gene/PD_$gene.fam temp.fam
cp $gene/PD_$gene.SETID temp.SETID
echo "" > temp.SSD
echo "" > temp.info
echo "" > temp.results.skato
echo "" > temp.results.burden

module load StdEnv/2020
module load r/4.2.2

R < SKAT.R --no-save

mv temp.results.skato $gene/PD_$gene.results.skato
mv temp.results.burden $gene/PD_$gene.results.burden
done < gene_list.txt

# ----------- DLB ----------- #
while read -r gene; do
cp $gene/DLB_$gene.bed temp.bed
cp $gene/DLB_$gene.bim temp.bim
cp $gene/DLB_$gene.fam temp.fam
cp $gene/DLB_$gene.SETID temp.SETID
echo "" > temp.SSD
echo "" > temp.info
echo "" > temp.results.skato
echo "" > temp.results.burden

module load StdEnv/2020
module load r/4.2.2

R < SKAT.R --no-save

mv temp.results.skato $gene/DLB_$gene.results.skato
mv temp.results.burden $gene/DLB_$gene.results.burden
done < gene_list.txt


## ---------------- Analysis End ----------------- ##

