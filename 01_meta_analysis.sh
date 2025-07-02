######################################################
###                                                ###
### **Script**  01_meta_analysis.sh                ###
###                                                ###
### **Project** Genome-wide association study of   ### 
###             REM-sleep behaviour disorder       ###
###             identifies new risk loci           ###
###                                                ###
### **Last update** June, 2025                     ###
###                                                ###
######################################################


# Meta-analysis of GWAS summary statistics
## Using all-RBD analysis as example


## ---------------- Analysis Start --------------- ##

###### Run METAL ######
## Metal code to perform meta-analysis of fullSTATS 
module load metal
metal 
    SCHEME STDERR
    FREQLABEL freq
    AVERAGEFREQ ON
    MINMAXFREQ ON
    MARKER SNP
    ALLELE A1 A2
    EFFECT beta
    STDERR se
    PVALUE p
    WEIGHT N
    PROCESS META_sumstats_LK_flipped.txt
    PROCESS omnix_sumstats_flipped.tab 
    PROCESS nba_sumstats_flipped.tab
    OUTFILE temp .tbl
    ANALYZE HETEROGENEITY
    QUIT    


###### Clean up metal output ######
## Extract chromosome and position from SNP ID in the first column (format: Chr:Pos:A1:A2)
awk 'BEGIN{FS=OFS="\t"}{split($1,snp,":"); print snp[1],snp[2]}' temp1.tbl > chrbp.txt
sed 's/chr//g' chrbp.txt > chrbp2.txt
paste chrbp2.txt temp1.tbl > results/META_nba_omnix.tbl 
mv temp1.tbl.info results/META_nba_omnix.tbl.info 

###### Generate summary statistics and plots ######
cp results/META_nba_omnix.tbl meta.tab

## Run script sumstats_from_assoc.R
R < ~/runs/emsom/scripts/sumstats_from_meta.R --no-save

## Save plots
mv QQ.meta.tiff results/QQ_meta_META_nba_omnix.tiff
mv ManH.meta.tiff results/ManH_meta_META_nba_omnix.tiff

## Save summary statistics as fullSTATS
mv fullSTATS.meta.tab results/fullSTATS_meta_META_nba_omnix.txt


###### Annotate ######
head -n1 results/fullSTATS_meta_META_nba_omnix.txt > header.txt 

## For annovar compatability from full summary statistics keep columns: Chr, Pos, Pos, Ref, Alt
awk '{print $1,$2,$2,$5,$4}' results/fullSTATS_meta_META_nba_omnix.txt | sed '1d' > input.annovar2.txt

## Run ANNOVAR
perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar2.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg38 -out annotatedtable2 -remove -protocol refGene,avsnp150,clinvar_20210501,gnomad_genome  \
    -operation g,f,f,f -nastring . -polish
paste results/fullSTATS_meta_META_nba_omnix.txt annotatedtable2.hg38_multianno.txt > temp.txt
tr ' ' '\t' < temp.txt > temp2.txt
sort -gk11 temp2.txt > results/ANNO_meta_META_nba_omnix.txt

## Remove files we won't need
rm header.txt 
rm annotatedtable.hg38_multianno.txt 
rm input.annovar.txt

## ---------------- Analysis End ----------------- ##
