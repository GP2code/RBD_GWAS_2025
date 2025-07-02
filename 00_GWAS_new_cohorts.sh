######################################################
###                                                ###
### **Script**  00_GWAS_new_cohorts.sh             ###
###                                                ###
### **Project** Genome-wide association study of   ### 
###             REM-sleep behaviour disorder       ###
###             identifies new risk loci           ###
###                                                ###
### **Last update** June, 2025                     ###
###                                                ###
######################################################


# GWASs of new cases/controls  
### Performed for NeuroBooster and OmniExpress genotype datasets

## ===================== NBA ===================== ##

## ---------------- Analysis Start --------------- ##

## Regression
plink --bfile nba_topmed_sc_updated --logistic beta --maf 0.01 --ci .95 --covar covar.txt --covar-name Age,Sex,PC1,PC2,PC3,PC4,PC5 --out nba_age_sex_pcs

## Keep rows where P-value is not NA & test is 'ADD' (Additive) & sort P-values in ascending order
awk '{if ($12!="NA") print $0}' nba_age_sex_pcs.assoc.logistic | grep 'ADD' | sort -gk12 > logistic/p.nba_age_sex_pcs.assoc.logistic

## Manually add rare variant regression results
plink --bfile nba_topmed_sc_updated --extract extract_rare_vars.txt --logistic beta --ci .95 --covar covar.txt --covar-name Age,Sex,PC1,PC2,PC3,PC4,PC5 --out nba_age_sex_pcs_GBA

## Keep rows where P-value is not NA & test is 'ADD' (Additive) & sort P-values in ascending order
awk '{if ($12!="NA") print $0}' nba_age_sex_pcs_GBA.assoc.logistic | grep 'ADD' | sort -gk12 > logistic/p.nba_age_sex_pcs_GBA.assoc.logistic
cat logistic/p.nba_age_sex_pcs_GBA.assoc.logistic logistic/p.nba_age_sex_pcs.assoc.logistic > temp 
mv temp logistic/p.nba_age_sex_pcs.assoc.logistic

## Generate summary statistics and plots 
echo nba_age_sex_PCs > marker.txt 
cp nba.info info
cp logistic/p.nba_age_sex_pcs.assoc.logistic assoc

## Run script sumstats_from_assoc.R
R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
paste lambda.txt >> lambdas.txt

## Save plots
mv QQ.tiff plots/QQ_nba_age_sex_pcs.tiff
mv ManH.tiff plots/ManH_nba_age_sex_pcs.tiff
## Save summary statistics as fullSTATS
mv Metal.tab sumstats/Metal_nba_age_sex_pcs.txt
mv COJO.tab sumstats/COJO_nba_age_sex_pcs.txt
mv fullSTATS.tab sumstats/fullSTATS_nba_age_sex_pcs.txt

## Annotate 
## For annovar compatability from full summary statistics keep columns Chr, Pos, Pos, Ref, Alt
awk '{print $1,$2,$2,$5,$4}' sumstats/fullSTATS_nba_age_sex_pcs.txt | sed '1d' > input.annovar.txt
## Run ANNOVAR
perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210501,gnomad_genome \
    -operation g,f,f,f -nastring . -polish
paste sumstats/fullSTATS_nba_age_sex_pcs.txt annotatedtable.hg38_multianno.txt > temp.txt
tr ' ' '\t' < temp.txt > temp2.txt
sort -gk10 temp2.txt > annos/ANNO_nba_age_sex_pcs.txt

## ---------------- Analysis End ----------------- ##




## ==================== Omnix ==================== ##

## ---------------- Analysis Start --------------- ##


###### Regression ######
plink --bfile omnix_topmed_sc_updated --logistic beta --maf 0.01 --ci .95 --covar covar.txt --covar-name Age,Sex,PC1,PC2,PC3,PC4,PC5 --out omnix_age_sex_PCs

## Keep rows where P-value is not NA & test is 'ADD' (Additive) & sort P-values in ascending order
awk '{if ($12!="NA") print $0}' omnix_age_sex_PCs.assoc.logistic | grep 'ADD' | sort -gk12 > logistic/p.omnix_age_sex_PCs.assoc.logistic

## Manually add rare variant regression results
plink --bfile omnix_topmed_sc_updated --extract extract_rare_vars.txt --logistic beta --ci .95 --covar covar.txt --covar-name Age,Sex,PC1,PC2,PC3,PC4,PC5 --out omnix_age_sex_pcs_GBA

## Keep rows where P-value is not NA & test is 'ADD' (Additive) & sort P-values in ascending order
awk '{if ($12!="NA") print $0}' omnix_age_sex_pcs_GBA.assoc.logistic | grep 'ADD' | sort -gk12 > logistic/p.omnix_age_sex_pcs_GBA.assoc.logistic
cat logistic/p.omnix_age_sex_pcs_GBA.assoc.logistic logistic/p.omnix_age_sex_PCs.assoc.logistic > temp 
mv temp logistic/p.omnix_age_sex_pcs.assoc.logistic


###### Generate summary statistics and plots ######
echo omnix_age_sex_PCs > marker.txt 
cp omnix.info info
cp logistic/p.omnix_age_sex_pcs.assoc.logistic assoc

## Run script sumstats_from_assoc.R
R < ~/runs/emsom/scripts/sumstats_from_assoc.R --no-save
paste lambda.txt >> lambdas.txt

## Save plots
mv QQ.tiff plots/QQ_omnix_age_sex_pcs.tiff
mv ManH.tiff plots/ManH_omnix_age_sex_pcs.tiff
## Save summary statistics as fullSTATS
mv Metal.tab sumstats/Metal_omnix_age_sex_pcs.txt
mv COJO.tab sumstats/COJO_omnix_age_sex_pcs.txt
mv fullSTATS.tab sumstats/fullSTATS_omnix_age_sex_pcs.txt


###### Annotate ######
## For annovar compatability from full summary statistics keep columns: Chr, Pos, Pos, Ref, Alt
awk '{print $1,$2,$2,$5,$4}' sumstats/fullSTATS_omnix_age_sex_pcs.txt | sed '1d' > input.annovar.txt
## Run ANNOVAR
perl ~/runs/emsom/softwares/annovar/table_annovar.pl input.annovar.txt ~/runs/emsom/softwares/annovar/humandb/ \
    -buildver hg38 -out annotatedtable -remove -protocol refGene,avsnp150,clinvar_20210501,gnomad_genome \
    -operation g,f,f,f -nastring . -polish
paste sumstats/fullSTATS_omnix_age_sex_pcs.txt annotatedtable.hg38_multianno.txt > temp.txt
tr ' ' '\t' < temp.txt > temp2.txt
sort -gk10 temp2.txt > annos/ANNO_omnix_age_sex_pcs.txt

## ---------------- Analysis End ----------------- ##

