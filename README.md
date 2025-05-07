# RBD_GWAS_2025

Rapid-eye-movement (REM)-sleep behaviour disorder (RBD) is a prodromal condition of α-synucleinopathies, including Parkinson’s disease (PD) and dementia with Lewy bodies (DLB), and provides a unique opportunity to study and understand early-stage neurodegeneration. Previous research on RBD has identified common risk loci in several PD and DLB genes, but the full extent of its genetic architecture and disease-specific associations remain unknown. 

In the present study, we aimed to identify novel loci associated with RBD by conducting an updated European RBD genome-wide association study (GWAS). We performed a meta-analysis using a previous RBD GWAS and newly genotyped RBD patients and healthy controls, resulting in a total of 3,564 RBD cases and 140,393 controls analyzed. We further assessed the functional impact of new associations using fine-mapping and colocalization analyses with adult brain eQTLs. To confirm if novel RBD loci are implicated in additional synucleinopathy phenotypes, we examined common variant associations in GWASs of PD risk, DLB risk, PD with RBD risk, PD age at onset, and various measures of progression. Additionally,  rare variants were analyzed with SKAT-O in PD and DLB whole genome sequencing cohorts. Lastly, we analyzed known PD- and DLB-risk loci for associations with RBD. 

We identified RCOR1 as a novel and unique risk locus for RBD, with fine-mapping suggesting the association to be driven by a variant in intron 2. Colocalization revealed no evidence for shared genetic signals with brain eQTLs, suggesting the association is not mediated by altered expression in adult brain tissue. Common and rare variants in the RCOR1 locus were not associated with additional α-synucleinopathy phenotypes. We observed several PD and DLB risk loci to be associated with RBD, including MAPT, STK39, and SIPA1L2. The MAPT associated variant tags the protective H2 haplotype, consistent with previous findings in PD. Several SNCA variants previously associated with PD or DLB were also associated with RBD, but with differing directions of effect, highlighting the complex genetic landscape underlying α-synucleinopathy subtypes. This study identifies new RBD loci and strengthens our understanding of the genetic modifiers underlying RBD. Further replication and functional studies will be required to validated our findings and explore their biological implications. 

### GWAS_new_cohorts.sh
Contains code used to run basic GWASs in our new iRBD datasets. 

### meta_analysis.sh
Contains code used to run meta-analyses with the summary statistics from the previous RBD GWAS. 

### RCOR1_PD_DLB.sh
Contains code used to assess common and rare variant associations in the RCOR1 region with PD, DLB, PD with RBD, and more. 
