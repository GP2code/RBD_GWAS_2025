# European RBD genome-wide association study (GWAS)

`GP2 ❤️ Open Science 😍`

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.15794660.svg)](https://doi.org/10.5281/zenodo.15794660)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**Last Updated:** June 2025

## Summary
Rapid-eye-movement (REM)-sleep behaviour disorder (RBD) is a prodromal condition of α-synucleinopathies, including Parkinson’s disease (PD) and dementia with Lewy bodies (DLB), and provides a unique opportunity to study and understand early-stage neurodegeneration. Previous research on RBD has identified common risk loci in several PD and DLB genes, but the full extent of its genetic architecture and disease-specific associations remain unknown. 

In the present study, we aimed to identify novel loci associated with RBD by conducting an updated European RBD genome-wide association study (GWAS). We performed a meta-analysis using a previous RBD GWAS and newly genotyped RBD patients and healthy controls, resulting in a total of 3,564 RBD cases and 140,393 controls analyzed. We further assessed the functional impact of new associations using fine-mapping and colocalization analyses with adult brain eQTLs. To confirm if novel RBD loci are implicated in additional synucleinopathy phenotypes, we examined common variant associations in GWASs of PD risk, DLB risk, PD with RBD risk, PD age at onset, and various measures of progression. Additionally,  rare variants were analyzed with SKAT-O in PD and DLB whole genome sequencing cohorts. Lastly, we analyzed known PD- and DLB-risk loci for associations with RBD. 

We identified RCOR1 as a novel and unique risk locus for RBD, with fine-mapping suggesting the association to be driven by a variant in intron 2. Colocalization revealed no evidence for shared genetic signals with brain eQTLs, suggesting the association is not mediated by altered expression in adult brain tissue. Common and rare variants in the RCOR1 locus were not associated with additional α-synucleinopathy phenotypes. We observed several PD and DLB risk loci to be associated with RBD, including MAPT, STK39, and SIPA1L2. The MAPT associated variant tags the protective H2 haplotype, consistent with previous findings in PD. Several SNCA variants previously associated with PD or DLB were also associated with RBD, but with differing directions of effect, highlighting the complex genetic landscape underlying α-synucleinopathy subtypes. This study identifies new RBD loci and strengthens our understanding of the genetic modifiers underlying RBD. Further replication and functional studies will be required to validated our findings and explore their biological implications. 


---

## **Citation**

Preprint link: *pending*


---
## Helpful Links 
- [GP2 Website](https://gp2.org/)
    - [GP2 Cohort Dashboard](https://gp2.org/cohort-dashboard-advanced/)
- [Introduction to GP2](https://movementdisorders.onlinelibrary.wiley.com/doi/10.1002/mds.28494)
    - [Other GP2 Manuscripts (PubMed)](https://pubmed.ncbi.nlm.nih.gov/?term=%22global+parkinson%27s+genetics+program%22)

## Data Statement

* Data for iRBD samples are in the process of being added to GP2, and will be available in future releases.
* 23andMe summary statistic data can be obtained by qualified researchers from 23andMe upon completion of a data access request (https://research.23andme.com/dataset-access/). 
* The brain eQTL datasets used for colocalization analysis were obtained from the GTEx Portal on 04/28/2025 and can be downloaded from https://gtexportal.org/home/.


## Figures and Supplementary Figures
*(pending publication)*

## Tables and Supplementary Tables
*(pending publication)*

---

## Repository Orientation
- This is the online repository for **"European RBD genome-wide association study (GWAS)"** study, and contains code for all analyses excluding fine-mapping code which can be found at https://github.com/daria-nikanorova/Fine-mapping_SCARB2_CTSB.

```
THIS_REPO
├── README.md
├── LICENSE
├── 00_GWAS_new_cohorts.sh
├── 01_meta_analysis.sh
├── 02_RCOR1_PD_DLB.sh
└── 03_colocalization.R
```

---

## Script Details

* Languages: Bash, R


| Script | Description |
| ------ | ----------- |
| 00_GWAS_new_cohorts.sh | Perform basic GWAS on new iRBD cohorts| | 
| 01_meta_analysis.sh | Run meta-analysis on summary statistics from new cohorts and 23andMe cohort |
| 02_RCOR1_PD_DLB.sh | Run common and rare variant association investigations for PD, DLB, PD with RBD, and more | 
| 03_colocalization.R | Run colocalization analyses with GTEx eQTL data |


---
## Software 


| Software | Version | Resource URL | RRID | Notes |
|----------|---------|--------------|------|-------|
| PLINK | 1.9 | https://www.cog-genomics.org/plink/ | RRID:SCR_001757 | Used for QC and various genetic analyses |
| R | 4.2 | http://www.r-project.org/ | RRID:SCR_001905 | Used for data processing, colocalization, and fine-mapping. Packages: data.table, qqman, dplyr, SKAT, coloc, susieR, magrittr, ggplot |
| perl | 5.36.1 | https://www.perl.org/get.html | RRID:SCR_018313 | Used with ANNOVAR for variant annotation |
| metal | 2011-03-25 | https://csg.sph.umich.edu/abecasis/metal/download/ | RRID:SCR_002013 | Used for performing meta-analyses of GWAS summary statistics |
| FINEMAP | 1.4.2 | http://www.christianbenner.com | NA | Used for fine-mapping |
| gcta64 | 1.94.1 | https://yanglab.westlake.edu.cn/software/gcta/ | NA | Used for conditional and joint analyses |
| ANNOVAR | NA | https://annovar.openbioinformatics.org/en/latest/ | RRID:SCR_012821 | Used with perl for variant annotation |

