# AAD_t1d

Code related to analysis examining whether any variants associated with type 1 diabetes (T1D) have  differential effect sizes between those diagnosed early (<7 years) or late (>13 years).

# Data

The underlying data is sensitive and thus not publically available, though results will be posted on here and paths to underlieing data on the DIL servers will be placed here. 
This repository has been created in an attempt to make it clear exactly how analyses and data processing were carried out.

# Prerequisites

Many R packages and other software are used throughout this analysis. 
R packages:

```
snpStats: https://bioconductor.org/packages/release/bioc/html/snpStats.html
annotSnpStats: https://github.com/chr1swallace/annotSnpStats
tidyverse: https://github.com/tidyverse/
gridExtra: https://cran.r-project.org/web/packages/gridExtra/index.html
multinomRob: https://cran.r-project.org/web/packages/multinomRob/index.html 
GUESSFM: https://github.com/chr1swallace/GUESSFM
R2GUESS: https://cran.r-project.org/web/packages/R2GUESS/index.html
HiBag: https://bioconductor.org/packages/release/bioc/html/HIBAG.html
```

Other software:
```
PLINK (version 1.9 and 2.0 both used): https://www.cog-genomics.org/plink/2.0/
KING version 2.1.6: http://people.virginia.edu/~wc9c/KING/manual.html
GUESS: http://www.bgx.org.uk/software/guess.html
IMPUTE2: http://mathgen.stats.ox.ac.uk/impute/impute_v2.html 
SNPTEST: https://mathgen.stats.ox.ac.uk/genetics_software/snptest/snptest.html
QCTOOL: https://www.well.ox.ac.uk/~gav/qctool_v2/
```

# Driver scripts and running order

Now outlining the analysis pipeline so if data were made available to a collaborator, they should be able to exactly reproduce results.


### HLA analysis:
```
readin_ichip.R	-> Reads in immunoChip data from all available collections, removes related individuals, inferred from SNP data (to 2nd degree), merges in other phenotype data including principal components generated from SNP data
hla_imputation_redo.R	-> Uses HIBAG to impute HLA classical alleles for DRB1, DQA1, DAB1, A, B and C
hla_check.R	-> Compares imputed classical HLA alleles to directly genotyped on a subset of individuals with classical genotyping data available
hla_readin.R	-> Reads in the imputed classical HLA genotypes and converts to class II haplotypes that are associated with T1D and codes up correctly for analysis. Removes alleles called with <0.5 posterior probability.
hla_multinomial_n.R	-> Performs multinomial logistic regressions, from which heterogeneity tests are performed comparing constrained and unconstrained models with regards to the effect size between <7s and >13s.
```

### non-HLA analysis (primary):
```
impute_missing_snps.R	-> Where SNPs of interest were removed due to QC, imputing them so we can use them in the analyses. Using IMPUTE2.
imputation_readin.R	-> Reads in the imputation results, keeping only the index SNP of interest once checking it passes all QC checks
multinomial_redo_n.R	-> Carries out multinomial regression on all 55 T1D-associted variants to test for heterogneity of effect size between <7s and <13s using all individuals in analysis.
multinomial_readin.R	-> Reads in the results form the above script, which were saved for examiniation in this script
```
### non-HLA analyses (sensitivity):
```
multinomial_split_iter.R	-> samples 50% of population 100 times and performs multinomial regression on all 200 datasets to test stability og results in non-HLA regions
multinomial_split_iter_readin.R	-> reads in results of above analysis and produces Figure to show stability across iterations.
multinomial_uk_only.R	-> multinomial regression of the 55 regions using only individuals from the UK/northern Ireland to ensure associations aren't due to population structure
multinomial_readin_uk_only.R	-> Reads in the	results	form the above script, which were saved	for examiniation in this script
multinomial_sensitivity_to_cutoff.R -> multinomial regression of the 55 regions using all individuals but cutting strata at <6 and <5 instead of <7 and re-examining results
```

### fine-mapping:
```
guessfm_setup_uk_only.R	-> sets up GUESSFM procudure by imptuting 0.5Mb around the regions of interest, removing low quality imputed variants then running stochastic search. UK/NI individuals only included in analysis.
imputed_test_uk.R	-> gets summary stats for the credible SNPs, these are then used in the next script.
guessfm_readin_uk_only.R	-> reads in the results of the stochastic search and does post processing including calculating ABFs for each model including those that were pruned out for stoachstic search
```

### other:
```
table_1.R	-> generates Table 1 in manuscript, baseline characteristics
```

# Data locations of DIL server

If anyone needs to recreate this analysis, hopefully this will be useful to anyone working at the DIL with regards to where the raw data underlying this analysis is saved.

### Immunochip data
UK case-controls: 
```
/fs3/projects/todd/ipswich/Immunochip/PLINK/distribution/t1d-cases-preqc.(bed,bim,fam)
/fs3/projects/todd/ipswich/Immunochip/PLINK/distribution/uva-controls-preqc.(bed,bim,fam)
/fs3/projects/todd/ipswich/Immunochip/PLINK/distribution/sanger-controls-preqc.(bed,bim,fam)
```

T1DGC families:
```
/fs3/projects/todd/ipswich/Immunochip/PLINK/distribution/t1dgc-asp-preqc.(bed,bim,fam)
```

Warren, IDDMGEN, T1DGEN, NI GRID collections:
```
/fs3/projects/todd/ipswich/Immunochip/PLINK/distribution/finn-preqc.(bed,bim,fam)
```

### Age-at-diagnosis data:
UK cases (GRID):
```
/fs3/projects/todd/ipswich/Immunochip/support/casecontrol/t1d-subject-lookup-2015-08-06.csv
```

T1DGC families:
```
/fs3/projects/todd/chiswick/store/t1dgc-datasets/T1DGC.2011.03_Golden_Pedigree/T1DGC.2011.03_Golden_Pedigree/Resources/T1DGC.2011.03_Resources.csv
```

Warren, IDDMGEN, T1DGEN, NI GRID:
```
/fs3/projects/todd/ipswich/shared/researcher-data-requests/jamie-2016-01-21/
```

### Direct HLA typing data:
```
/fs3/projects/todd/chiswick/store/hla/summary-data/cc-hla-2015-11-13.tab
```


