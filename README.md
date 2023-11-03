# GWAS Catalog analytics

-----
This is a shell script to find the intersected GWAS SNPs or genes within the input BED file or gene list. There are two modes for this script. The default is to read the designated group.tsv to search the keywords of diseases to obtain the SNPs or genes. While -a mode is to traversal all diseases in AllDiseases.txt. 

 * Bed2GWASCatalog.sh|GeneGWASCatalog.sh: These scripts will intersect the input BED file in hg38 or the gene list with GWAS catalog SNPs or genes, and calculate the enrichment using *binomial test* for the overlaps.

 * GwasCatalog.bed|GwasCatalog.genes : Source file of all GWAS associtations, SNP coordinates and 1bp plus, or gene offical names, TAB delimited. Converted from [GWAS](https://www.ebi.ac.uk/gwas/), version [1.0.2](https://www.ebi.ac.uk/gwas/api/search/downloads/alternative) 
 * group.tsv: Grouped diseases and keywords to search in default mode, TAB delimited.
 * AllDiseases.txt: All diseases in GWAS catalog.

> Requirements:
awk, sed, bedtools, R, R packages: ggplot2,ggrepel,ggsci

[![996.icu](https://img.shields.io/badge/link-996.icu-red.svg)](https://996.icu) [![LICENSE](https://img.shields.io/badge/license-Anti%20996-blue.svg)](https://github.com/996icu/996.ICU/blob/master/LICENSE)

----

#### Input

BED file, chr, pos1, pos2, TAB delimited, or gene list in TEXT file, one per row.

#### Usage

help message can be shown by `Bed2GWASCatalog.sh -h`

```shell
./Bed2GWASCatalog.sh
Usage: Bed2GWASCatalog.sh <options> <input.bed>

### INPUT: BED file ###
This script will intersect the input BED file in hg38 with GWAS catalog SNPs and calculate the enrichment for the overlaps.
### bedtools required for running, R pacakges ggplot2/ggrepel/ggsci requried for plotiing ###

Options:
-a Running intersection with all the diseases in the catalog rather than grouped diseases only by default
-n [int] Diseases included for at lease N SNPs in the GWAS file, 50 as default. Only work with -a is ON
-g [str] Custom GWAS file, three columns: chr,pos,disease, TAB delimited. Will use default GwasCatalog.bed under the script path if not designated.
-p [str] Custom grouped disease file, two columns: keywords,diseases, TAB delimited. Only work with -a is OFF
-h Print this help message
```

For default mode, edit `group.tsv` or provde a new one to the script by `-g`:

```shell
./Bed2GWASCatalog.sh TCF21_peaks_hg38.bed
#Custom group.tsv
./Bed2GWASCatalog.sh -p group.tsv TCF21_peaks_hg38.bed
```

For default mode, add `-a`:

```shell
./Bed2GWASCatalog.sh -a TCF21_peaks_hg38.bed
```

Custom GWAS file is also supported:

```shell
./Bed2GWASCatalog.sh -g your_gwas.bed TCF21_peaks_hg38.bed
```

> Usage of GeneGWASCatalog.sh is identical to Bed2GWASCatalog.sh



#### Output

- data.tsv: The text output of the enrichment analysis, one disease per line, TAB delimited.
- output.png: data.tsv ploting with ggplot2, only generated in default mode.

![output.png](https://raw.githubusercontent.com/zhaoshuoxp/GWASanalytics/main/output.png)

------

Author [@zhaoshuoxp](https://github.com/zhaoshuoxp)  
Nov 3 2023  

