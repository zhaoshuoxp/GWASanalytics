#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status.

# Default settings
SHELL_FOLDER=$(dirname "$0")
gwas="${SHELL_FOLDER}/GwasCatalog.bed"
mod="group"
num=50
group="${SHELL_FOLDER}/group.tsv"

help(){
	cat <<-EOF
	Usage: Bed2GWASCatalog.sh <options> <input.bed>

	### INPUT: BED file ###
	This script will intersect the input BED file in hg38 with GWAS catalog SNPs and calculate the enrichment for the overlaps.
	### bedtools required for running, R packages ggplot2/ggrepel/ggsci required for plotting ###

	Options:
	-a Running intersection with all the diseases in the catalog rather than grouped diseases only by default
	-n [int] Diseases included for at least N SNPs in the GWAS file, 50 as default. Only work with -a is ON
	-g [str] Custom GWAS file, three columns: chr,pos,disease, TAB delimited. Will use default GwasCatalog.bed under the script path if not designated.
	-p [str] Custom grouped disease file, two columns: keywords,diseases, TAB delimited. Only work with -a is OFF
	-h Print this help message
EOF
	exit 0
}

all_catalog(){
	disease=$(cut -f4 "$gwas" | sort | uniq -c | sort -k1,1nr | awk -v val="$num" '$1+0>=val' | awk '{print $2}')
	echo "Grepping diseases..."
	for i in ${disease[@]}; do
		grep -w -i "$(echo "$i" | sed 's/-/\\-/g')" "$gwas" | awk -v OFS="\t" '{print $1,$2,$2+1}' > "${i}.gwascatalog.bed" &
	done
	wait
	echo "Grepping done."
	
	echo "Intersecting overlaps..."
	for i in *.gwascatalog.bed; do
		bedtools intersect -a "$i" -b "$1" | sort -u > "${pre}_${i/.gwascatalog.bed/}.overlap" &
	done
	wait
	echo "Intersecting done."
}

group_catalog(){
	echo "Grepping diseases..."
	while read -r line; do
		keyword=$(echo "$line" | awk '{print $1}')
		disease=$(echo "$line" | awk '{print $2}')
		grep "$keyword" "$gwas" | awk -v OFS="\t" '{print $1,$2,$2+1}' >> "${disease}.gwascatalog.bed" 
	done < "$group"
	
	for i in $(cut -f2 "$group" | sort -u); do
		sort -u "${i}.gwascatalog.bed" > tmp_sort && mv tmp_sort "${i}.gwascatalog.bed"
	done
	echo "Grepping done."
	
	echo "Intersecting overlaps..."
	for i in *.gwascatalog.bed; do
		bedtools intersect -a "$i" -b "$1" | sort -u > "${pre}_${i/.gwascatalog.bed/}.overlap" &
	done
	wait
	echo "Intersecting done."
}

enrich(){
	# Cache chrom sizes to prevent redundant downloads
	if [ ! -f "hg38.chrom.sizes" ]; then
		echo "Downloading hg38 chrom sizes..."
		curl -sS https://hgdownload.cse.ucsc.edu/goldenpath/hg38/bigZips/hg38.chrom.sizes > hg38.chrom.sizes
	fi
	hg38=$(awk '{ sum+=($2)} END {print sum}' hg38.chrom.sizes)
	
	# Calculate coverage
	sort -k1,1V -k2,2n "$1" > tmp_sorted.bed
	bedtools merge -i tmp_sorted.bed -c 1 -o count > tmp_input.bed
	rm tmp_sorted.bed
	
	cov=$(awk '{ sum+=($3-$2)} END {print sum}' tmp_input.bed)
	echo "Coverage of BED file: $cov"
	fra=$(awk -v hg38="$hg38" '{ sum+=($3-$2)} END {print sum/hg38}' tmp_input.bed)
	echo "Input fraction of hg38: $fra"
	
	echo -e "Disease\tOverlaps\tCatalogNumber\tFraction\tFoldEnrichment" > data.tsv
	
	for i in *.gwascatalog.bed; do
		hits=$(cat "${pre}_${i/.gwascatalog.bed/}.overlap" | sort -u | wc -l | awk '{print $1}')
		if [ "$hits" -gt 0 ]; then
			cat_total=$(wc -l < "$i" | awk '{print $1}')
			fold=$(awk -v h="$hits" -v c="$cat_total" -v f="$fra" 'BEGIN {print (h / (c * f))}')
			echo -e "${i/.gwascatalog.bed/}\t$hits\t$cat_total\t$fra\t$fold" >> data.tsv
		fi
	done
}

# No ARGs error
if [ $# -lt 1 ]; then
	help
	exit 1
fi

while getopts "hag:n:p:" arg; do
	case $arg in
		a) mod="all";;
		g) gwas="$OPTARG";;
		n) num="$OPTARG";;
		p) group="$OPTARG";;
		h) help ;;
		?) help; exit 1;;
	esac
done

shift $((OPTIND - 1))

main(){
	file=$(basename "$1")
	pre=${file%.*}
	
	if [ "$mod" != 'all' ]; then 
		echo "Running on grouped catalog by default..."
		group_catalog "$1"
		enrich "$1"
		
		cat > plot.R <<-EOF
		#!/usr/bin/env Rscript
		library("ggplot2")
		library("ggrepel")
		library("ggsci")

		data <- read.table('data.tsv', header=TRUE)
		
		data\$Pval <- pbinom(data\$Overlaps - 1, data\$CatalogNumber, data\$Fraction, lower.tail=FALSE)
		data\$Padj <- p.adjust(data\$Pval, method = "BH")
		
		data <- data[order(data\$Pval, decreasing=TRUE), ]
		data\$rank <- nrow(data):1

		# Render plot with transparent background
		png(file='output.png', height = 7, width = 8, res=600, units = "in", family="sans", bg="transparent")
		
		ggplot(data, aes(x=FoldEnrichment, y=-log10(Padj), label=Disease)) + 
			geom_point(shape=19, alpha=0.6, aes(size=Overlaps, color=Disease)) + 
			xlab("Fold Enrichment (Observed / Expected)") + ylab("-log10(Padj)") + 
			ggtitle("GWAS SNPs Enrichment (Binomial Test)") + 
			scale_size(range = c(5, 20)) +
			scale_color_npg() +
			theme_classic(base_size = 14) +
			theme(
				axis.text = element_text(size=14), 
				axis.title = element_text(size=16, face="bold"),
				plot.title = element_text(size=16, face="bold", hjust=0.5),
				panel.background = element_rect(fill = "transparent", colour = NA),
				plot.background = element_rect(fill = "transparent", colour = NA),
				panel.border = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = NA)
			) +
			guides(color = guide_legend(override.aes = list(size = 5))) +
			geom_text_repel(aes(label=Disease, color = Disease), size=5, show.legend=FALSE, box.padding = 0.5)
			
		invisible(dev.off())

		write.table(data, 'data.tsv', sep='\t', quote=FALSE, row.names=FALSE, col.names=TRUE)
		EOF
		
		chmod 755 plot.R 
		./plot.R
		mv output.png "${pre}.png"
	else
		echo "Running on all diseases..."
		all_catalog "$1"
		enrich "$1"
	fi
	
	# Cleanup temp files
	rm -f *.gwascatalog.bed "${pre}_"*.overlap tmp_input.bed plot.R
	mv data.tsv "${pre}.tsv"
}

main "$1"

echo "Run succeeded."