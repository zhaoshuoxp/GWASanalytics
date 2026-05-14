#!/bin/bash
set -e # Exit immediately if a command exits with a non-zero status

# Default settings
SHELL_FOLDER=$(dirname "$0")
gwas="${SHELL_FOLDER}/GwasCatalog.genes"
mod="group"
num=50
group="${SHELL_FOLDER}/group.tsv"

help(){
	cat <<-EOF
	Usage: GeneGWASCatalog.sh <options> <input.genelist>

	### INPUT: Gene list in TEXT format, one gene per row ###
	This script will intersect the input Gene list with GWAS catalog Genes and calculate the enrichment for the overlaps.
	### R packages ggplot2/ggrepel/ggsci required for plotting ###

	Options:
	-a Running intersection with all the diseases in the catalog rather than grouped diseases only by default
	-n [int] Diseases included for at least N Genes in the GWAS file, 50 as default. Only work with -a is ON
	-g [str] Custom GWAS file, three columns: chr,pos,disease, TAB delimited. Will use default GwasCatalog.genes under the script path if not designated.
	-p [str] Custom grouped disease file, two columns: keywords,diseases, TAB delimited. Only work with -a is OFF
	-h Print this help message
EOF
	exit 0
}

all_catalog(){
	disease=$(cut -f2 "$gwas" | sort | uniq -c | sort -k1,1nr | awk -v val="$num" '$1+0>=val' | awk '{print $2}')
	echo "Grepping diseases..."
	for i in ${disease[@]}; do
		grep -w -i "$(echo "$i" | sed 's/-/\\-/g')" "$gwas" | awk -v OFS="\t" '{print $1}' > "${i}.gwascatalog.gene" &
	done
	wait 
	echo "Grepping done."
	
	echo "Intersecting overlaps..."
	for i in *.gwascatalog.gene; do
		# Robust intersection using comm
		comm -12 <(sort -u "$i") <(sort -u "$1") > "${pre}_${i/.gwascatalog.gene/}.overlap" &
	done
	wait
	echo "Intersecting done."
}

group_catalog(){
	echo "Grepping diseases..."
	while read -r line; do
		keyword=$(echo "$line" | awk '{print $1}')
		disease=$(echo "$line" | awk '{print $2}')
		grep "$keyword" "$gwas" | awk -v OFS="\t" '{print $1}' >> "${disease}.gwascatalog.gene" 
	done < "$group"
	
	for i in $(cut -f2 "$group" | sort -u); do
		sort -u "${i}.gwascatalog.gene" > tmp_sort && mv tmp_sort "${i}.gwascatalog.gene"
	done
	echo "Grepping done."
	
	echo "Intersecting overlaps..."
	for i in *.gwascatalog.gene; do
		comm -12 <(sort -u "$i") <(sort -u "$1") > "${pre}_${i/.gwascatalog.gene/}.overlap" &
	done
	wait # Added missing wait to ensure background intersections finish
	echo "Intersecting done."
}

enrich(){
	# Assuming 61256 is the total number of annotated human genes (e.g., Gencode)
	fra=$(sort -u "$1" | wc -l | awk '{print $1/61256}')
	echo -e "Disease\tOverlaps\tCatalogNumber\tFraction\tFoldEnrichment" > data.tsv
	
	total=$(cat *.gwascatalog.gene | sort -u | wc -l | awk '{print $1}')
	hits_total=$(cat *.overlap | sort -u | wc -l | awk '{print $1}')
	
	for i in *.gwascatalog.gene; do
		hits=$(cat "${pre}_${i/.gwascatalog.gene/}.overlap" | sort -u | wc -l | awk '{print $1}')
		if [ "$hits" -gt 0 ]; then
			cat_total=$(wc -l < "$i" | awk '{print $1}')
			fold=$(awk -v h="$hits" -v ht="$hits_total" -v c="$cat_total" -v t="$total" 'BEGIN {print ((h/ht)/(c/t))}')
			echo -e "${i/.gwascatalog.gene/}\t$hits\t$cat_total\t$fra\t$fold" >> data.tsv
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
	file=$(basename "$1") # FIXED: properly defined the file variable
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
		# Cumulative right-tail binomial probability
		data\$Pval <- pbinom(data\$Overlaps - 1, data\$CatalogNumber, data\$Fraction, lower.tail=FALSE)
		data\$Padj <- p.adjust(data\$Pval, method = "BH")
		data <- data[order(data[,6]), ]
		data\$rank <- nrow(data):1

		# Render plot with transparent background
		png(file='output.png', height = 7, width = 8, res=600, units = "in", family="Arial", bg="transparent")
		
		ggplot(data, aes(x=FoldEnrichment, y=-log10(Padj), label=Disease)) + 
			geom_point(shape=19, alpha=0.5, aes(size=Overlaps, color=Disease)) + 
			xlab("Fold Enrichment") + ylab("-log10(Padj)") + 
			ggtitle("GWAS Gene Enrichment - Binomial Test") + 
			scale_size(range = c(5, 20)) +
			scale_color_npg() +
			theme_classic() +
			theme(
				axis.text = element_text(size=18), 
				axis.title = element_text(size=18),
				panel.background = element_rect(fill = "transparent", colour = NA),
				plot.background = element_rect(fill = "transparent", colour = NA),
				panel.border = element_blank(),
				legend.background = element_rect(fill = "transparent", colour = NA)
			) +
			guides(color = guide_legend(override.aes = list(size = 5))) +
			geom_text_repel(aes(label=Disease, color = Disease), size=6, show.legend=FALSE)
			
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
	
	# Cleanup
	rm -f *.gwascatalog.gene "${pre}_"*.overlap plot.R
	mv data.tsv "${pre}.tsv"
}

main "$1"

echo "Run succeeded."