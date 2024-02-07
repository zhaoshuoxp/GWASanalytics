#!/bin/bash

#Default settings
SHELL_FOLDER=$(dirname "$0")
# clean the GWAS file, replace all special characters to "_", including \/()[]<>+|%$&,;:'" and space. Otherwise will disrupt grep command
gwas=${SHELL_FOLDER}/GwasCatalog.genes
mod="group"
# CardiogramplusC4D is 52, so 50 as default to make sure including it
num=50
group=${SHELL_FOLDER}/group.tsv

help(){
	cat <<-EOF
	Usage: GeneGWASCatalog.sh <options> <input.genelist>

	### INPUT: Gene list in TEXT format, one gene per row ###
	This script will intersect the input Gene list with GWAS catalog Genes and calculate the enrichment for the overlaps.
	### R pacakges ggplot2/ggrepel/ggsci requried for plotiing ###

	Options:
	-a Running intersection with all the diseases in the catalog rather than grouped diseases only by default
	-n [int] Diseases included for at lease N Genes in the GWAS file, 50 as default. Only work with -a is ON
	-g [str] Custom GWAS file, three columns: chr,pos,disease, TAB delimited. Will use default GwasCatalog.gene under the script path if not designated.
	-p [str] Custom grouped disease file, two columns: keywords,diseases, TAB delimited. Only work with -a is OFF
	-h Print this help message
EOF
	exit 0
}

all_catalog(){
	# when comparing int/double in awk command, use -v to assign external vars to internal vars, otherwise the numbers will be treated as empty string
	disease=$(cut -f2 $gwas |sort |uniq -c|sort -k1,1nr|awk -v val=$num '$1+0>=val' |awk '{print $2}')
	echo "Grepping diseases"
	for i in ${disease[@]};do
		# didn't replace - to _, grep treats "-" as a range symbol, so have to add "\" before it
		grep -w -i "$(echo $i|sed 's/-/\\-/g')" $gwas |awk -v OFS="\t" '{print $1}' > ${i}.gwascatalog.gene &
	done
	wait 
	echo "Grepping Done"
	
	echo "Intersecting overlaps"
	for i in *.gwascatalog.gene;do
		cat $i $1|sort |uniq -d  > ${pre}_${i/.gwascatalog.gene/}.overlap &
	done
	wait
	echo "Intersecting done"
}

group_catalog(){
	echo "Grepping diseases"
	while read line; do
		keyword=$(echo $line |awk '{print $1}')
		disease=$(echo $line |awk '{print $2}')
		# make sure to remove ${disease}.gwascatalog.gen after running because of >>
		grep "$keyword" $gwas |awk -v OFS="\t" '{print $1}' >> ${disease}.gwascatalog.gene 
	done < $group
	for i in $(cut -f2 $group|sort -u);do
		sort -u ${i}.gwascatalog.gene > 1 && mv 1 ${i}.gwascatalog.gene
	done
	echo "Grepping Done"
	
	echo "Intersecting overlaps"
	for i in *.gwascatalog.gene;do
		cat $i $1 |sort|uniq -d  > ${pre}_${i/.gwascatalog.gene/}.overlap &
	done
	echo "Intersecting done"
}

enrich(){
	# enrich analysis, make sure use sort -u to remove duplicates
	fra=$(sort -u $1|awk 'END{print NR/61256}')
	echo -e "Disease\tOverlaps\tCatalogNumber\tFraction\tOdds">data.tsv
	total=$(cat *.gwascatalog.gene|sort -u |wc -l |awk '{print $1}')
	hits_total=$(cat *.overlap|sort -u |wc -l |awk '{print $1}')
	for i in *.gwascatalog.gene;do
		hits=$(cat ${pre}_${i/.gwascatalog.gene/}.overlap|sort -u|wc -l|awk '{print $1}')
		if [ $hits -gt 0 ];then
			cat_total=$(wc -l $i|awk '{print $1}')
			fold="$(awk 'BEGIN {print (("'"$hits"'"/"'"$hits_total"'")/("'"$cat_total"'"/"'"$total"'"))}')"
			echo -e "${i/.gwascatalog.gene/}\t$hits\t$cat_total\t$fra\t$fold" >> data.tsv
		fi
	done
}


# no ARGs error
if [ $# -lt 1 ];then
	help
	exit 1
fi

while getopts "hag:n:p:" arg
do
	case $arg in
		# only include grouped catalog if no -a
		a) mod="all";;
		# Custom GWAS file
		g) gwas=$OPTARG;;
		# Custom N number cutoff for catalog including
		n) num=$OPTARG;;
		# Custom grouped disease file
		p) group=$OPTARG;;
		h) help ;;
		?) help
			exit 1;;
	esac
done

shift $(($OPTIND - 1))

main(){
	if [ $mod != 'all' ];then 
		echo "Running on grouped catalog by default"
		# don't forget $1 to run
		# take basename of $1 to keep output, otherwise errors when using input with path 
		pre=$(basename $1)
		group_catalog $1
		enrich $1
		
		cat >plot.r<<-EOF
		#!/usr/bin/env Rscript
		library("ggplot2")
		library("ggrepel")
		library("ggsci")

		data <- read.table('data.tsv', header=T)
		data\$Pval<-dbinom(data\$Overlaps,data\$CatalogNumber,data\$Fraction,log=F)
		data\$Padj<-p.adjust(data\$Pval)
		data[order(data[,6]),]->data
		data\$rank<-nrow(data):1

		png(file='output.png',height = 7, width = 8, res=600, units = "in", family="Arial")
		ggplot(data, aes(x=Odds, y=-log10(Padj),label=Disease)) + 
			geom_point(shape=19, alpha=0.5, aes(size=Overlaps,color=Disease)) + 
			xlab("Odds Ratio") + ylab("-logPval") + 
			ggtitle ("GWAS SNPs enrichment - binomial test") + 
			scale_size(range = c(5, 20))+
			theme(axis.text=element_text(size=18), 
				axis.title=element_text(size=18))+
			guides(color = guide_legend(override.aes = list(size = 5)))+
			scale_color_npg()+
			geom_text_repel(aes(label=Disease, color = Disease), size=6,show.legend=F)
		dev.off()

		write.table(data,'data.tsv',sep='\t',quote=F,row.names=F, col.names=T)
	
		EOF
		chmod 755 plot.r 
		./plot.r
	else
		echo "Running on all diseases"
		pre=$(basename $1)
		all_catalog $1
		enrich $1
	fi
	rm *.gwascatalog.gene ${pre}_*.overlap 
}

main $1

# check running status
if [ $? -ne 0 ]; then
	help
	exit 1
else
	echo "Run succeed"
fi

################ END ################
#          Created by Aone          #
#        quanyiz@stanford.edu       #
################ END ################