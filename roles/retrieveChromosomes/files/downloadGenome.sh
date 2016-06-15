#!/bin/bash
#FTPROOT="ftp://ftp.ensembl.org/pub/release-81/fasta/mus_musculus/dna"
#INPREFIX=Mus_musculus.GRCm38.dna
#NUMCHROMS=19
#OUTFOLDER="/home/ubuntu/genomes/MmGRCm38dna"
#OUTPREFIX="Mm_GRC38_dna"

##### get the options
while getopts R:I:N:F:P: o
do	case "$o" in
	R) FTPROOT="$OPTARG";;
	I) INPREFIX="$OPTARG";;
	N) NUMCHROMS="$OPTARG";;
	F) OUTFOLDER="$OPTARG";;
	P) OUTPREFIX="$OPTARG";;
	[?]) print >&2 "Usage: $0 [-R FTPROOT] [-I INPREFIX [-N NUMCHROMS] [-F OUTFOLDER] [-P OUTPREFIX]"
		exit 1;;
	esac
done

outfileOnlyChromosomes="${OUTPREFIX}_onlyChromosomes.fa"
outfileEverything="${OUTPREFIX}.fa"

cd "$OUTFOLDER"
for chromNumber in $(seq ${NUMCHROMS}) MT X Y
do chromFile="${INPREFIX}.chromosome.${chromNumber}.fa"
wget -O "${chromFile}.gz" "${FTPROOT}/${chromFile}.gz"
zcat "${chromFile}.gz" >> "$outfileOnlyChromosomes"
rm "${chromFile}.gz"
#gunzip "${chromFile}.gz"
#sed -i "0,/>.*/s//>chr${chromNumber}/" "$chromFile"
#cat "$chromFile" >> "$outfileOnlyChromosomes"
#rm "$chromFile"
done

# add the non-chromosomal regions 
chromFile="${INPREFIX}.nonchromosomal.fa"
wget -O "${chromFile}.gz" "${FTPROOT}/${chromFile}.gz"
cp "$outfileOnlyChromosomes" "$outfileEverything"
zcat "${chromFile}.gz" >> "$outfileEverything"
rm "${chromFile}.gz"
#gunzip "${chromFile}.gz"
#cat "$chromFile" >> "$outfileEverything"
#rm "$chromFile"
