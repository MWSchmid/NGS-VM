#!/bin/bash
me=$(basename "$0")
OUTDIR=$(dirname "$0")
declare -a NAMES
declare -a LINKS
declare -a FILES
# samtools 1.2
NAMES[1]="samtools"
LINKS[1]="https://github.com/samtools/samtools/releases/download/1.2/samtools-1.2.tar.bz2"
FILES[1]="samtools-1.2.tar.bz2"
# bamtools 2.4.0
NAMES[2]="bamtools"
LINKS[2]="https://github.com/pezmaster31/bamtools/archive/master.zip"
FILES[2]="bamtools-2.4.0.zip"
# bedtools 2.24.0
NAMES[3]="bedtools"
LINKS[3]="http://github.com/arq5x/bedtools2/archive/master.zip"
FILES[3]="bedtools-2.24.0.zip"
# bowtie 1.1.2
NAMES[4]="bowtie1"
LINKS[4]="http://sourceforge.net/projects/bowtie-bio/files/bowtie/1.1.2/bowtie-1.1.2-src.zip"
FILES[4]="bowtie-1.1.2.zip"
# bowtie2 2.2.5
NAMES[5]="bowtie2"
LINKS[5]="http://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.2.5/bowtie2-2.2.5-source.zip"
FILES[5]="bowtie-2.2.5.zip"
# bwa 0.7.12-r1039
NAMES[6]="bwa"
LINKS[6]="https://github.com/lh3/bwa/archive/master.zip"
FILES[6]="bwa-0.7.12.zip"
# soap 2.20
NAMES[7]="soapAligner"
LINKS[7]="http://soap.genomics.org.cn/down/SOAPaligner-v2.20-src.tar.gz"
FILES[7]="soapAligner-2.20.tar.gz"
# soap builder 2.20
NAMES[8]="soapBuilder"
LINKS[8]="http://soap.genomics.org.cn/down/SOAPaligner-v2.20-src_builder.tar.gz"
FILES[8]="soapBuilder-2.20.tar.gz"
# star 2.4.2a
NAMES[9]="star"
LINKS[9]="https://github.com/alexdobin/STAR/archive/STAR_2.4.2a.tar.gz"
FILES[9]="star-2.4.2a.tar.gz"
# subread 1.4.6-p4
NAMES[10]="subread"
LINKS[10]="http://sourceforge.net/projects/subread/files/subread-1.4.6-p4/subread-1.4.6-p4-source.tar.gz"
FILES[10]="subread-1.4.6-p4.tar.gz"
# HTSeq 0.6.1
NAMES[11]="HTSeq"
LINKS[11]="https://pypi.python.org/packages/source/H/HTSeq/HTSeq-0.6.1.tar.gz"
FILES[11]="HTSeq-0.6.1.tar.gz"
# RSEM 1.2.22
NAMES[12]="RSEM"
LINKS[12]="http://deweylab.biostat.wisc.edu/rsem/src/rsem-1.2.22.tar.gz"
FILES[12]="RSEM-1.2.22.tar.gz"
# Rcount 1.0
NAMES[13]="Rcount"
LINKS[13]="https://github.com/MWSchmid/Rcount/blob/master/other/linux_64bit.zip?raw=true"
FILES[13]="Rcount-1.0.zip"
# SPP 1.14
NAMES[14]="SPP"
LINKS[14]="https://github.com/hms-dbmi/spp/archive/1.14.tar.gz"
FILES[14]="spp-1.14.tar.gz"
# phantompeakqualtools 1.1
NAMES[15]="phantompeakqualtools"
LINKS[15]="https://phantompeakqualtools.googlecode.com/files/ccQualityControl.v.1.1.tar.gz"
FILES[15]="phantompeakqualtools-1.1.tar.gz"
# MAnorm 1.0
NAMES[16]="MAnorm"
LINKS[16]="http://bcb.dfci.harvard.edu/~gcyuan/MAnorm/MAnorm_data/MAnorm_Linux_R_Package.zip"
FILES[16]="MAnorm-1.0.zip"
# MACS 2.1.0
NAMES[17]="MACS"
LINKS[17]="https://github.com/taoliu/MACS/archive/master.zip"
FILES[17]="MACS-2.1.0.zip"
# CEAS 1.0.2
NAMES[18]="CEAS"
LINKS[18]="http://liulab.dfci.harvard.edu/CEAS/src/CEAS-Package-1.0.2.tar.gz"
FILES[18]="CEAS-1.0.2.tar.gz"
# GEM 2.6
NAMES[19]="GEM"
LINKS[19]="http://groups.csail.mit.edu/cgs/gem/download/gem.v2.6.tar.gz"
FILES[19]="GEM-2.6.tar.gz"
# multovl 1.3
NAMES[20]="multovl"
LINKS[20]="http://sourceforge.net/projects/multovl/files/multovl-1.3.8203e7d-Source.tar.gz"
FILES[20]="multovl-1.3.tar.gz"
# SICER 1.1
NAMES[21]="SICER"
LINKS[21]="http://home.gwu.edu/~wpeng/SICER_V1.1.tgz"
FILES[21]="SICER-1.1.tar.gz"
# PeakAnalyzer 1.4
NAMES[22]="PeakAnalyzer"
LINKS[22]="http://www.ebi.ac.uk/sites/ebi.ac.uk/files/groups/bertone/software/PeakAnalyzer_1.4.tar.gz"
FILES[22]="PeakAnalyzer-1.4.tar.gz"
# jChIP 1.0
NAMES[23]="jChIP"
LINKS[23]="http://sourceforge.net/projects/jchip/files/jChIP.jar"
FILES[23]="jChIP-1.0.jar"
# PAPST 1.0
NAMES[24]="PAPST"
LINKS[24]="https://github.com/paulbible/papst/archive/master.zip"
FILES[24]="PAPST-1.0.zip"
# fastqc 0.11.4
NAMES[25]="fastqc"
LINKS[25]="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/fastqc_v0.11.4.zip"
FILES[25]="fastqc-0.11.4.zip"
# trimGalore 0.4.0
NAMES[26]="trimGalore"
LINKS[26]="http://www.bioinformatics.babraham.ac.uk/projects/trim_galore/trim_galore_v0.4.0.zip"
FILES[26]="trimGalore-0.4.0.zip"
# bismark 0.14.3
NAMES[27]="bismark"
LINKS[27]="http://www.bioinformatics.babraham.ac.uk/projects/bismark/bismark_v0.14.3.tar.gz"
FILES[27]="bismark-0.14.3.tar.gz"
# picard 1.140 
NAMES[28]="picard"
LINKS[28]="https://github.com/broadinstitute/picard/archive/master.zip"
FILES[28]="picard-1.140.tar.gz"
# htsjdk 
NAMES[29]="htsjdk"
LINKS[29]="https://github.com/samtools/htsjdk/archive/master.zip"
FILES[29]="htsjdk-ukn.tar.gz"
# trimmomatic 0.33
NAMES[30]="trimmomatic"
LINKS[30]="http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.33.zip"
FILES[30]="trimmomatic-0.33.zip"
# fqtrim 0.9.4
NAMES[31]="fqtrim"
LINKS[31]="http://ccb.jhu.edu/software/fqtrim/dl/fqtrim-0.9.4.tar.gz"
FILES[31]="fqtrim-0.9.4.tar.gz"
#something
#NAMES[]=""
#LINKS[]=""
#FILES[]=""

for NUM in $(seq 31)
do echo "   - ${NAMES[$NUM]}"
done

# ask for confirmation
read -p "Do you want to download and replace all (y/n)?" -n 1 -r
echo
if [[ ! $REPLY =~ ^[Yy]$ ]]
    then
    exit 1
fi

# download everything
mkdir -p $OUTDIR
for NUM in $(seq 31)
do echo "downloading ${NAMES[$NUM]}"
curl -L -o $OUTDIR/${FILES[$NUM]} ${LINKS[$NUM]}
done

# repack to get rid of weird folder names within the archive
mkdir -p $OUTDIR/repacked
for NUM in $(seq 31)
do echo "checking ${NAMES[$NUM]}"
if [[ ${FILES[$NUM]} =~ \.t?gz$ ]] || [[ ${FILES[$NUM]} =~ \.zip$ ]] || [[ ${FILES[$NUM]} =~ \.t?bz2$ ]]
then
	echo "repacking ${NAMES[$NUM]}"
	mkdir $OUTDIR/${NAMES[$NUM]}
	tar -xf $OUTDIR/${FILES[$NUM]} --strip-components=1 -C $OUTDIR/${NAMES[$NUM]}
	tar -czf $OUTDIR/repacked/${NAMES[$NUM]}.tar.gz -C $OUTDIR ${NAMES[$NUM]}
	rm -r $OUTDIR/${NAMES[$NUM]}
else
	echo "skip repacking ${NAMES[$NUM]} - just packing it"
	tar -czf $OUTDIR/repacked/${NAMES[$NUM]}.tar.gz -C $OUTDIR ${FILES[$NUM]}
fi
done


