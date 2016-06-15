#!/bin/bash
#
# ChIPseq_SICER.sh -- runs SICER to identify broad peaks (histone mods).
# with or without reference sample.
#
# Authors: Marc W. Schmid <marcschmid@gmx.ch>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
#

me=$(basename "$0")

## defaults
inputDir=""
outputDir=""
prefix=""
inputFile=""
inputFileReference=""
winSize=200
gapSize=200
fragSize=170
genomeID=""
redundancyThreshold=1 # no duplicate reads
effectiveGenomeFraction="0.9"
evalue=100
threads=1
memory=4
maxMem=4

## Exit status codes (following <sysexits.h>)
EX_OK=0			# successful termination
EX_USAGE=64		# command line usage error
EX_DATAERR=65		# data format error
EX_NOINPUT=66		# cannot open input
EX_NOUSER=67		# addressee unknown
EX_NOHOST=68		# host name unknown
EX_UNAVAILABLE=69	# service unavailable
EX_SOFTWARE=70		# internal software error
EX_OSERR=71		# system error (e.g., can't fork)
EX_OSFILE=72		# critical OS file missing
EX_CANTCREAT=73		# can't create (user) output file
EX_IOERR=74		# input/output error
EX_TEMPFAIL=75		# temp failure; user is invited to retry
EX_PROTOCOL=76		# remote error in protocol
EX_NOPERM=77		# permission denied
EX_CONFIG=78		# configuration error

## helper functions

function die () {
    rc="$1"
    shift
    (echo -n "$me: ERROR: ";
        if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
    exit $rc
}

## usage info

usage () {
    cat <<__EOF__
Usage:
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [BAMFILE_TEST] [BAMFILE_CONTROL] [GENOME-ID]
Run SICER on a test and reference/control sample to identify broad regions (histone modifications).
Arguments:
INDIR: Directory with the input file (<extension/type>).
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. The output file will be named <OUTPREFIX>.SICER.<winSize>.<gapSize>.bed
BAMFILE_TEST: Name of the bam file of the test sample.
BAMFILE_CONTROL: Name of the bam file of the control sample.
GENOME-ID: The genome-ID used in the "GenomeData.py" script.
Options:
  -v            Enable verbose logging (no effect)
  -h            Print this help text
  -t		Number of available threads (no effect)
  -m            Amount of memory to be allocated (per core, in GB)
  -w		winSize passed on to SICER (default 200)
  -g		gapSize passed on to SICER (default 200 - ok for H3K4me3, take 600 for H3K27me3)
  -f		fragSize passed on to SICER (default 170)
  -e		effectiveGenomeFraction (default 0.9) - fraction of the genome which is alignable
Dependencies:
you maybe need to add your genome to:
/usr/local/bin/SICER/lib/GenomeData.py
__EOF__
}

warn () {
  (echo -n "$me: WARNING: ";
      if [ $# -gt 0 ]; then echo "$@"; else cat; fi) 1>&2
}

have_command () {
  type "$1" >/dev/null 2>/dev/null
}

require_command () {
  if ! have_command "$1"; then
    die $EX_UNAVAILABLE "Could not find required command '$1' in system PATH. Aborting."
  fi
}

is_absolute_path () {
    expr match "$1" '/' >/dev/null 2>/dev/null
}

input_exists () {
  echo -n "Checking input file ${1}..."
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_NOINPUT "Could not find input file: ${1}."
  fi
}

output_exists () {
  echo -n "Checking output file ${1}... "
  if [ -e $1 ]; then
    echo "[ok]"
  else
    echo "[failed]"
    die $EX_OSFILE "Could not find output file: ${1}."
  fi
}

remove_if_present () {
  if [ -e $1 ]; then
    rm $1
  fi
}

## parse command-line

short_opts='hvt:m:w:g:f:e:'
long_opts='help,verbose,threads,memory,winSize,gapSize,fragSize,effGenFrac'

getopt -T > /dev/null
rc=$?
if [ "$rc" -eq 4 ]; then
    # GNU getopt
    args=$(getopt --name "$me" --shell sh -l "$long_opts" -o "$short_opts" -- "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    # use 'eval' to remove getopt quoting
    eval set -- $args
else
    # old-style getopt, use compatibility syntax
    args=$(getopt "$short_opts" "$@")
    if [ $? -ne 0 ]; then
        die $EX_USAGE "Type '$me --help' to get usage information."
    fi
    set -- $args
fi

while [ $# -gt 0 ]; do
    case "$1" in
        --winSize|-w)    shift; winSize=$1 ;;
        --gapSize|-g)    shift; gapSize=$1 ;;
        --fragSize|-f)   shift; fragSize=$1 ;;
        --effGenFrac|-e) shift; effectiveGenomeFraction=$1 ;;
	--threads|-t)  shift; threads=$1 ;;
	--memory|-m)   shift; memory=$1 ;;	
        --verbose|-v)  verbose='--verbose' ;;
        --help|-h)     usage; exit 0 ;;
        --)            shift; break ;;
    esac
    shift
done

## sanity checks

# the required arguments must be present 
if [ $# -lt 6 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

inputDir=$1
shift
outputDir=$1
shift
prefix=$1
shift
inputFile=$1
shift
inputFileReference=$1
shift
genomeID=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command SICER.sh
require_command bedtools

# checking input
input_exists ${inputDir}/${inputFile}
input_exists ${inputDir}/${inputFileReference}

# convert bam to bed
echo "=== converting bam to bed"
bedtools bamtobed -i ${inputDir}/${inputFile} > ${outputDir}/${inputFile//.bam}.bed
bedtools bamtobed -i ${inputDir}/${inputFileReference} > ${outputDir}/${inputFileReference//.bam}.bed

# run sicer
command="SICER.sh ${inputDir} ${outputDir}/${inputFile//.bam}.bed ${outputDir}/${inputFileReference//.bam}.bed ${outputDir} ${genomeID} ${redundancyThreshold} ${winSize} ${fragSize} ${effectiveGenomeFraction} ${gapSize} ${evalue}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

awk -v OFS="\t" -v PREFIX=${prefix} '{if ($8 < 0.01) {print $1,$2,$3,PREFIX"_peak_"NR,$6,"."}}' "${outputDir}/${inputFile//.bed}-W${winSize}-G${$gapSize}-islands-summary-FDR100" | sort -k1,1 -k2,2n > ${outputDir}/${prefix}.SICER.${winSize}.${gapSize}.bed

## Checking output
output_exists "${outputDir}/${prefix}.SICER.${winSize}.${gapSize}.bed"

## remove unused output
rm "${outputDir}/${inputFile//.bed}-*"
remove_if_present ${outputDir}/${inputFile//.bam}.bed
remove_if_present ${outputDir}/${inputFileReference//.bam}.bed

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
