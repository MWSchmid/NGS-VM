#!/bin/bash
#
# mergeAndSortSplitPEreads.sh -- re-pair separated PE reads.
# NOTE: subread utility has sth like this
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
inputFileReverse=""
solidOption=""
pyScript="$HOME/mergeAndSortSplitPEreads.py"
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
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [FASTQ-forward] [FASTQ-reverse] [FASTQ-unpaired] [...]
This script will merge two PE read files of different size and order:
mergeAndSortSplitPEreads.py
Output:
<OUTPREFIX>_R<1/2>.repaired.fq.gz
<OUTPREFIX>_MIXED.unpairable.fq.gz
Arguments:
INDIR: Directory with the input files.
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. The output files will be named <OUTPREFIX>_<MIXED/R<1/2>>.<unpairable/repaired>.fq.gz
FASTQ-forward: Name of the fastq(.gz) file with the forward reads.
FASTQ-reverse: Name of the fastq(.gz) file with the reverse reads.
FASTQ-unpaired: Name of a fastq(.gz) file with unpaired reads. One can supply as many files as wanted.
The reads in these files will be added to the *.unpairable.fq.gz.
Options:
  -v            Enable verbose logging (no effect)
  -h            Print this help text
  -t		Number of available threads (no effect)
  -m            Amount of memory to be allocated (per core, in GB, no effect)
  -s		Path to the python script (default: ${pyScript})
  -c		Specifies that the read come from a SOLID machine
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

short_opts='hvt:m:s:c'
long_opts='help,verbose,threads,memory,pyScript,solid'

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
	--solid|-c)    solidOption=" --solid" ;;
        --pyScript|-s) shift; pyScript=$1 ;;
	--threads|-t)  shift; threads=$1 ;;
	--memory|-m)   shift; memory=$1 ;;	
        --verbose|-v)  verbose='--verbose' ;;
        --help|-h)     usage; exit 0 ;;
        --)            shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 5 ]; then
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
inputFileReverse=$1
shift
inputFilesUnpaired=("$@")
numUnpairedInputFiles=${#inputFilesUnpaired[@]}

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command python

# checking input
input_exists ${inputDir}/${inputFile}
input_exists ${inputDir}/${inputFileReverse}
input_exists $pyScript

if [ $numUnpairedInputFiles > 0 ]; then
  for curFile in "${inputFilesUnpaired[@]}"; do
  input_exists ${inputDir}/${curFile}
  done
fi

# merge the reads
command="python ${pyScript} ${inputDir}/${inputFile} ${inputDir}/${inputFileReverse}\
	 ${outputDir}/${prefix}_R1.repaired.fq.gz ${outputDir}/${prefix}_R2.repaired.fq.gz\
	 ${outputDir}/${prefix}_TEMP.unpairable.fq.gz --maxReads 5000000${solidOption}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Script ended with exit code $rc"

# merge the unpairable and the ab initio unpaired files (if given)
if [ $numUnpairedInputFiles > 0 ]; then
inputFilesUnpairedList=$(printf "${inputDir}/%s " "${inputFilesUnpaired[@]}")
inputFilesUnpairedList=${inputFilesUnpairedList::-1}
zcat ${outputDir}/${prefix}_TEMP.unpairable.fq.gz ${inputFilesUnpairedList} | pigz -p ${compThreadsMerge} > "${outputDir}/${prefix}_MIXED.unpairable.fq.gz"
else
mv "${outputDir}/${prefix}_TEMP.unpairable.fq.gz ${outputDir}/${prefix}_MIXED.unpairable.fq.gz"
fi

# checking output
output_exists "${outputDir}/${prefix}_R1.repaired.fq.gz" 
output_exists "${outputDir}/${prefix}_R2.repaired.fq.gz"
output_exists "${outputDir}/${prefix}_MIXED.unpairable.fq.gz"
remove_if_present "${outputDir}/${prefix}_TEMP.unpairable.fq.gz"
remove_if_present "${outputDir}/sorted_${inputFile}"
remove_if_present "${outputDir}/sorted_${inputFileReverse}"

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
