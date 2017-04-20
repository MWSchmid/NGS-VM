#!/bin/bash
#
# preprocessIlluminaFastq_SE.sh -- script for preprocessing of Illumina reads (SE)
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
threads=1
memory=4
maxMem=4
remDupReads=""

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
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [FASTQ]

NOTE: not yet tested

This script will run:
1) fastqc: standard quality checks
2) trimGalore: Illumina adapter removal and low quality clipping
3) fqtrim: low complexity filter and duplicate collapsing if requested
Arguments:
INDIR: Directory with the input file.
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. The output files will be named <OUTPREFIX>.tr(.fi).fq.gz.
FASTQ: Name of the fastq(.gz) file.
Options:
  -v            Enable verbose logging (no effect)
  -h            Print this help text
  -t		Number of available threads
  -m            Amount of memory to be allocated (per core, in GB)
  -d		Collapse duplicate reads
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

## parse command-line

short_opts='hvt:m:d'
long_opts='help,verbose,threads,memory,remDupReads'

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
    	--remDupReads|-d)	       remDupReads=" -C" ;;
	--threads|-t)  		shift; threads=$1 ;;
	--memory|-m)  		shift; memory=$1 ;;	
        --verbose|-v) 		verbose='--verbose' ;;
        --help|-h)    		usage; exit 0 ;;
        --)           		shift; break ;;
    esac
    shift
done

maxMem=$(($threads * $memory))

## sanity checks

# the required arguments must be present 
if [ $# -lt 4 ]; then
    die $EX_USAGE "Missing required arguments. Type '$me --help' to get usage help."
fi

inputDir=$1
shift
outputDir=$1
shift
prefix=$1
shift
inputFile=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command fastqc
require_command trim_galore
require_command fqtrim

# checking input
input_exists ${inputDir}/${inputFile}

# run fastqc
command="fastqc -t ${threads} -o ${outputDir} ${inputDir}/${inputFile}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

# run trim galore
command="trim_galore --illumina -o ${outputDir} ${inputDir}/${inputFile}"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
mv "${outputDir}/${inputFile%%.*}_trimmed.fq.gz" "${outputDir}/${prefix}.tr.fq.gz"  

# run fqtrim
input_exists "${outputDir}/${prefix}.tr.fq.gz"  
command="fqtrim -A${remDupReads} -D -l 30 -r ${outputDir}/${prefix}_fqtrimReport.txt -p ${threads} -o fi.fq.gz ${outputDir}/${prefix}.tr.fq.gz" 
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

# moving output files from fqtrim
curFile="${prefix}.tr.fi.fq.gz"
output_exists "$curFile"
mv ${curFile} ${outputDir}/${curFile}
output_exists "${outputDir}/${curFile}"

# run fastqc again
command="fastqc -t ${threads} -o ${outputDir} ${outputDir}/${prefix}.tr.fi.fq.gz"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
