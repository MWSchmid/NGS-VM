#!/bin/bash
#
# ChIPseq_runsTest.sh -- compare two samples with a custom made runsTest.
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
curDir=$(dirname "$0")

## defaults
inputDir=""
outputDir=""
prefix=""
inputFile=""
inputFileReference=""
runsTestScript="$curDir/edgeTestChXPseq.py"
postProcScript="$curDir/processChXPedgeTest.R"
flankSize=100000
baseSteps=100
fPositive=0.7
numReps=5
randomized=""
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
  $me [options] [INDIR] [OUTDIR] [OUTPREFIX] [BIGWIGFILE_TEST] [BIGWIGFILE_CONTROL]
Compare two normalized (!) bigWigs with a custom-made edge-test. IMPORTANT:
It's a one-sided test.
Output:
<OUTPREFIX>.ET.<flankSize>_<baseSteps>_<fPositive>_<numReps>.txt: chrom, start, end, leftBorder, rightBorder, positiveFraction
<OUTPREFIX>.ET.<flankSize>_<baseSteps>_<fPositive>_<numReps>.bed: chrom, start, end, peak_<number>, positiveFraction*abs(LB)*abs(RB)/max(...)*100
Arguments:
INDIR: Directory with the input file (<extension/type>).
OUTDIR: Directory in which all output will be store.
OUTPREFIX: Prefix for output. The output file will be named <OUTPREFIX>.ET.<flankSize>_<baseSteps>_<fPositive>_<numReps>.bed/txt
BIGWIGFILE_TEST: Name of the normalized bigWig file of the test sample.
BIGWIGFILE_CONTROL: Name of the normalized bigWig file of the control sample.
Options:
  -v        enable verbose logging (no effect)
  -h        print this help text
  -t        number of available threads (no effect)
  -m        amount of memory to be allocated (no effect)
  -s        path to the edgeTestChXPseq.py script
  -p        path to the processChXPedgeTest.R script
  -b        take only every Xth base (default: 100)
  -f        minimal fraction of positive differences within a candidate region (default: 0.7)
  -n        number of random sets (default: 5)
  -l        size of the regions left and right of the border to check (default 100'000)
  -r        do the same as always, but randomize the test set as well (default: off) - NOT VISIBLE IN THE FILE NAME PER SE

Dependencies:
sudo easy_install ngslib
sudo pip install skidmarks
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

short_opts='hvt:m:s:p:b:f:n:l:r'
long_opts='help,verbose,threads,memory,script,postProcScript,baseSteps,fPositive,numReps,flankSize,randomized'

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
        --postProcScript|-p)   shift; postProcSript=$1 ;;
        --script|-s)   shift; runsTestScript=$1 ;;
        --baseSteps|-b)shift; baseSteps=$1 ;;
        --fPositive|-f)shift; fPositive=$1 ;;
        --numReps|-n)  shift; numReps=$1 ;;
        --flankSize|-l)shift; flankSize=$1 ;;
        --randomized|-r) randomized=" --randomized" ;;
        --threads|-t)  shift; threads=$1 ;;
        --memory|-m)   shift; memory=$1 ;;	
        --verbose|-v)  verbose="--verbose" ;;
        --help|-h)     usage; exit 0 ;;
        --)            shift; break ;;
    esac
    shift
done

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
inputFileReference=$1

## main
echo "=== ${me}: Starting at `date '+%Y-%m-%d %H:%M:%S'`"

require_command Rscript
require_command ${runsTestScript} 

# checking input
input_exists ${inputDir}/${inputFile}
input_exists ${inputDir}/${inputFileReference}
input_exists ${postProcScript}

# run script
outfileName="${prefix}.ET.${flankSize}_${baseSteps}_${fPositive}_${numReps}"
command="${runsTestScript} ${inputDir}/${inputFile} ${inputDir}/${inputFileReference} --baseSteps ${baseSteps} --fPositive ${fPositive} --numReps ${numReps} --flankSize ${flankSize}${randomized} > ${outputDir}/${outfileName}.txt"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"

command="Rscript ${postProcScript} ${outputDir}/${outfileName}.txt ${outputDir}/${outfileName}.bed"
echo "=== ${me}: Running: ${command}"
eval $command
rc=$?
echo "=== ${me}: Command ended with exit code $rc"
#awk -v OFS="\t" -v TST=${inputFile//.bw} -v CNT=${inputFileReference//.bw} '{if ($5 < 0) {print $1,$2,$3,CNT,$5} else {print $1,$2,$3,TST,$5} }' ${outputDir}/${prefix}.RT.${winLen}.txt > ${outputDir}/${prefix}.RT.${winLen}.bed

## Checking output
output_exists "${outputDir}/${outfileName}.txt"
output_exists "${outputDir}/${outfileName}.bed"

## All done.
echo "=== ${me}: Script done at `date '+%Y-%m-%d %H:%M:%S'`."

exit $rc
