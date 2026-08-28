#!/bin/bash
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"

FASTA="genome.fa"
PSAURON="psauron_score.csv"
SCOREFILE="scores.txt"
MULT=`perl -e 'print exp(1)'`
OUTDEV="out.err"
MIN_CDS=200

GC=
RC=
NC=
ALT_ENABLED=0
ALT_OUTPUT=
ALT_ARGS=()
LOCAL_K_ENABLED=0
LOCAL_K_OUTPUT=
LOCAL_K_ARGS=()
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
kill -9 0
exit 1
}

log () {
    dddd=$(date)
    echo -e "${GC}[$dddd]${NC} $@"
}

function error_exit {
    dddd=$(date)
    echo -e "${RC}[$dddd]${NC} $1" >&2
    exit "${2:-1}"
}

function usage {
echo "Usage:"
echo "uniann.sh [arguments]"
echo "-f file sith a single fasta sequence"
echo "-m multiplier to rescale probabilities before applying log, must be a number between 1 and 100, default: exp(1)"
echo "-s file with scores for starts, stops and splice sites, 1-based coordinates"
echo "-p psauron score file"
echo "-n (flag) do not save Viterbi matrix in out.err"
echo "--alternative-splicing enable anchored alternative-transcript search"
echo "--alt-output file write the additional UniAnn-format GFF3 annotation"
echo "--alt-report file write the optional diagnostic TSV"
echo "--alt-max-distance bp maximum search distance from each anchor (default 1000000)"
echo "--alt-min-score-delta score filter complete alternatives below this delta"
echo "--alt-include-incomplete include incomplete candidates in the TSV only"
echo "--alt-debug print alternative-search diagnostics"
echo "--local-k-best [K] enable gene-local k-best transcript search (default K=5)"
echo "--local-k-output file write the k-best GFF3 annotation"
echo "--local-k-report file write the diagnostic TSV"
}

#parsing arguments
if [[ $# -eq 0 ]];then
usage
exit 1
fi

while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -f|--fasta)
            FASTA="$2"
            shift
            ;;
        -p|--psauron)
            PSAURON="$2"
            shift
            ;;
        -m|--mult)
            MULT="$2"
            shift
            ;;
        -n|--noviterbi)
            OUTDEV="/dev/null"
            ;;
        -s|--scores)
            SCOREFILE="$2"
            shift
            ;;
        --alternative-splicing)
            ALT_ENABLED=1
            ;;
        --alt-output|--alt-report|--alt-max-distance|--alt-min-score-delta)
            if [[ $# -lt 2 ]]; then
                error_exit "Missing value for $1"
            fi
            ALT_ARGS+=("$1" "$2")
            if [[ "$1" == "--alt-output" ]]; then
                ALT_OUTPUT="$2"
            fi
            shift
            ;;
        --alt-include-incomplete|--alt-debug)
            ALT_ARGS+=("$1")
            ;;
        --local-k-best)
            LOCAL_K_ENABLED=1
            LOCAL_K_ARGS+=("$1")
            if [[ $# -gt 1 && "$2" =~ ^[0-9]+$ ]]; then
                LOCAL_K_ARGS+=("$2")
                shift
            fi
            ;;
        --local-k-output|--local-k-report)
            if [[ $# -lt 2 ]]; then
                error_exit "Missing value for $1"
            fi
            LOCAL_K_ARGS+=("$1" "$2")
            if [[ "$1" == "--local-k-output" ]]; then
                LOCAL_K_OUTPUT="$2"
            fi
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            usage
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

if [[ ! -s $FASTA ]];then
echo "Input file $FASTA not found or not specified!"
usage
exit 1
fi

if [[ ! -s $SCOREFILE ]];then
echo "Input training annotation file $TRAINING_A not found or not specified!"
usage
exit 1
fi

if [[ ! -s $PSAURON ]];then
echo "Input file of PSAURON scores $PSAURON not found or not specified!"
usage
exit 1
fi

#this produces out.ps.txt
log "Preprocessing inputs" && \
if [ ! -s out.ps.txt ];then
  $MYPATH/preprocess_psauron_scores.pl $FASTA $PSAURON
fi

FACTOR=`cat $SCOREFILE |perl -ane '{$F[5]=$F[6] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}'|perl -ane 'BEGIN{$max_score=0}{if($F[3] eq "donor"){$score=log($F[5]*'$MULT'+1e-10);$max_score=$score if($score>$max_score);}}END{die("Incorrect scores in the input file: must be between 0 and 1!") if($max_score<=0);print int(1000/$max_score+0.5)}'` && \
log "Multiplier is $MULT, Factor is $FACTOR" && \
#this produces out.atg.txt out.gt.txt and out.ag.txt out.stop
cat $SCOREFILE | \
  perl -ane '{$F[5]=$F[-1] if($#F>5);;print join("\t",@F),"\n" if($F[2] eq "+");}' | \
  tee >(perl -ane '{if($F[3] eq "donor"){$score=log($F[5]*'$MULT'+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.gt.txt) |\
  tee >(perl -ane '{if($F[3] eq "acceptor"){$score=log($F[5]*'$MULT'+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.ag.txt) |\
  tee >(perl -ane '{if($F[3] eq "start"){$score=log($F[5]*('$MULT')+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.atg.txt) | \
  perl -ane '{if($F[3] eq "stop"){$score=log($F[5]*('$MULT')+1e-10)*'$FACTOR'; print $F[1]-1,"\t$score\n"}}' > out.stop.txt && \

if [[ $ALT_ENABLED -eq 1 ]]; then
  ALT_ARGS=("--alternative-splicing" "${ALT_ARGS[@]}")
  if [[ -z "$ALT_OUTPUT" ]]; then
    ALT_ARGS+=("--alt-output" "$FASTA.uniann.alternative_splicing.gff")
  fi
elif [[ ${#ALT_ARGS[@]} -gt 0 ]]; then
  error_exit "Alternative-splicing options require --alternative-splicing"
fi

if [[ $LOCAL_K_ENABLED -eq 1 ]]; then
  # -n suppresses the Viterbi matrix; without it the DP/BT dump is kept
  if [[ "$OUTDEV" == "/dev/null" ]]; then
    LOCAL_K_ARGS+=("--no-dp-dump")
  fi
  if [[ -z "$LOCAL_K_OUTPUT" ]]; then
    LOCAL_K_ARGS+=("--local-k-output" "$FASTA.uniann.local_k_best.gff")
  fi
elif [[ ${#LOCAL_K_ARGS[@]} -gt 0 ]]; then
  error_exit "Local k-best options require --local-k-best"
fi

log "Building gene models" && \
#also enforce MIN_CDS
$MYPATH/uniann "$FASTA" out.ps.txt out.gt.txt out.ag.txt out.atg.txt out.stop.txt "${ALT_ARGS[@]}" "${LOCAL_K_ARGS[@]}" 2>$OUTDEV | \
  gffread --tlf |\
  perl -F'\t' -ane '{
    if($F[8]=~/exonCount=(1|2);exons=(\S+);CDS=(\d+):(\d+);CDSphase=\d/){
      $cdslen=0;
      @exons=split(/,/,$2);
      foreach $e(@exons){
        $cdslen+=$2-$1+1 if($e=~/(\d+)-(\d+)/);
      }
      print if($cdslen>'$MIN_CDS');
    }else{
      print;
    }
  }' |\
  gffread > $FASTA.gff.tmp && \
mv $FASTA.gff.tmp $FASTA.uniann.gff

echo "Output gff file is $FASTA.uniann.gff"
