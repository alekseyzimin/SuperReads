#!/bin/bash
#this pipeline generates genome annotation using hisat2, Stringtie2 and maker
GENOME="genome.fa"
PROT="proteins.fa"
GENOMEFILE="genome.fa"
RNASEQ_PAIRED="paired"
RNASEQ_UNPAIRED="unpaired"
BATCH_SIZE=1000000
UNIPROT="uniprot.fa"
MYPATH="`dirname \"$0\"`"
MYPATH="`( cd \"$MYPATH\" && pwd )`"
PID=$$
export PATH=$MYPATH:$PATH;
set -o pipefail
NUM_THREADS=1
usage="Usage: annotate_maker.sh -t <number of threads> -g <genome fasta file with full path> -p <file containing list of filenames paired Illumina reads from RNAseq experiment, one pair of filenames per line> -u <file containing list of filenames of unpaired Illumina reads from RNAseq experiment, one filename per line> -r <fasta file of protein sequences> -s <uniprot proteins> -v <verbose flag>"
GC=
RC=
NC=
if tty -s < /dev/fd/1 2> /dev/null; then
    GC='\e[0;32m'
    RC='\e[0;31m'
    NC='\e[0m'
fi

trap abort 1 2 15
function abort {
log "Aborted"
rm -rf /dev/shm/tmp_$PID
kill 0
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
if [ $# -lt 1 ];then
  error_exit "$usage"
fi

#parsing arguments
while [[ $# > 0 ]]
do
    key="$1"

    case $key in
        -t|--threads)
            NUM_THREADS="$2"
            shift
            ;;
        -g|--genome)
            GENOMEFILE="$2"
            shift
            ;;
        -p|--paired)
            RNASEQ_PAIRED="$2"
            shift
            ;;
        -r|--proteins)
            PROT="$2"
            shift
            ;;
        -s|--swissprot)
            UNIPROT="$2"
            shift
            ;;
        -u|--unpaired)
            RNASEQ_UNPAIRED="$2"
            shift
            ;;
        -v|--verbose)
            set -x
            ;;
        -h|--help|-u|--usage)
            echo "$usage"
            exit 0
            ;;
        *)
            echo "Unknown option $1"
            exit 1        # unknown option
            ;;
    esac
    shift
done

GENOME=`basename $GENOMEFILE`
#checking is dependencies are installed
for prog in $(echo "ufasta hisat2 stringtie maker gffread gff3_merge snap");do
  which $prog
  if [ $? -gt 0 ];then error_exit "$prog not found the the PATH";fi
done

export SNAP_PATH=`which maker`
export SNAP_PATH=`dirname $SNAP_PATH`/../exe/snap/

#checking inputs
if [ ! -s $RNASEQ_PAIRED ] && [ ! -s $RNASEQ_UNPAIRED ];then
  error_exit "Must specify at lease one non-empty file with filenames of RNAseq reads with -p or -u"
fi

if [ ! -s $UNIPROT ];then
  error_exit "File with uniprot sequences is missing, please supply it with -s <uniprot_file.fa>"
fi

if [ ! -s $GENOMEFILE ];then
  error_exit "File with genome sequence is missing, please supply it with -g <genome_file.fa>"
fi

if [ ! -s $PROT ];then
  error_exit "File with proteing sequences for related species is missing, please supply it with -r <proteins_file.fa>"
fi

#first we align
if [ ! -e align-build.success ];then
  log "building HISAT2 index"
  hisat2-build $GENOMEFILE $GENOME.hst 1>hisat2-build.out 2>&1 && touch align-build.success && rm -f align.success
fi

if [ ! -e align.success ];then
  log "aligning RNAseq reads"
  echo "#!/bin/bash" >hisat2.sh
  if [ -s $RNASEQ_PAIRED ];then
    awk 'BEGIN{n=0}{print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' -1 "$1" -2 "$2" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++}' $RNASEQ_PAIRED >> hisat2.sh
  fi
  if [ -s $RNASEQ_UNPAIRED ];then
    awk 'BEGIN{n=0}{print "if [ ! -s tissue"n".bam ];then hisat2 '$GENOME'.hst --dta -p '$NUM_THREADS' -U "$1" 2>tissue"n".err | samtools view -bhS /dev/stdin > tissue"n".bam.tmp && mv tissue"n".bam.tmp tissue"n".bam; fi"; n++}' $RNASEQ_UNPAIRED >> hisat2.sh
  fi
  bash ./hisat2.sh && touch align.success && rm -f sort.success
fi

if [ ! -e sort.success ];then
  log "sorting alignment files"
  for f in $(ls tissue*.bam);do
    if [ ! -s $f.sorted.bam ];then
      samtools sort -@ $NUM_THREADS -m 1G $f $f.sorted.tmp && mv $f.sorted.tmp.bam $f.sorted.bam
    fi
  done
  touch sort.success && rm -f stringtie.success
fi

if [ ! -e stringtie.success ] && [ -e sort.success ];then
  for f in $(ls tissue*.bam.sorted.bam);do
    if [ ! -s $f.gtf ];then
      stringtie -p $NUM_THREADS $f -o $f.gtf.tmp && mv $f.gtf.tmp $f.gtf
    fi
  done
  INCOUNT=`ls tissue*.bam.sorted.bam|wc -l`
  OUTCOUNT=`ls tissue*.bam.sorted.bam.gtf|wc -l`
  if [ $INCOUNT -eq $OUTCOUNT ];then
    touch stringtie.success
  else
    error_exit "one or more stringtie jobs failed"
  fi
  rm -f split.success
fi

if [ ! -e split.success ] && [ -e stringtie.success ];then 
#now to set up parallel maker we need to split fasta files and transcript files split filenames into batches batch.number
  ufasta sizes -H $GENOMEFILE | perl -ane 'BEGIN{$fn="batch";$index=1;$batch_size=int("'$BATCH_SIZE'");open(FILE,">$fn.$index");$n=0;}{if($n > $batch_size){close(FILE);$index++;open(FILE,">$fn.$index");$n=$F[1];}else{$n+=$F[1];}print FILE $F[0],"\n";}' 
 NUM_BATCHES=`ls batch.* |wc -l`
#now we extract all batches from gtf files and set up directories
  for f in $(seq 1 $NUM_BATCHES);do
    mkdir -p $f.dir
    rm -rf $f.dir/$f.transcripts.fa*
    ufasta extract -f batch.$f $GENOMEFILE > $f.dir/$f.fa
    for t in $(ls tissue*.bam.sorted.bam.gtf);do
      perl -ane 'BEGIN{open(FILE,"batch.'$f'");while($line=<FILE>){chomp($line);$h{$line}=1;}}{print if defined($h{$F[0]})}' $t | gffread -g $f.dir/$f.fa -w $f.dir/$f.transcripts.fa.$t /dev/stdin
    done
    cat $f.dir/$f.transcripts.fa.* | awk 'BEGIN{n=0}{if($1 ~ /^>/){print $1"_"n;n++}else{print $0}}' > $f.dir/$f.transcripts.all.fa.tmp && \
    rm $f.dir/$f.transcripts.fa.* && \
    mv $f.dir/$f.transcripts.all.fa.tmp $f.dir/$f.transcripts.fa
  done
  touch split.success && rm -f maker1.success
fi

if [ ! -e maker1.success ] && [ -e split.success ] && [ -e stringtie.success ];then
  log "running maker"
#first we create a script to run
  maker -CTL
  mkdir -p /dev/shm/tmp_$PID
  echo "#!/bin/bash" > run_maker.sh
  echo "cd \$1 && echo \"running maker in \$PWD\" && maker -cpus 4 -base $GENOME 1>maker.log 2>&1" >> run_maker.sh && \
  chmod 0755 run_maker.sh
  NUM_BATCHES=`ls batch.* |wc -l`
#let's create the appropriate maker ctl files
  for f in $(seq 1 $NUM_BATCHES);do
    sed s,^genome=,genome=$PWD/$f.dir/$f.fa, maker_opts.ctl | \
    sed s,^est=,est=$PWD/$f.dir/$f.transcripts.fa, | \
    sed s,^est2genome=0,est2genome=1, | \
    sed s,^protein2genome=0,protein2genome=1, | \
    sed s,^single_exon=0,single_exon=1, | \
    sed s,^TMP=,TMP=/dev/shm/tmp_$PID, | \
    sed s,^max_dna_len=100000,max_dna_len=1000000, | \
    sed s,^cpus=1,cpus=4, | \
    sed s,^min_contig=1,min_contig=1000, |\
    sed s,^protein=,protein=$PROT, > $f.dir/maker_opts.ctl && \
    cp maker_exe.ctl $f.dir && \
    cp maker_bopts.ctl $f.dir
  done
#running maker
  ls -d *.dir | xargs -P $NUM_THREADS  -I % ./run_maker.sh % && touch maker1.success && rm -rf /dev/shm/tmp_$PID
fi

if [ -e maker1.success ] && [ ! -e functional.success ];then 
  log "concatenating outputs"
  NUM_BATCHES=`ls batch.* |wc -l`
  for f in $(seq 1 $NUM_BATCHES);do
    (cd $f.dir/$GENOME.maker.output && gff3_merge -g -d ${GENOME}_master_datastore_index.log -o ../$f.gff && fasta_merge -d ${GENOME}_master_datastore_index.log -o ../$f.fasta )
  done
  gff3_merge -o $GENOME.gff.tmp *.dir/*.gff && mv $GENOME.gff.tmp $GENOME.gff && \
  cat *.dir/*.all.maker.proteins.fasta > $GENOME.proteins.fasta.tmp  && mv $GENOME.proteins.fasta.tmp $GENOME.proteins.fasta && \
  cat *.dir/*.all.maker.transcripts.fasta > $GENOME.transcripts.fasta.tmp && mv $GENOME.transcripts.fasta.tmp $GENOME.transcripts.fasta && \
  log "maker output is in $GENOME.gff $GENOME.proteins.fasta $GENOME.transcripts.fasta" && \
  log "Performing functional annotation" && \
  makeblastdb -in $UNIPROT -input_type fasta -dbtype prot -out uniprot && \
  blastp -db uniprot -query $GENOME.proteins.fasta -out  $GENOME.maker2uni.blastp -evalue 0.000001 -outfmt 6 -num_alignments 1 -seg yes -soft_masking true -lcase_masking -max_hsps 1 -num_threads $NUM_THREADS && \
  maker_functional_gff $UNIPROT $GENOME.maker2uni.blastp $GENOME.gff > $GENOME.functional_note.gff.tmp && mv $GENOME.functional_note.gff.tmp $GENOME.functional_note.gff && \
  maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.proteins.fasta > $GENOME.functional_note.proteins.fasta.tmp  && mv $GENOME.functional_note.proteins.fasta.tmp $GENOME.functional_note.proteins.fasta && \
  maker_functional_fasta $UNIPROT $GENOME.maker2uni.blastp $GENOME.transcripts.fasta > $GENOME.functional_note.transcripts.fasta.tmp  && mv $GENOME.functional_note.transcripts.fasta.tmp $GENOME.functional_note.transcripts.fasta && \
  touch functional.success
fi

if [ -e functional.success ];then
  log "Output annotation is in $GENOME.functional_note.gff $GENOME.functional_note.proteins.fasta $GENOME.functional_note.transcripts.fasta" 
fi



