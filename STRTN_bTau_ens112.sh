#!/bin/bash

PROGNAME="$( basename $0 )"

# Usage
function usage() {
    cat <<EOS >&2
Usage: ${PROGNAME} [-o <output>] [-g <genome (required)>] [-a <annotation>] [-b <path (required)>] [-i <path (required)>] [-w <path (required)>] [-p <path (required)>]
Options:
  -o, --out               Output file name. (default: OUTPUT)
  -g, --genome            Genome (hg19/hg38/mm9/mm10/mm39/canFam3/canFam6/bosTau9). Required!
  -a, --annotation        Gene annotation (ref{RefSeq}/ens{Ensembl}/kg{UCSC KnownGenes}/wgEncodeGencodeBasis*) for QC and counting. Default : ref. NOTE: no Ensembl for hg38 & mm10 & mm39 & canFam6, no KnownGenes for canFam3 & canFam6 & bosTau9, no Gencode for mm9 & canFam3 & canFam6 & bosTau9.
  -b, --basecalls         /PATH/to/the Illumina basecalls directory. Required!
  -i, --index             /PATH/to/the directory and basename of the HISAT2 index. Fasta file has to be 'basename.fasta'. Required!
  -w, --working           /PATH/to/the working directory. Required!
  -p, --picardhome	  /PATH/to/the picard.jar. Required!
  -c, --center            The name of the sequencing center that produced the reads. (default: CENTER)
  -r, --run               The barcode of the run. Prefixed to read names. (default: RUNBARCODE)
  -s, --structure         Read structure (default: 8M3S75T6B)
  -d, --dta               Downstream-transcriptome-assembly for HISAT2, which is useful for TFE-based analysis but leads to fewer alignments with short-anchors.
  -h, --help              Show usage.
  -v, --version           Show version.
EOS
    exit 1
}

function version() {
      cat <<EOS >&2
STRT-N NextSeq-automated-pipeline ver2022
EOS
  exit 1
}

# Default parameters
OUTPUT_NAME=OUTPUT
run_VALUE=RUNBARCODE
center_VALUE=FUGU
READ_STRUCTURE=8M3S75T6B
IF_DTA=true

# Parameter settings
PARAM=()
for opt in "$@"; do
    case "${opt}" in
	'-o' | '--out' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    OUTNAME=true
	    OUTPUT_NAME="$2"
	    shift 2
	    ;;
	'-g' | '--genome' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    elif [[ "$2" =  "hg19" ]]; then
		GENOME=true
		GENOME_VALUE="hg19"
		shift 2
	    elif [[ "$2" =  "hg38" ]]; then
		GENOME=true
		GENOME_VALUE="hg38"
		shift 2
	    elif [[ "$2" =  "mm9" ]]; then
		GENOME=true
		GENOME_VALUE="mm9"
		shift 2
	    elif [[ "$2" =  "mm10" ]]; then
		GENOME=true
		GENOME_VALUE="mm10"
		shift 2
	    elif [[ "$2" = "mm39" ]]; then
		GENOME=true
		GENOME_VALUE="mm39"
		shift 2
	    elif  [[ "$2" =  "canFam3" ]]; then
		GENOME=true
		GENOME_VALUE="canFam3"
		shift 2
	    elif  [[ "$2" =  "canFam6" ]]; then
		GENOME=true
		GENOME_VALUE="canFam6"
		shift 2
	    elif  [[ "$2" =  "bosTau9" ]]; then
		GENOME=true
		GENOME_VALUE="bosTau9"
		shift 2
	    else
		usage
		exit 1
	    fi
	    ;;
	'-a' | '--annotation' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    elif [[ "$2" =  "ref" ]]; then
		ANNO=true
		ANNO_VALUE="ref"
		shift 2
	    elif [[ "$2" =  "kg" ]]; then
		ANNO=true
		ANNO_VALUE="kg"
		shift 2
	    elif [[ "$2" =  "ens" ]]; then
		ANNO=true
		ANNO_VALUE="ens"
		shift 2
	    elif [[ "$2" =  wgEncodeGencodeBasic* ]]; then
		ANNO=true
		ANNO_VALUE=$2
		shift 2
	    else
		usage
		exit 1
	    fi
	    ;;
	'-b' | '--basecalls' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    BaseCallsDir=true
	    BaseCallsDir_PATH="$2"
	    shift 2
	    ;;
	'-i' | '--index' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    Index=true
	    Index_PATH="$2"
	    shift 2
	    ;;
	'-w' | '--working' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    WorkingDir=true
	    WorkingDir_PATH="$2"
	    shift 2
	    ;;
	'-p' | '--picardhome' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    PICARDHOME=true
	    PICARDHOME_PATH="$2"
	    shift 2
	    ;;
	'-c' | '--center' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    center_VALUE="$2"
	    shift 2
	    ;;
	'-r' | '--run' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    run_VALUE="$2"
	    shift 2
	    ;;
	'-s' | '--structure' )
	    if [[ -z "$2" ]] || [[ "$2" =~ ^-+ ]]; then
		echo "${PROGNAME}: option requires an argument -- $( echo $1 | sed 's/^-*//' )" 1>&2
		exit 1
	    fi
	    READ_STRUCTURE="$2"
	    shift 2
	    ;;
	'-d' | '--dta' )
	    IF_DTA=true; shift
	    ;;
	'-h' | '--help' )
	    usage
	    ;;
	'-v' | '--version' )
	    version
	    ;;
	'--' | '-' )
	    shift
	    PARAM+=( "$@" )
	    break
	    ;;
	-* )
	    echo "${PROGNAME}: illegal option -- '$( echo $1 | sed 's/^-*//' )'" 1>&2
	    exit 1
	    ;;
    esac
done

if [[ -n "${PARAM[@]}" ]]; then
    usage
fi

[ "${GENOME}" != "true" ] && usage
[ "${BaseCallsDir}" != "true" ] && usage
[ "${Index}" != "true" ] && usage
[ "${ANNO}" != "true" ] && ANNO_VALUE=ref
[ "${WorkingDir}" != "true" ] && usage
[ "${PICARDHOME}" != "true" ] && usage

# Make temporary and output directory
mkdir tmp
mkdir out

# To correctly separate and analyze sequencing data, preparation of barcodes and creation of parameters file are important.
# Preparation for barcodes and creating necessary files for demultiplexing sequencing data using barcodes
ALL_LINES=`cat src/barcode.txt | wc -l`
NLINES=`expr $ALL_LINES \- 1`

for i in `seq 1 $NLINES`
do
    echo -e ${OUTPUT_NAME}_${i}_Lane1.bam"\t"${OUTPUT_NAME}_${i}_Lane1"\t"${OUTPUT_NAME}_${i}_Lane1 >> tmp/out
done
paste tmp/out <(awk 'NR>1{print $1}' src/barcode.txt) | cut -f 1-4 > tmp/out2 && rm tmp/out
echo -e ${OUTPUT_NAME}_non-indexed_Lane1.bam"\t"${OUTPUT_NAME}_non-indexed_Lane1"\t" ${OUTPUT_NAME}_non-indexed_Lane1"\t"N >> tmp/out2
echo -e OUTPUT"\t"SAMPLE_ALIAS"\t"LIBRARY_NAME"\t"BARCODE_1  | cat - tmp/out2 > library.param.lane1 && rm tmp/out2

# Determining the number of lanes and creating parameter files for each lane
nlanes=`ls -l ${BaseCallsDir_PATH} | grep ^d | wc -l`
for i in `seq 2 $nlanes`
do
    sed -e "s/Lane1/Lane${i}/g" library.param.lane1 > library.param.lane${i}
done

# Convert BCL files to BAM files to separate the reads from each sample based on unique barcodes 
for i in `seq 1 $nlanes`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar ExtractIlluminaBarcodes \
	 BASECALLS_DIR=${BaseCallsDir_PATH}/ \
	 LANE=${i} \
	 READ_STRUCTURE=${READ_STRUCTURE} \
	 BARCODE_FILE=src/barcode.txt  \
	 METRICS_FILE=metrics_output_lane${i}.txt ;
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar IlluminaBasecallsToSam \
	 BASECALLS_DIR=${BaseCallsDir_PATH}/ \
	 LANE=${i} \
	 READ_STRUCTURE=${READ_STRUCTURE} \
	 RUN_BARCODE=${run_VALUE} \
	 IGNORE_UNEXPECTED_BARCODES=true \
	 LIBRARY_PARAMS=library.param.lane${i} \
	 SEQUENCING_CENTER=${center_VALUE} \
	 INCLUDE_NON_PF_READS=false
done

rm library.param.lane*
mkdir out/ExtractIlluminaBarcodes_Metrics && mv metrics_output_lane*.txt out/ExtractIlluminaBarcodes_Metrics

# Make the fasta reference / sequence dictionary if they do not exist.
if [[ ! -e ${Index_PATH}.fasta ]]; then
    hisat2-inspect ${Index_PATH} > ${Index_PATH}.fasta
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar CreateSequenceDictionary R=${Index_PATH}.fasta O=${Index_PATH}.dict
fi
if [[ ! -e ${Index_PATH}.dict ]]; then
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar CreateSequenceDictionary R=${Index_PATH}.fasta O=${Index_PATH}.dict
fi

# Mapping by HISAT2 and merging with the original unaligned BAM files to generate UMI-annotated BAM files
mkdir tmp/UMI
mkdir out/HISAT2_Metrics

if [ $IF_DTA = true ]; then
    for file in *.bam
    do
	name=$(basename $file .bam)
	echo $name >> out/HISAT2_Metrics/Alignment-summary.txt
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SortSam \
	       I=$file \
	       O=tmp/.unmapped.sorted.bam \
	       SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SamToFastq \
	       I=tmp/.unmapped.sorted.bam\
	       F=tmp/.tmp.fastq \
	    | hisat2 -p 8 --dta -x ${Index_PATH} \
		     -U tmp/.tmp.fastq -S /dev/stdout \
		     2>> out/HISAT2_Metrics/Alignment-summary.txt  \
	    | java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SortSam \
		     I=/dev/stdin \
		     O=tmp/.mapped.sorted.sam \
		     SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar MergeBamAlignment \
	       ATTRIBUTES_TO_RETAIN=XS \
	       UNMAPPED=tmp/.unmapped.sorted.bam  \
	       ALIGNED=tmp/.mapped.sorted.sam \
	       O=tmp/UMI/$name.umi.bam \
	       R=${Index_PATH}.fasta
    done
else
    for file in *.bam
    do
	name=$(basename $file .bam)
	echo $name >> out/HISAT2_Metrics/Alignment-summary.txt
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SortSam \
	       I=$file \
	       O=tmp/.unmapped.sorted.bam \
	       SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SamToFastq \
	       I=tmp/.unmapped.sorted.bam \
	       F=tmp/.tmp.fastq \
	    | hisat2 -p 8 -x ${Index_PATH} \
		     -U tmp/.tmp.fastq -S /dev/stdout \
		     2>> out/HISAT2_Metrics/Alignment-summary.txt  \
	    | java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar SortSam \
		     I=/dev/stdin \
		     O=tmp/.mapped.sorted.sam \
		     SORT_ORDER=queryname;
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar MergeBamAlignment \
	       ATTRIBUTES_TO_RETAIN=XS \
	       UNMAPPED=tmp/.unmapped.sorted.bam  \
	       ALIGNED=tmp/.mapped.sorted.sam \
	       O=tmp/UMI/$name.umi.bam \
	       R=${Index_PATH}.fasta
    done
fi

rm tmp/.unmapped.sorted.bam
rm tmp/.mapped.sorted.sam
mkdir tmp/merged
mkdir tmp/Unaligned_bam
mv *.bam tmp/Unaligned_bam

# Merging all lanes that UMI-annotated BAM files corresponding to each sample derived from four lanes
for i in `seq 1 $NLINES`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar MergeSamFiles \
	 $(printf "I=%s " tmp/UMI/${OUTPUT_NAME}_${i}_Lane*.umi.bam) \
	 O=/dev/stdout |
	java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar AddOrReplaceReadGroups \
	     I=/dev/stdin \
	     O=tmp/merged/${OUTPUT_NAME}_${i}.merged.bam \
	     RGLB=${OUTPUT_NAME}_${i} RGPL=NextSeq RGPU=${i} RGSM=${i}
done

rm -rf tmp/UMI

# Mark potential PCR duplicates in each of the demultiplexed BAM files created from the merged BAM file, where each file represents a sample with a unique barcode
mkdir out/MarkDuplicates_Metrics
for i in `seq 1 $NLINES`
do
    java -Xmx5g -Djava.io.tmpdir=tmp -jar ${PICARDHOME_PATH}/picard.jar MarkDuplicates \
	 INPUT=tmp/merged/${OUTPUT_NAME}_${i}.merged.bam \
	 OUTPUT=out/${OUTPUT_NAME}_${i}.output.bam \
	 METRICS_FILE=out/MarkDuplicates_Metrics/${OUTPUT_NAME}_${i}.metrics.txt \
	 BARCODE_TAG=RX
done

rm -rf tmp/merged

# Preparation for annotation, extraction the information on genomic regions from annotation files and quality check
if [[ ${GENOME_VALUE} = "hg38" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq, KnownGenes, or Gencode for hg38"
    exit 1
elif [[ ${GENOME_VALUE} = "mm10" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for mm10"
    exit 1
elif [[ ${GENOME_VALUE} = "mm39" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "ens" ]]; then
    echo "No Ensembl gene annotations!! Please use RefSeq or KnownGenes, or Gencode for canFam6"
    exit 1   
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "No KnownGenes annotations!! Please use RefSeq or Ensembl for bosTau9"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam3" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for canFam3"
    exit 1
elif [[ ${GENOME_VALUE} = "canFam6" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for canFam6"
    exit 1
elif [[ ${GENOME_VALUE} = "mm9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq, KnownGenes, or Ensembl for mm9"
    exit 1
elif [[ ${GENOME_VALUE} = "bosTau9" ]] && [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "No Gencode annotations!! Please use RefSeq or Ensembl for bosTau9"
    exit 1
elif [[ ${ANNO_VALUE} =  "ens" ]]; then
    # Download gtfToGenePred binary (UCSC) and make it executable
    wget https://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/gtfToGenePred
    chmod +x gtfToGenePred

    # Download and decompress Ensembl GTF
    wget https://ftp.ensembl.org/pub/release-112/gtf/bos_taurus/Bos_taurus.ARS-UCD1.3.112.gtf.gz
    gunzip Bos_taurus.ARS-UCD1.3.112.gtf.gz

    # Generate ensGene.txt
    # Add an extra leading column (0) to match UCSC's 16-column format
    ./gtfToGenePred -genePredExt -ignoreGroupsWithoutExons \
    Bos_taurus.ARS-UCD1.3.112.gtf /dev/stdout \
    | awk 'BEGIN { OFS="\t" } { print 0, $0 }' \
    > src/ensGene.txt

    # Generate ensemblToGeneName.txt
    # Maps Ensembl transcript IDs to gene names
    ./gtfToGenePred -genePredExt -geneNameAsName2 -ignoreGroupsWithoutExons \
    Bos_taurus.ARS-UCD1.3.112.gtf /dev/stdout \
    | awk 'BEGIN { OFS="\t" } { print 0, $0 }' \
    > Bos_taurus.ARS-UCD1.3.112.gtf.refflat

    cut -f2,13 Bos_taurus.ARS-UCD1.3.112.gtf.refflat \
    | awk '$1 != "" && $2 != ""' \
    > src/ensemblToGeneName.txt

    # Cleanup intermediate files
    rm Bos_taurus.ARS-UCD1.3.112.gtf Bos_taurus.ARS-UCD1.3.112.gtf.refflat gtfToGenePred

    join -1 1 -2 2 -t $'\t' <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) > src/common.txt
    join -1 1 -2 2 -t $'\t' -v 2 <(sort -k 1,1 src/ensemblToGeneName.txt) <(sort -k 2,2 src/ensGene.txt) | awk 'BEGIN{OFS="\t"}{print $1,$13,$2,$1=$2="",$0}' | cut -f 1-3,7- > src/no-genename.txt
    cat src/common.txt src/no-genename.txt > src/ens-genes.txt
    rm src/common.txt && rm src/no-genename.txt
    ruby bin/ENSEMBL-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  "kg" ]]; then
    echo "Downloading the UCSC KnownGenes annotation data..."
    curl -o src/knownGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/knownGene.txt.gz
    curl -o src/kgXref.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/kgXref.txt.gz
    gunzip src/knownGene.txt.gz
    gunzip src/kgXref.txt.gz
    join  -1 1 -2 1 -t $'\t' <(sort -k 1,1 src/kgXref.txt | cut -f 1-5) <(sort -k 1,1 src/knownGene.txt) > src/knowngene-names.txt
    rm src/knownGene.txt && rm src/kgXref.txt
    ruby bin/KnownGenes-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  "ref" ]]; then
    echo "Downloading the NCBI RefSeq annotation data..."
    curl -o src/refGene.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/refGene.txt.gz
    gunzip src/refGene.txt.gz
    ruby bin/RefSeq-extract.rb
    shift 2
elif [[ ${ANNO_VALUE} =  wgEncodeGencodeBasic* ]]; then
    echo "Downloading the Gencode annotation data..."
    curl -o src/Gencode.txt.gz http://hgdownload.cse.ucsc.edu/goldenPath/${GENOME_VALUE}/database/${ANNO_VALUE}.txt.gz
    gunzip src/Gencode.txt.gz
    ruby bin/GENCODE-extract.rb
    shift 2
else
    echo "Something is wrong with the annotation data file."
    exit 1
fi

echo "Processing the chromosome size data..."


# Download and decompress the Bos taurus genome FASTA Ensembl v112
wget https://ftp.ensembl.org/pub/release-112/fasta/bos_taurus/dna/Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz
gunzip Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.gz

# Clean up chromosome headers: keep only 1–29, X, and MT (remove extra description)
sed '/^>/ s/>\([0-9XMT]\+\).*/>\1/' Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa \
    > cleaned_Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa

# Create FASTA index and extract chromosome sizes
samtools faidx cleaned_Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa
cut -f1,2 cleaned_Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.fai > bosTau9.chrom.sizes

# Place chrom.sizes in src/ and remove original uncleaned FASTA to save space
cp bosTau9.chrom.sizes src/
rm Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa cleaned_Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa.fai cleaned_Bos_taurus.ARS-UCD1.3.dna_sm.toplevel.fa

cat src/${GENOME_VALUE}.chrom.sizes | awk '{print $1"\t"1"\t"$2}' | sortBed -i > src/chrom.size.bed
cat src/proxup.bed | grep -v _alt  | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^NKLS | sortBed -i stdin | intersectBed -a stdin -b src/chrom.size.bed > src/proxup_trimmed.bed
cat src/5utr.bed src/proxup_trimmed.bed | grep -v _alt | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^NKLS | sortBed -i stdin | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - | grep -v , > src/coding_5end.bed
cat src/exon.bed src/proxup_trimmed.bed | grep -v _alt | grep -v _hap | grep -v _fix | grep -v _random | grep -v ^NKLS | sortBed -i stdin | mergeBed -s -o distinct,distinct,distinct -c 4,5,6 -i - > src/coding.bed
cat src/ERCC.bed src/coding_5end.bed | awk '{print $4 "\t" $1 "\t" $2+1 "\t" $3 "\t" $6}' > src/5end-regions.saf


rm src/5utr.bed
rm src/exon.bed
rm src/proxup.bed
rm src/proxup_trimmed.bed

# Quality check
cd out
echo -e Barcode"\t"Qualified_reads"\t"Total_reads"\t"Redundancy"\t"Mapped_reads"\t"Mapped_rate\
     "\t"Spikein_reads"\t"Spikein-5end_reads"\t"Spikein-5end_rate"\t"Coding_reads"\t"Coding-5end_reads"\t"Coding-5end_rate > ${OUTPUT_NAME}-QC.txt

for file in *.output.bam
do
    name=$(basename $file .output.bam)
    samtools index $file
    QR=$(samtools view -F 256 $file | wc -l)
    Total=$(samtools view -F 256 -F 1024 $file | wc -l)
    Redundancy=$(echo "scale=2;$QR/$Total" | bc)
    Map=$(samtools view -F 256 -F 1024 -F 4 $file | wc -l)
    Rate=$(echo "scale=1;$Map*100/$Total" | bc)
    Spike=$(samtools view -F 256 -F 1024 -F 4 $file |grep -e ERCC -e NIST| wc -l)
    spikein_5end_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/ERCC.bed | cut -f 4 | sort -u | wc -l)
    spikein_5end_rate=$(echo "scale=1;$spikein_5end_reads*100/$Spike" | bc)
    coding_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/coding.bed | cut -f 4 | sort -u | wc -l)
    coding_5end_reads=$(samtools view -u -F 256 -F 1024 -F 4 $file | intersectBed -abam stdin -wa -bed -b ../src/coding_5end.bed | cut -f 4 | sort -u | wc -l)
    coding_5end_rate=$(echo "scale=1;$coding_5end_reads*100/$coding_reads" | bc)
    echo -e $name"\t"$QR"\t"$Total"\t"$Redundancy"\t"$Map"\t"$Rate"\t"$Spike"\t"$spikein_5end_reads"\t"$spikein_5end_rate"\t"$coding_reads"\t"$coding_5end_reads"\t"$coding_5end_rate >> ${OUTPUT_NAME}-QC.txt
done

# Counting reads align to 5'end of genes with parameters by featureCounts
featureCounts -T 8 -s 1 --largestOverlap --ignoreDup --primary -a ../src/5end-regions.saf -F SAF -o ${OUTPUT_NAME}_byGene-counts.txt *.bam

mkdir Output_bai && mv *.bam.bai Output_bai
mkdir Output_bam && mv *.bam Output_bam

# Quality check plotting
# Run the R script
R CMD BATCH --slave --vanilla  ../bin/QC-plot.R QC-plot.R.log 
cd ..

# Removing reads that have flags 
for file in ${WorkingDir_PATH}/out/Output_bam/*.output.bam
do
    name=$(basename $file .output.bam)
    samtools view -u -F 256 -F 1024 -F 4 $file > ${name}.output.bam
    samtools index ${name}.output.bam
done

# Create spike-ins size file and chromosomal size file
wget https://tsapps.nist.gov/srmext/certificates/documents/SRM2374_putative_T7_products_NoPolyA_v2.FASTA
seqtk seq -L 5 SRM2374_putative_T7_products_NoPolyA_v2.FASTA > spikeins_size.fasta
samtools faidx spikeins_size.fasta
cut -f1-2 spikeins_size.fasta.fai > spikeins_sizes


cat spikeins_sizes >> ${GENOME_VALUE}.chrom.sizes
cp ${GENOME_VALUE}.chrom.sizes ${GENOME_VALUE}.chrom.sizes_with_spike_ins
rm SRM2374_putative_T7_products_NoPolyA_v2.FASTA
rm spikeins_*

cd out

# Run the R script to find scale factor for each sample
R CMD BATCH --slave --vanilla  ${WorkingDir_PATH}/bin/calculationScaleFactor.R scaleFactor.R.log 

cd ..

# Create bedGraph files and bigWig files for forward and reverse strand for each sample
for file in *.output.bam
do
    name=$(basename $file .output.bam)
    scale_factor=$(grep $file ${WorkingDir_PATH}/out/scale_factor_df.txt | awk '{print $2}')
    bedtools genomecov -5 -strand + -ibam $file -bg -scale $scale_factor | sort -k1,1 -k2,2n > ${WorkingDir_PATH}/${name}_plus.bedgraph
    bedtools genomecov -5 -strand - -ibam $file -bg -scale $scale_factor | sort -k1,1 -k2,2n > ${WorkingDir_PATH}/${name}_minus.bedgraph
    bedGraphToBigWig ${WorkingDir_PATH}/${name}_plus.bedgraph ${GENOME_VALUE}.chrom.sizes_with_spike_ins ${WorkingDir_PATH}/${name}_plus.bw
    bedGraphToBigWig ${WorkingDir_PATH}/${name}_minus.bedgraph ${GENOME_VALUE}.chrom.sizes_with_spike_ins ${WorkingDir_PATH}/${name}_minus.bw
done

# Create coding_5´end bigBed file
cp ${WorkingDir_PATH}/src/coding_5end.bed ${WorkingDir_PATH}/
sort -k1,1 -k2,2n coding_5end.bed > sorted_coding_5end.bed
wget http://hgdownload.soe.ucsc.edu/admin/exe/linux.x86_64/bedToBigBed
chmod a+x bedToBigBed

cat ${GENOME_VALUE}.chrom.sizes > ${GENOME_VALUE}.chrom.sizes_without_spike_ins
./bedToBigBed sorted_coding_5end.bed ${GENOME_VALUE}.chrom.sizes_without_spike_ins coding_5end.bb

rm coding_5end.bed
rm *.bedgraph
rm ${GENOME_VALUE}.chrom.sizes_with_spike_ins
rm out/scale_factor_df.txt
rm sorted_coding_5end.bed
rm ${GENOME_VALUE}.chrom.sizes_without_spike_ins
rm bedToBigBed
