#!/bin/sh

flag_manifest='manifest.tsv'
flag_outdir='Kneaded'
flag_trimmer='trimmomatic'
flag_host='human'
flag_hostfiltering='bt2_strict'
flag_threads=10
inputerror='false'

echo ''

while getopts ':hi:o:a:c:f:t:' opt
do
    case ${opt} in
	h)
	    echo -e '\tPerforms quality trimming and host DNA removal on Illumina-sequenced metagenomic data'
        echo -e '\tRun in kneaddata conda env at /home/cwwalsh/miniconda3/envs/kneaddata/'
        echo -e ''
        echo -e '\tOptions:'
        echo -e ''
	    echo -e '\t -i Manifest File (headerless TSV file with sample name and absolute paths to demultiplexed forward and reverse reads) [default: manifest.tsv]'
		echo -e '\t -o Output Directory [default: Kneaded/]'
        echo -e '\t -a Trimming Algorithm [default: trimmomatic]'
            echo -e '\t\t Currently Implemented: trimmomatic,fastp'
	    echo -e '\t -c Host Name For Decomtamination [default: human]'
            echo -e '\t\t Currently Accepted Values: mouse,human'
        echo -e '\t -f Host DNA Filtering Method [default: bt2_strict]'
            echo -e '\t\t bt2_strict : Will remove read pairs if Bowtie2 aligns EITHER read host genome'
            echo -e '\t\t bt2_lenient : Will remove read pairs if Bowtie2 aligns BOTH reads to host genome'
            echo -e '\t\t srascrubber : Use ncbi::sra-human-scrubber instead of Bowtie2 - currently only works with human samples'
		echo -e '\t -t Threads/CPUs To Use [default: 10]'
	    echo ''
	    exit 0
	    ;;
	i) flag_manifest=$OPTARG ;;
	o) flag_outdir=$OPTARG ;;
    a) flag_trimmer=$OPTARG ;;
	c) flag_host=$OPTARG ;;
	f) flag_hostfiltering=$OPTARG ;;
	t) flag_threads=$OPTARG ;;
	\?) echo -e '\t Usage: knead_for_speed.sh -i ManifestFile \n\tOR\n\tHelp and Optional Arguments: knead_for_speed.sh -h\n' >&2
	    exit 1
	    ;;
	:) echo -e '\t Error: Use -h for full options list\n'
	   exit 1
    esac
done

# CONFIRM THAT MANIFEST FILE EXISTS
if [ ! -f $flag_manifest ]
then
	echo 'Manifest File Not Found'
    inputerror='true'
fi

# CONFIRM THAT FILES LISTED IN MANIFEST EXIST
while read sample read1 read2
do
	if [ ! -f $read1 ]
	then
		echo 'File '$read1' does not exist'
		inputerror='true'
	fi

	if [ ! -f $read2 ]
	then
		echo 'File '$read2' does not exist'
		inputerror='true'
	fi
done < $flag_manifest

# CONFIRM THAT SUPPORTED TRIMMER IS SPECIFIED
if [ $flag_trimmer == 'trimmomatic' ] || [ $flag_trimmer == 'fastp' ]
then
    :
else
    echo 'Trimmer name not recognised'
    inputerror='true'
fi

# CONFIRM THAT SUPPORTED HOST IS SPECIFIED
if [ $flag_host == 'mouse' ]
then
    bt2db='/home/cwwalsh/Databases/Kneaddata_HostRemoval/GRCm39'
elif [ $flag_host == 'human' ]
then
    bt2db='/home/cwwalsh/Databases/Kneaddata_HostRemoval/hg37dec_v0.1'
else
    echo 'Host database name not recognised'
    inputerror='true'
fi

# CONFIRM THAT SUPPPORTED HOST REMOVAL METHOD IS SPECIFIED
if [ $flag_hostfiltering == 'bt2_strict' ] || [ $flag_hostfiltering == 'bt2_lenient' ] || [ $flag_hostfiltering == 'srascrubber' ]
then
    :
else
    echo 'Host filtering method not recognised'
    inputerror='true'
fi

# CONFIRM THAT SRASCRUBBER IS USED ONLY ON HUMAN SAMPLES
if [ $flag_hostfiltering == 'srascrubber' ] && [ $flag_host == 'mouse' ]
then
    echo 'Sra-scrubber currently incompatible with mouse, only human'
    echo 'Please use Bowtie2 or hassle Calum until he builds a mouse database'
    inputerror='true'
fi

# IF ANY OF THE ABOVE CONDITIONS ARE VIOLATED, PRINT MESSAGES AND EXIT
if [ $inputerror = 'true' ]
then
    echo ''
    exit 1
fi

# IF OUTPUT DIRECTORY ALREADY EXISTS
# WARN THAT CONTENTS WILL BE OVERWRITTEN
# WAIT 5 SECONDS TO GIVE USER TIME TO CANCEL
# CONTINUE
if [ ! -d $flag_outdir ]
then
	mkdir -p $flag_outdir
else
    echo ''
	echo 'Output Directory '$flag_outdir' Already Exists: Contents Will Be Overwritten'
    sleep 5
fi

# MAIN LOOP
while read sample read1 read2
do
    if [ $flag_trimmer == 'trimmomatic' ]
    then
        trimmomatic PE \
            "$read1" "$read2" \
            "$flag_outdir"/"$sample"_R1_paired.fastq.gz \
            "$flag_outdir"/"$sample"_R1_unpaired.fastq.gz \
            "$flag_outdir"/"$sample"_R2_paired.fastq.gz \
            "$flag_outdir"/"$sample"_R2_unpaired.fastq.gz \
            MINLEN:60 ILLUMINACLIP:/home/cwwalsh/miniconda3/envs/kneaddata/lib/python3.9/site-packages/kneaddata/adapters/NexteraPE-PE.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:50 \
            -threads $flag_threads

        rm -f "$flag_outdir"/"$sample"_R1_unpaired.fastq.gz "$flag_outdir"/"$sample"_R2_unpaired.fastq.gz
    else
        fastp \
            --in1 "$read1" \
            --in2 "$read2" \
            --out1 "$flag_outdir"/"$sample"_R1_paired.fastq.gz \
            --out2 "$flag_outdir"/"$sample"_R2_paired.fastq.gz \
            --detect_adapter_for_pe \
            --length_required 50 \
            --thread $flag_threads \
            --html "$flag_outdir"/"$sample"_fastp.html \
            --json "$flag_outdir"/"$sample"_fastp.json
    fi

    if [ $flag_hostfiltering == 'srascrubber' ]
    then
        # SRASCRUBBER REQUIREs UNCOMPRESSED INTERLEAVEs FASTQ INPUTS
        # INTERLEAVING INPUT, AND DEINTERLEAVING OUTPUT, WITH BBTOOLS
        reformat.sh \
            in1="$flag_outdir"/"$sample"_R1_paired.fastq.gz \
            in2="$flag_outdir"/"$sample"_R2_paired.fastq.gz \
            out=stdout.fq | scrub.sh \
                -p $flag_threads \
                -x | reformat.sh \
                    in=stdin.fq \
                    int=t \
                    out1="$flag_outdir"/"$sample"_R1.fastq.gz \
                    out2="$flag_outdir"/"$sample"_R2.fastq.gz

        rm -f "$flag_outdir"/"$sample"_R1_paired.fastq.gz
        rm -f "$flag_outdir"/"$sample"_R2_paired.fastq.gz

    else
        bowtie2 \
            --threads $flag_threads \
            --seed 42 \
            -x "$bt2db" \
            -1 "$flag_outdir"/"$sample"_R1_paired.fastq.gz \
            -2 "$flag_outdir"/"$sample"_R2_paired.fastq.gz | samtools view -b > "$flag_outdir"/"$sample".bam 

        # BAM FLAG FILTERING TAKEN FROM HERE ( https://gist.github.com/darencard/72ddd9e6c08aaff5ff64ca512a04a6dd )
        if [ $flag_hostfiltering == 'bt2_strict' ]
        then
            # 'strict' - DISCARD READ PAIRS WHERE AT LEAST ONE READ MAPS TO HOST GENOME
            samtools view -f 12 -F 256 "$flag_outdir"/"$sample".bam > "$flag_outdir"/"$sample"_microbial.bam
        else
            # 'lenient' - ONY DISCARD READ PAIRS IF BOTH READS MAP TO HOST GENOME
            # R1 & R2 UNMAPPED
            samtools view -u -f 12 -F 256 "$flag_outdir"/"$sample".bam > "$flag_outdir"/"$sample"_unmapped.bam
            # R1 MAPPED, R2 UNMAPPED
            samtools view -u -f 8 -F 260 "$flag_outdir"/"$sample".bam > "$flag_outdir"/"$sample"_R1mapped.bam
            # R1 UNMAPPED, R2 MAPPED
            samtools view -u -f 4 -F 264 "$flag_outdir"/"$sample".bam > "$flag_outdir"/"$sample"_R2mapped.bam
            # COMBINE THESE
            samtools merge \
                -o "$flag_outdir"/"$sample"_microbial.bam \
                "$flag_outdir"/"$sample"_unmapped.bam \
                "$flag_outdir"/"$sample"_R1mapped.bam \
                "$flag_outdir"/"$sample"_R2mapped.bam
        fi

        # CONVERT UNMAPPED READ PAIRS TO FASTQ FORMAT
        samtools fastq \
            -1 "$flag_outdir"/"$sample"_R1.fastq.gz \
            -2 "$flag_outdir"/"$sample"_R2.fastq.gz \
            -0 /dev/null \
            -s /dev/null \
            -n \
            "$flag_outdir"/"$sample"_microbial.bam

        rm -f "$flag_outdir"/"$sample".bam
        rm -f "$flag_outdir"/"$sample"_R1mapped.bam
        rm -f "$flag_outdir"/"$sample"_R2mapped.bam
        rm -f "$flag_outdir"/"$sample"_R1_paired.fastq.gz 
        rm -f "$flag_outdir"/"$sample"_R2_paired.fastq.gz 
        rm -f "$flag_outdir"/"$sample"_microbial.bam
    fi

done < $flag_manifest

# CREATE NEW MANIFEST FILE WITH SAMPLE_ID AND ABSOLUTE FILEPATHS TO POST-QC FORWARD AND REVERSE READS
while read sample read1 read2
do
    echo $sample "$PWD"/"$flag_outdir"/"$sample"_R1.fastq.gz "$PWD"/"$flag_outdir"/"$sample"_R2.fastq.gz
done < $flag_manifest > manifest_kneaded.tsv
