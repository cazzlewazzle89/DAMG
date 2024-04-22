#!/bin/sh

# USAGE: sh Scripts/ont_backup_stats_share.sh runID renamingTSV

set -e
RUNID=$1
RENAMINGTSV=$2

INPUTDIRECTORY=/home/damg/for_transfer/"$RUNID"
OUTPUTDIRECTORY=/home/damg/data/"$RUNID"

# CONFIRM THAT INPUT DIR EXISTS
if [ ! -d "$INPUTDIRECTORY" ]
then
    echo "Input Directory Does Not Exist"
    exit 1
else
    echo "Input Directory:" "$INPUTDIRECTORY"
fi

# CONFIRM THAT OUTPUT DIR DOES NOT EXIST
if [ -d "$OUTPUTDIRECTORY" ]
then
    echo "Output Directory Already Exists"
    exit 1
fi

# IF WRITE PERMISSIONS ARE SET CORRECTNLY ON INPUT DIR
# MAKE OUTPUT DIR
# EMPTY OUTPUT FASTQ_PASS DIR
# MOVE POD5, REPORT, AND ANALYSIS DIRS TO OUTPUT DIR
if [ -w "$INPUTDIRECTORY" ]
then
    echo "Ouput Directory:" "$OUTPUTDIRECTORY"
    mkdir -p "$OUTPUTDIRECTORY"/fastq_pass/
    
    if [ -d "$INPUTDIRECTORY"/fast5_pass ]
    then
        mv "$INPUTDIRECTORY"/fast5_pass/ "$OUTPUTDIRECTORY"
    fi

    if [ -d "$INPUTDIRECTORY"/pod5_pass ]
    then
        mv "$INPUTDIRECTORY"/pod5_pass/ "$OUTPUTDIRECTORY"
    fi

    if [ -d "$INPUTDIRECTORY"/fast5 ]
    then
        mv "$INPUTDIRECTORY"/fast5/ "$OUTPUTDIRECTORY"
    fi
    
    if [ -d "$INPUTDIRECTORY"/pod5 ]
    then
        mv "$INPUTDIRECTORY"/pod5/ "$OUTPUTDIRECTORY"
    fi

    if [ -d "$INPUTDIRECTORY"/reports ]
    then
        mv "$INPUTDIRECTORY"/reports/ "$OUTPUTDIRECTORY"
    fi

    if [ -d "$INPUTDIRECTORY"/other_reports ]
    then
        mv "$INPUTDIRECTORY"/other_reports/ "$OUTPUTDIRECTORY"
    fi

    if [ -d "$INPUTDIRECTORY"/analysis ]
    then
        mv "$INPUTDIRECTORY"/analysis/ "$OUTPUTDIRECTORY"
    fi
else
    echo 'Permission Denied'
    exit 1
fi

# IF INPUT FASTQ_PASS DIR EXISTS
# PRINT SOME MESSAGES TO SCREEN
# CONCATENATE MULTIPLE FASTQ FILES PER SAMPLES INTO A SINGLE ONE
if [ -d "$INPUTDIRECTORY"/fastq_pass ]
then
    SAMPLECOUNT=$(ls -d "$INPUTDIRECTORY"/fastq_pass/* | wc -l)
    echo "$SAMPLECOUNT" "Samples Found"

    for i in $(ls -d "$INPUTDIRECTORY"/fastq_pass/* | sed 's/.*fastq_pass\///')
    do
        echo "Concatenating Reads From" "$i"
        pigz -cd -p 10 "$INPUTDIRECTORY"/fastq_pass/"$i"/*.fastq.gz | pigz -p 10 > "$OUTPUTDIRECTORY"/fastq_pass/"$i".fastq.gz
    done
fi

# SAME AS ABOVE, BUT IF INPUT FASTQ DIR IS CALLED FASTQ/
if [ -d "$INPUTDIRECTORY"/fastq ]
then
    SAMPLECOUNT=$(ls -d "$INPUTDIRECTORY"/fastq/* | wc -l)
    echo "$SAMPLECOUNT" "Samples Found"

    for i in $(ls -d "$INPUTDIRECTORY"/fastq/* | sed 's/.*fastq\///')
    do
        echo "Concatenating Reads From Sample" "$j" "(""$i"")"
        pigz -cd -p 10 "$INPUTDIRECTORY"/fastq/"$i"/*.fastq.gz | pigz -p 10 > "$OUTPUTDIRECTORY"/fastq/"$i".fastq.gz
    done
fi

# MOVE RENAMING INFO FILE TO OUTPUT DIR IF PROVIDED
if [ -z "$2" ]
then
    echo "Barcode-To-Name File Not Provided. Will Not Rename Samples"
else
    if [ ! -f "$RENAMINGTSV" ]
    then
        echo "Barcode-To-Name File Does Not Exist"
        exit 1
    else
        echo "Barcode-To-Name File:" "$RENAMINGTSV"
        cp "$RENAMINGTSV" "$OUTPUTDIRECTORY"/renaming.tsv
    fi
fi

# MOVE TO OURPUT FASTQ_PASS DIR FOR RENAMING STEP
if [ -d "$INPUTDIRECTORY"/fastq_pass ]
then
    cd "$OUTPUTDIRECTORY"/fastq_pass/
fi

# MOVE TO OURPUT FASTQ DIR FOR RENAMING STEP
if [ -d "$INPUTDIRECTORY"/fastq ]
then
    cd "$OUTPUTDIRECTORY"/fastq/
fi

# RENAME SAMPLES
if [ -n "$2" ]
then
    echo "Renaming Samples"
    while read i j
    do
        mv "$i".fastq.gz "$j".fastq.gz
    done < ../renaming.tsv
fi

# CREATE SEQ STATS
echo "Generating Read Length and Quality Statistics"
seqkit stats -a * > ../seqkit_stats.tsv

# SHARE ALL FILES FROM RUN
ls *.fastq.gz | sed 's/.fastq.gz$//' > names.tsv
mdu share --input_file names.tsv --source . -t /home/cwwalsh/share_template.txt
rm -f names.tsv

# PRINT SEQ STATS TO SCREEN
cat ../seqkit_stats.tsv


