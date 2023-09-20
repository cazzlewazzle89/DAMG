
mkdir FILTERED_FASTQ/

while read i j
do
    
    python /home/cwwalsh/Scripts/DAMG/ONT-16S/utils/filter_fastq.py \
	   --input_file "$j" \
	   --output_file FILTERED_FASTQ/"$i".fastq.gz \
	   --min_length 1400 \
	   --max_length 1700
    emu abundance \
	--type map-ont \
	--db /home/cwwalsh/Databases/Emu/ \
	--keep-counts \
	--output-dir EmuResults/ \
	--output-basename "$i" \
	--threads 10 \
	FILTERED_FASTQ/"$i".fastq.gz
	
done < manifest.txt

rm -f EmuResults/*_rel-abundance-threshold-0.0001.tsv

emu combine-outputs EmuResults/ species
emu combine-outputs --counts EmuResults/ species
