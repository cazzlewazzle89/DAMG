#!/bin/sh

i_flag=''
m_flag=''
o_flag=${PWD}
l_flag='Unspecified'
s_flag='false'
r_flag='Unspecified'
d_flag='false'
p_flag='skip'
a_flag='skip'
t_flag='8'
inputerror='false'

echo ''

while getopts ':hi:m:o:l:sr:d:p:a:t:' opt
do
    case ${opt} in
	h)
	    echo -e 'Options:'
	    echo -e '\t -i Manifest File (TSV file with sample name and absolute paths to demultiplexed forward and reverse reads) [Mandatory]'
	    echo -e '\t -m Metadata File [Mandatory]'
		echo -e '\t -o Output Directory [Optional]'
			echo -e '\t\t If Unspecified: Will Output Files To Current Working Directory'
	    echo -e '\t -l Trim Length [Optional]'
			echo -e '\t\t If Unspecified: Pipeline Will Pause And Wait for Prompt'
			echo -e '\t\t This Is To Allow User To Evaluate Read Quality Metrics In demux.qzv'
	    echo -e '\t -s Deblur Only Run Mode [Optional] [Default: False]'
            echo -e '\t\t Specifying This Option Will Skip Taxonomic Classification And Phylogenetic Placement Of sOTUs, Calculation Of Alpha And Beta Diversity, And ANCOM Analysis'
            echo -e '\t\t Useful When Samples Are Split Across Multiple Sequencing Runs'
	    echo -e '\t -r Rarefaction Depth [Optional]'
		echo -e '\t -d Run Decontam [Optional] [Default: False]'
			echo -e '\t\t If Desired: Specify Name Of Column In Metadata Sheet Describing Whether Each Sample Is A Negative Control'
			echo -e '\t\t Column Must Be Boolean With Values TRUE Or FALSE'
	    echo -e '\t -p Metadata Variable(s) For PERMANOVA [Optional] [Default: Skip]'
	    	echo -e '\t\t For Multiple Variables Use A Space Separated List Within Double Quotes'
	    	echo -e '\t\t eg. "meta1 meta2 meta3"'
	    echo -e '\t -a Metadata Variable(s) For ANCOM [Optional] [Default: Skip]'
	    	echo -e '\t\t eg. "meta1 meta2 meta3"'
	    	echo -e '\t\t For Multiple Variables Use A Space Separated List Within Double Quotes'
		echo -e '\t -t Threads Or Parallel Jobs To Use When Possible [Optional] [Default: 8]'
	    echo ''
	    exit 0
	    ;;
	i) i_flag=$OPTARG ;;
	m) m_flag=$OPTARG ;;
	o) o_flag=$OPTARG ;;
	l) l_flag=$OPTARG ;;
	r) r_flag=$OPTARG ;;
	d) d_flag=$OPTARG ;;
	p) p_flag=$OPTARG ;;
	a) a_flag=$OPTARG ;;
	s) s_flag='true' ;;
	t) t_flag=$OPTARG ;;
	\?) echo -e '\t Usage: damg_16S_manifestformat.sh -i PathToManifestFile -m MetadataFile\n\tOR\n\tHelp and Optional Arguments: damg_16S_manifestformat.sh -h\n' >&2
	    exit 1
	    ;;
	:) echo -e '\t Error: Use -h for full options list\n'
	   exit 1
    esac
done

if [ "$i_flag" = '' ] || [ "$m_flag" = '' ]
then
    echo 'ERROR'
    echo 'Usage: damg_16S_manifestformat.sh -i PathToManifestFile -m MetadataFile'
    echo 'Use -h for full options list'
    echo ''
    exit 1
else
    echo 'Manifest file: '$i_flag
    echo 'Metadata File: '$m_flag
	echo 'Output Directory: '$o_flag
	echo 'Threads or Parallel Jobs To Use: '$t_flag
fi

if [ $l_flag = 'Unspecified' ]
then
    echo 'Trim Length: Unspecified - User Input Will Be Required'
else
    echo 'Trim Length: '$l_flag'bp'
fi

if [ $d_flag = 'false' ]
then
    echo 'Skipping Decontam'
else
    echo 'Decontam metadata column: '$d_flag
fi

if [ $s_flag = 'false' ]
then
    echo 'Deblur Only Run Mode: False (Running Full Pipeline)'

    if [ $r_flag = 'Unspecified' ]
    then
		echo -e 'Rarefaction Depth: Unspecified - User Input Will Be Required'
    else
		echo -e 'Rarefaction Depth: '$r_flag
    fi

    if [[ $p_flag = 'skip' ]]
    then
		echo -e 'Metadata Variable(s) for PERMANOVA: Unspecified - Skipping'
    else
		echo -e 'Metadata Variable(s) for PERMANOVA: '$p_flag
    fi

    if [[ $a_flag = 'skip' ]]
    then
		echo -e 'Metadata Variable(s) for ANCOM: Unspecified - Skipping'
    else
		echo -e 'Metadata Variable(s) for ANCOM: '$a_flag
    fi
elif [ $s_flag = 'true' ]
then
    echo 'Deblur Only Run Mode: True (Pipeline Will Finish After sOTU generation)'
else
    echo 'Deblur Only Run Mode: Incorrectly Specified - must be either "false" or "true"'
    inputerror='true'
fi

if [ ! -f $i_flag ] || [ ! -f $b_flag ] || [ ! -f $m_flag ]
then
    echo ''
    
    if [ ! -f $i_flag ]
    then
		echo 'Manifest File Not Found'
		inputerror='true'
    fi

    if [ ! -f $b_flag ]
    then
		echo 'Barcodes File Not Found'
		inputerror='true'
    fi

    if [ ! -f $m_flag ]
    then
		echo 'Metadata File Not Found'
		inputerror='true'
    fi
fi

if [ ! -d $o_flag ]
then
	mkdir -p $o_flag
else
	echo 'Output Directory '$o_flag' Already Exists: Contents Will Be Overwritten'
fi

csvtk headers -t $i_flag > temp_manifestheaders

if $(sed -n '1p' temp_manifestheaders | grep -x -q 'sample-id')
then
	:
else
	echo 'First column in Manifest File should be named sample-id'
	inputerror='true'
fi

if $(sed -n '2p' temp_manifestheaders | grep -x -q 'forward-absolute-filepath')
then
	:
else
	echo 'Second column in Manifest File should be named forward-absolute-filepath'
	inputerror='true'
fi

if $(sed -n '3p' temp_manifestheaders | grep -x -q 'reverse-absolute-filepath')
then
	:
else
	echo 'Third column in Manifest File should be named reverse-absolute-filepath'
	inputerror='true'
fi

rm -f temp_manifestheaders

sed '1d' $i_flag > manifest.temp

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
done < manifest.temp

rm -f manifest.temp
  
echo ''

if [ $inputerror = 'true' ]
then
    exit 1
fi

qiime tools import \
  --type 'SampleData[PairedEndSequencesWithQuality]' \
  --input-path $i_flag \
  --output-path demux.qza \
  --input-format PairedEndFastqManifestPhred33V2

qiime demux summarize \
	--i-data demux.qza \
	--o-visualization demux.qzv

qiime quality-filter q-score \
	--i-demux demux.qza \
	--o-filtered-sequences demux-filtered.qza \
	--o-filter-stats demux-filterstats.qza

qiime metadata tabulate \
	--m-input-file demux-filterstats.qza \
	--o-visualization demux-filterstats.qzv

if [ $l_flag = 'Unspecified' ]
then
    read -p 'Trim Length: ' trimlength_var
    qiime deblur denoise-16S \
		--i-demultiplexed-seqs demux-filtered.qza \
		--p-trim-length $trimlength_var \
		--p-sample-stats \
		--p-jobs-to-start $t_flag \
		--p-no-hashed-feature-ids \
		--o-table feature-table.qza \
		--o-representative-sequences rep-seqs.qza \
		--o-stats deblur-stats.qza
else
    qiime deblur denoise-16S \
		--i-demultiplexed-seqs demux-filtered.qza \
		--p-trim-length $l_flag \
		--p-sample-stats \
		--p-jobs-to-start $t_flag \
		--p-no-hashed-feature-ids \
		--o-table feature-table.qza \
		--o-representative-sequences rep-seqs.qza \
		--o-stats deblur-stats.qza
fi

if [ $d_flag != 'false' ]
then
	mv feature-table.qza feature-table-predecontam.qza

	Rscript /home/cwwalsh/Scripts/DAMG/decontam.R $m_flag $d_flag

	qiime tools import \
		--input-path feature-table-decontam.biom \
		--type 'FeatureTable[Frequency]' \
		--input-format BIOMV100Format \
		--output-path feature-table.qza

	rm feature-table-decontam.biom
fi

if [ $s_flag = 'true' ]
then
	if [ ! $o_flag = $PWD ]
	then
		mv deblur.log deblur-stats.qza $o_flag
		mv demux-filtered.qza demux-filterstats.qza demux-filterstats.qzv demux.qza demux.qzv $o_flag
		mv feature-table.qza rep-seqs.qza $o_flag
		rm -f raw-seqs.qza

		if [ $d_flag != 'false' ]
		then
			mv contaminant_otus.csv feature-table-predecontam.qza $o_fla
		fi
	fi

    echo 'Skipping Analysis Steps As Instructed'
    echo 'Pipeline Complete!'
    echo ''
    exit 0
fi

qiime feature-table summarize \
	--i-table feature-table.qza \
	--m-sample-metadata-file $m_flag \
	--o-visualization feature-table.qzv

qiime feature-table tabulate-seqs \
	--i-data rep-seqs.qza \
	--o-visualization rep-seqs.qzv

qiime fragment-insertion sepp \
	--i-representative-sequences rep-seqs.qza \
	--i-reference-database /home/cwwalsh/Scripts/DAMG/sepp-refs-gg-13-8.qza \
	--o-tree insertion-tree.qza \
	--o-placements insertion-placements.qza

qiime feature-classifier classify-sklearn \
	--i-reads rep-seqs.qza \
	--i-classifier /home/cwwalsh/Scripts/DAMG/gg-13-8-99-515-806-nb-classifier.qza \
	--o-classification taxonomy.qza

qiime metadata tabulate \
	--m-input-file taxonomy.qza \
	--o-visualization taxonomy.qzv

errcode=1
while [ $errcode -ne 0 ]
do
    if [ $r_flag = 'Unspecified' ]
    then
	
	read -p 'Rarefaction depth(s) for testing (can provide multiple depths separated by a space): ' rarefactiondepths_var
	
	for depth in $rarefactiondepths_var
	do
	    qiime diversity alpha-rarefaction \
			--i-table feature-table.qza \
			--i-phylogeny insertion-tree.qza \
			--p-max-depth $depth \
			--m-metadata-file $m_flag \
			--o-visualization alpha-rarefaction-"$depth".qzv
	done
	
	errcode=$?
	
	if [ $errcode -ne 0 ]
	then
	    echo 'Inappropriate Rarefaction Depth Specified'
	    echo 'Please Retry'
	    errcode=1
	else
	    echo 'Rarefaction testing complete!'
	    read -p 'Rarefaction depth for diversity analysis (entering no value will skip this step): ' diversityrarefactiondepth_var
	    if [ ! $diversityrarefactiondepth_var = '' ]
	    then
		qiime diversity core-metrics-phylogenetic \
			--i-table feature-table.qza \
			--i-phylogeny insertion-tree.qza \
			--p-sampling-depth $diversityrarefactiondepth_var \
			--m-metadata-file $m_flag \
			--p-n-jobs-or-threads $t_flag \
			--output-dir CoreMetricsPhylogenetic
	    fi
	fi
    else
	qiime diversity alpha-rarefaction \
		--i-table feature-table.qza \
		--i-phylogeny insertion-tree.qza \
		--p-max-depth $r_flag \
		--m-metadata-file $m_flag \
		--o-visualization alpha-rarefaction-"$r_flag".qzv
	qiime diversity core-metrics-phylogenetic \
		--i-table feature-table.qza \
		--i-phylogeny insertion-tree.qza \
		--p-sampling-depth $r_flag \
		--m-metadata-file $m_flag \
		--p-n-jobs-or-threads $t_flag \
		--output-dir CoreMetricsPhylogenetic
	errcode=$?
    fi
done

for i in $(ls CoreMetricsPhylogenetic/*_vector.qza | sed 's/_vector\.qza//')
do
    qiime diversity alpha-group-significance \
		--i-alpha-diversity "$i"_vector.qza \
		--m-metadata-file $m_flag \
		--o-visualization "$i"_groupsig.qzv
done

if [[ $p_flag = 'skip' ]]
then
	echo 'Skipping PERMANOVA As Instructed'
else
    for metric in $(ls CoreMetricsPhylogenetic/*_distance_matrix.qza | sed 's/_distance_matrix.qza//')
    do
		for variable in $p_flag
		do
	    	qiime diversity beta-group-significance \
				--i-distance-matrix "$metric"_distance_matrix.qza \
				--m-metadata-file $m_flag \
				--m-metadata-column $variable \
				--p-pairwise \
				--o-visualization "$metric"_groupsig_"$variable".qzv
		done
    done
fi

qiime taxa barplot \
	--i-table CoreMetricsPhylogenetic/rarefied_table.qza \
	--i-taxonomy taxonomy.qza \
	--m-metadata-file $m_flag \
	--o-visualization taxa-bar-plots.qzv

qiime composition add-pseudocount \
	--i-table feature-table.qza \
	--o-composition-table comp-feature-table.qza

if [[ $a_flag = 'skip' ]]
then
    echo 'Skipping ANCOM As Instructed'
else
    for variable in $a_flag
    do
		qiime composition ancom \
			--i-table comp-feature-table.qza \
			--m-metadata-file $m_flag \
			--m-metadata-column $variable \
			--o-visualization ancom-"$variable".qzv
    done
fi

for metric in $(ls CoreMetricsPhylogenetic/*_distance_matrix.qza | sed 's/_distance_matrix.qza// ; s/CoreMetricsPhylogenetic\///')
do
    qiime tools export \
        --input-path CoreMetricsPhylogenetic/"$metric"_distance_matrix.qza \
        --output-path CoreMetricsPhylogenetic/ \
    
    mv CoreMetricsPhylogenetic/distance-matrix.tsv CoreMetricsPhylogenetic/"$metric"_distance_matrix.tsv
done

qiime tools export \
    --input-path taxonomy.qza \
    --output-path .

qiime tools export \
    --input-path insertion-tree.qza \
    --output-path .

sed -i 's/; /| /' tree.nwk

qiime tools export \
    --input-path feature-table.qza \
    --output-path .

biom convert \
    -i feature-table.biom \
    -o feature-table_json.biom \
    --table-type="OTU table" \
    --to-json

if [ ! $o_flag = $PWD ]
then
	mv deblur.log deblur-stats.qza $o_flag
	mv demux-filtered.qza demux-filterstats.qza demux-filterstats.qzv demux.qza demux.qzv $o_flag
	mv feature-table.qza feature-table.qzv rep-seqs.qza $o_flag
	mv rep-seqs.qzv insertion-tree.qza insertion-placements.qza taxonomy.qza taxonomy.qzv $o_flag
	mv alpha-rarefaction-*.qzv CoreMetricsPhylogenetic/ taxa-bar-plots.qzv comp-feature-table.qza $o_flag

	if [ $d_flag != 'false' ]
		then
			mv contaminant_otus.csv feature-table-predecontam.qza $o_fla
		fi
	
	if [ ! $a_flag = '' ] && [ ! $ancomvariables_var = '' ]
	then
		mv ancom-*.qzv $o_flag
	fi

	mv taxonomy.tsv tree.nwk feature-table.biom feature-table_json.biom $o_flag
	
	rm -f raw-seqs.qza 
fi
	
echo 'Pipeline Complete!'
echo ''
