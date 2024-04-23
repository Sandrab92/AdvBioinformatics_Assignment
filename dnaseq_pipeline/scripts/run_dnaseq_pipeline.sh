#!/bin/bash

# 2.2. Pre-Alignment QC
# 1. Perform quality assessment on raw fastq files

# Define directory containing fastq.gz files
fastq_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/untrimmed_fastq/"

# Change directory to fastq directory
cd "$fastq_dir" || exit

# Run FastQC on all fastq.gz files
fastqc *.fastq.gz

echo "Quality assessment completed"



# 2.2. Pre-Alignment QC
# 1. Perform trimming on raw fastq files

# Define input directories and files
input_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/untrimmed_fastq/"
output_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/trimmed_fastq/"
forward_read="NGS0001_R1.fastq.gz"
reverse_read="NGS0001_R2.fastq.gz"

# Run Trimmomatic (including path to trimmomatic package and adapters)
trimmomatic PE -threads 4 -phred33 "$input_dir/$forward_read" "$input_dir/$reverse_read" -baseout "$output_dir/NGS0001_trimmed_R" ILLUMINACLIP:/home/ubuntu/anaconda3/pkgs/trimmomatic-0.39-hdfd78af_2/share/trimmomatic/adapters/NexteraPE-PE.fa:2:30:10 TRAILING:25 MINLEN:50

echo "Trimmomatic completed"



# 2.2. Pre-Alignment QC
# 2. Run fastqc on paired read from forward strand

# Define input directory
input_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/trimmed_fastq/"

# Run FastQC on paired trimmed sequencing data
# These can also be run on unpaired trimmed sequencing data
fastqc "${input_dir}NGS0001_trimmed_R_1P" "${input_dir}NGS0001_trimmed_R_2P"

echo "Quality assessment completed"



# 2.3. Alignment 
# Align the paired trimmed fastq files using bwa mem and reference genome hg19 (edit your bwa mem step to include read group information in your BAM file)

# 0.0. Index the reference genome

# Define path to reference genome
reference_genome="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/reference/hg19.fa.gz"

# Index the reference genome with BWA
bwa index "$reference_genome"

echo "Reference genome indexed successfully"


# 1.0. Alignment of trimmed data with BWA

# Define paths and filenames
reference_genome="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/reference/hg19.fa.gz"
input_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/trimmed_fastq/"
output_dir="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/"
forward_read="NGS0001_trimmed_R_1P"
reverse_read="NGS0001_trimmed_R_2P"
output_sam="NGS0001.sam"

# Run BWA mem with read group information
bwa mem -t 4 -R '@RG\tID:HWI-D0011.50.H7AP8ADXX.1.NGS0001\tSM:NGS0001\tPL:ILLUMINA\tLB:nextera-NGS0001-blood\tDT:2017-02-23\tPU:HWI-D00119' -I 250,50 "$reference_genome" "${input_dir}${forward_read}" "${input_dir}${reverse_read}" > "${output_dir}${output_sam}"

echo "BWA alignment completed. SAM file saved to ${output_dir}${output_sam}"


# 2.0. Sort sam file, convert sam to bam, generate bam index

# Define input and output paths
sam_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001.sam"
bam_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001.bam"
sorted_bam_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam"

# Convert SAM to BAM
samtools view -b -S -h -o "$bam_file" "$sam_file"

echo "Conversion completed"

# Sort BAM

# Define input and output paths
input_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001.bam"
output_sorted_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam"

# Run samtools sort

samtools sort -o "$output_sorted_bam" "$input_bam"

echo "Sorting completed. Output saved to $output_sorted_bam"


# Index sorted BAM

# Define input and output paths
input_sorted_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam"
output_bam_index="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam.bai"

# Index the sorted BAM file
samtools index "$input_sorted_bam" "$output_bam_index"

echo "Indexing completed. Output saved to $output_bam_index"



# 2.3. Post Alignment QC

# 1.0. Perfom Mark duplicates

# Define input and output paths
input_sorted_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted.bam"
output_marked_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam"
metrics_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_marked_dup_metrics.txt"

# Mark duplicates using Picard Tools
picard MarkDuplicates I="$input_sorted_bam" O="$output_marked_bam" M="$metrics_file"

echo "Duplicates marked. Output saved to $output_marked_bam"

# Index the marked BAM file
samtools index "$output_marked_bam"

echo "Indexing completed for marked BAM file."


# 2.0. Quality Filter the duplicate marked BAM file

# Define input and output paths
input_marked_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_marked.bam"
output_filtered_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam"

# Filter sorted marked BAM file based on MAPQ and bitwise flag
samtools view -F 1796 -q 20 -o "$output_filtered_bam" "$input_marked_bam"

echo "Filtered BAM file created: $output_filtered_bam"

# Index the sorted filtered BAM file
samtools index "$output_filtered_bam"

echo "Indexing completed for filtered BAM file."


# 3.0. Generate standard alignment statistics

# Define input BAM file path
input_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam"

# Define output files for flag statistics, index statistics, and depth
flagstats_output="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered_flagstats.txt"
idxstats_output="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered_idxstats.txt"
depth_output="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered_depth.txt"

# Generate flag statistics
samtools flagstat "$input_bam" > "$flagstats_output"
echo "Flag statistics generated: $flagstats_output"

# Generate index statistics
samtools idxstats "$input_bam" > "$idxstats_output"
echo "Index statistics generated: $idxstats_output"

# Generate coverage depth statistics
samtools depth "$input_bam" > "$depth_output"
echo "Depth statistics generated: $depth_output"


# Generate Insert size metrics and Histrogram

# Define input and output paths
input_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam"
output_metrics="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered_insert_size_metrics.txt"
output_histogram="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered_insert_size_histogram.pdf"
picard_jar="/home/ubuntu/anaconda3/pkgs/picard-2.18.29-0/share/picard-2.18.29-0/picard.jar"

# Collect insert size metrics
java -jar "$picard_jar" CollectInsertSizeMetrics I="$input_bam" O="$output_metrics" H="$output_histogram"




# 2.4. Variant Calling

# Define input and output paths
reference_genome="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/reference/hg19.fa"
indexed_reference="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/reference/hg19.fa.fai"
input_bam="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/aligned_data/NGS0001_sorted_filtered.bam"
output_vcf="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001.vcf"

# Decompress reference genome
if [ ! -f "$reference_genome" ]; then
    zcat "$reference_genome.gz" > "$reference_genome"
fi

# Index the reference genome
if [ ! -f "$indexed_reference" ]; then
    samtools faidx "$reference_genome"
fi

# Call variants using FreeBayes
freebayes --bam "$input_bam" --fasta-reference "$reference_genome" --vcf "$output_vcf"

# Check if Freebayes ran successfully
if [ $? -eq 0 ]; then
    echo "Variant calling completed. VCF file saved in $output_vcf"
    echo "Error: Freebayes failed to generate VCF file."
fi


# Compress and index VCF file

# Define input and output paths
vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001.vcf"
compressed_vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001.vcf.gz"

# Compress VCF file using bgzip
bgzip "$vcf_file"

# Index compressed VCF file using tabix
tabix -p vcf "$compressed_vcf_file"

# Check if compression and indexing were successful
if [ $? -eq 0 ]; then
    echo "Compression and indexing completed. Compressed VCF file saved in $compressed_vcf_file"
else
    echo "Error: Compression and indexing failed."
fi


# 2.0. Quality filter variants AND filter VCF based on regions specified in BED file

# Define input and output paths
vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001.vcf.gz"
filtered_vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered.vcf"
bed_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/data/annotation.bed"
filtered_regions_vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.vcf"

# Filter VCF file using vcffilter
vcffilter -f "QUAL > 1 & QUAL / AO > 10 & SAF > 0 & SAR > 0 & RPR > 1 & RPL > 1" "$vcf_file" > "$filtered_vcf_file"

# Filter VCF file based on regions specified in BED file using bedtools intersect
bedtools intersect -header -wa -a "$filtered_vcf_file" -b "$bed_file" > "$filtered_regions_vcf_file"

# Check if filtering was successful
if [ $? -eq 0 ]; then
    echo "VCF filtering based on quality and regions completed. Filtered VCF file saved in $filtered_regions_vcf_file"
else
    echo "Error: VCF filtering failed."
fi


# Compress and Index filtered regions

# Define input and output paths
filtered_regions_vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.vcf"
compressed_vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.vcf.gz"

# Compress filtered regions VCF file with bgzip
bgzip "$filtered_regions_vcf_file"

# Index compressed VCF file with tabix
tabix -p vcf "$compressed_vcf_file"

# Check if compression and indexing were successful
if [ $? -eq 0 ]; then
    echo "Compression and indexing of filtered regions VCF completed. Compressed and indexed VCF file saved in $compressed_vcf_file"
else
    echo "Error: Compression and indexing failed."
fi




# 2.5. Variant Annotation (Annovar)

# Define input and output paths
vcf_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.vcf.gz"
avinput_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.avinput"
convert2annovar="/home/ubuntu/annovar/convert2annovar.pl"

# Perform variant annotation using convert2annovar.pl
"$convert2annovar" -format vcf4 "$vcf_file" > "$avinput_file"

# Check if the conversion was successful
if [ $? -eq 0 ]; then
    echo "Variant annotation completed. Annotated output saved in $avinput_file"
else
    echo "Error: Variant annotation failed."
fi



# 2.5. Variant Prioritization (Annovar)

# Define input and output paths
avinput_file="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.avinput"
output_prefix="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions"
annovar_script="/home/ubuntu/annovar/table_annovar.pl"
humandb_dir="/home/ubuntu/annovar/humandb"
build_ver="hg19"

# Perform variant prioritization using table_annovar.pl
"$annovar_script" "$avinput_file" "$humandb_dir" -buildver "$build_ver" -out "$output_prefix" -remove -protocol refGene,ensGene,clinvar_20180603,exac03,dbnsfp31a_interpro -operation g,g,f,f,f -otherinfo -nastring .

# Check if the prioritization was successful
if [ -e "$output_prefix".hg19_multianno.txt ]; then
    echo "Variant prioritization completed. Annotated output saved in $output_prefix.hg19_multianno.txt"



# 2.5. Variant Annotation (SnpEff)

# Define input and output file paths
INPUT_FILE="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions.vcf.gz"
OUTPUT_FILE="/home/ubuntu/AdvBioinformatics_Assignment/dnaseq_pipeline/results/NGS0001_filtered_regions_annotated_snpeff.vcf"

# Define the command using path to updated java 11 to run snpEff latest version for variant annotation using GRCh38.86 genome build
COMMAND="/usr/lib/jvm/java-11-openjdk-amd64/bin/java -jar /home/ubuntu/snpEff/snpEff.jar eff -hgvs -canon -v GRCh38.86 $INPUT_FILE"

# Execute the command and redirect output to the specified file
$COMMAND > $OUTPUT_FILE

# Check if the command was successful
if [ $? -eq 0 ]; then
    echo "Annotation completed successfully."
else
    echo "Error: Annotation failed."
fi

# snpEff_genes.txt and snpEff_summary_files.html also generated after command is successfully executed

# Move snpEff_genes.txt and snpEff_summary.html to the results folder
mv snpEff_genes.txt snpEff_summary.html ~/AdvBioinformatics_Assignment/dnaseq_pipeline/results

# Check if the files were moved successfully
if [ $? -eq 0 ]; then
    echo "All output files successfully moved to results directory."
else
    echo "Error: Failed to move files."
fi

