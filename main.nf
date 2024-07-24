#! /usr/bin/env nextflow

log.info """
	Chrombpnet Pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/*
 * merge sorted and indexed input bam files by sample
 */
 
 process merge_bam {
	tag "$meta.sample"
	conda "bioconda::samtools=1.20"
	container = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
	
	input:
	tuple val(meta), path(bam), path(bam_index)
	
	output:
	tuple val(meta), 
	path("${meta.sample}_merged.bam"), 
	path("${meta.sample}_merged.bam.bai")
	
	script:
    """
    samtools merge -o ${meta.sample}_merged.bam --threads $task.cpus $bam
    samtools index ${meta.sample}_merged.bam
    """
}
 
 
/*
 * Call peaks on merged samples
 */
 
  process call_peaks {
	tag "$meta.sample"
	conda "bioconda::macs3=3.01"
	container = "quay.io/biocontainers/macs3:3.0.1--py312he57d009_3"
	
	input:
	tuple val(meta), path(bam), path(bam_index)
	
	output:
	tuple val(meta), 
	path("${meta.sample}_peaks.narrowPeak")
	
	script:
    """
    macs3 callpeak \
    -t $bam \
    -n $meta.sample \
    -f BAMPE \
    --keep-dup all \
    -g $params.macs_genome_size \
    --call-summits
    """
}
 

/*
 * generate training, validation, and test chromosomes
 */

  process prep_splits {
	container = "kundajelab/chrombpnet:latest"
	
	input:
	path(chrom_sizes)
	
	output:
	path("fold0.json")
	
	script:
    """
    tr ' ' '\n' <<< "$params.chrom_test" >> subset_chroms.txt
    tr ' ' '\n' <<< "$params.chrom_val" >> subset_chroms.txt
    tr ' ' '\n' <<< "$params.chrom_train" >> subset_chroms.txt
    
    awk 'BEGIN {FS=" "} FNR==NR {key=\$1; arrayLookup[key]=\$1; next} {key=\$1; if (arrayLookup[key]) print \$0}' subset_chroms.txt $chrom_sizes > subset.chrom.sizes
    
    
    chrombpnet prep splits \
    -c subset.chrom.sizes \
    -tcr $params.chrom_test \
    -vcr $params.chrom_val \
    -op fold0
    """
}

/*
 * generate non-peak regions for training bias model
 */

  process prep_nonpeaks {
	tag "$meta.sample"
	container = "kundajelab/chrombpnet:latest"
	
	input:
	path(peaks)
	path(fasta)
	path(chrom_sizes)
	path(chrom_splits)
	path(exclude_regions)
	 
	
	output:
	path("background_negatives.bed")
	
	script:
    """
    chrombpnet prep nonpeaks \
     -g $fasta \
     -p $peaks \
     -c  $chrom_sizes \
     -fl $chrom_splits \
     -br $exclude_regions \
     -o background
    """
}

/*
 * train bias model
 */
 
 process train_bias {
	container = "kundajelab/chrombpnet:latest"
	
	input:
	path(bam)
	path(bam_index)
	path(peaks)
	path(background_regions)
	path(fasta)
	path(chrom_sizes)
	path(chrom_splits)
	
	
	output:
	tuple val(meta), path("${meta.sample}_peaks.narrowPeak")
	
	script:
    """
    chrombpnet bias pipeline \
        -ibam ~/chrombpnet_tutorial/data/downloads/merged.bam \
        -d "ATAC" \
        -g $genome_fasta \
        -c $chrom_sizes \
        -p $peaks \
        -n $background_regions \
        -fl $chrom_splits \
        -b 0.5 \
        -o . \
        -fp $meta.sample
    """
 
/*
 * train chrombpnet
 */

//  process train_chrombpnet {
// 	tag "$meta.sample"
// 	container = "kundajelab/chrombpnet:latest"
// 	
// 	input:
// 	tuple val(meta), path(bam), path(bam_index)
// 	
// 	output:
// 	tuple val(meta), 
// 	path("${meta.sample}_peaks.narrowPeak")
// 	
// 	script:
//     """
//     chrombpnet pipeline \
//         -ibam ~/chrombpnet_tutorial/data/downloads/merged.bam \
//         -d "ATAC" \
//         -g ~/chrombpnet_tutorial/data/downloads/hg38.fa \
//         -c ~/chrombpnet_tutorial/data/downloads/hg38.chrom.sizes \
//         -p ~/chrombpnet_tutorial/data/peaks_no_blacklist.bed \
//         -n ~/chrombpnet_tutorial/data/output_negatives.bed \
//         -fl ~/chrombpnet_tutorial/data/splits/fold_0.json \
//         -b ~/chrombpnet_tutorial/bias_model/ENCSR868FGK_bias_fold_0.h5 \
//         -o ~/chrombpnet_tutorial/chrombpnet_model/
//     """

/*
 * generate predicted bigWigs
 */

/*
 * generate contribution score bigwigs
 */

/*
 *run tf-modisco to discover motifs with high contribution scores
 */

/*
 * find motif instances
 */

workflow {

// 	bam_ch = Channel.fromPath(params.samplesheet)
// 	| splitCsv( header:true )
//     | map { bam_row ->
//         bam_meta = bam_row.subMap('sample', 'id')
//         [
//         	bam_meta, 
//         	file(bam_row.reads, checkIfExists: true),
//             file(bam_row.read_index, checkIfExists: true)]
//     }
//     | map { meta, bam, index -> [meta.subMap('sample'), bam, index]}
//     | groupTuple 
//     | merge_bam
//     | call_peaks
//     | view
	
	prep_splits("${baseDir}/${params.chrom_sizes}")
	prep_nonpeaks()
	
// 	fq_input_ch = Channel.fromPath(params.samplesheet)
// 	| splitCsv( header:true )
// 	| filter { it["reads"].endsWith("fastq.gz") }
//     | map { fq_row ->
//         fq_meta = fq_row.subMap('sample')
//         [
//         	fq_meta, 
//         	file(fq_row.reads, checkIfExists: true)]
//     }
// 
// 
//        
//     fq_ch = bam2fastq(bam_ch)
//     .mix(fq_input_ch)
// 
//     split_fq_ch = split_fastq(fq_ch)
//     .transpose()
//     
//     polylox_ch = parse_polylox_barcodes(split_fq_ch)
//     | groupTuple
//     | merge_polylox_barcodes
//     | compute_pgen
}
