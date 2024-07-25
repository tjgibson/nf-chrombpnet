#! /usr/bin/env nextflow

log.info """
	Chrombpnet Pipeline
	===================================
	samplesheet: ${params.samplesheet}
	results_dir: ${params.results_dir}
	"""
	.stripIndent()

/*
 * generate empty dummy files to use for optional parameters
 */


/*
 * merge sorted and indexed input bam files by sample
 */
 
 // process merge_bam {
// 	tag "$meta.sample"
// 	conda "bioconda::samtools=1.20"
// 	container = "quay.io/biocontainers/samtools:1.20--h50ea8bc_0"
// 	
// 	input:
// 	tuple val(meta), path(bam), path(bam_index)
// 	
// 	output:
// 	tuple val(meta), 
// 	path("${meta.sample}_merged.bam"), 
// 	path("${meta.sample}_merged.bam.bai")
// 	
// 	script:
//     """
//     samtools merge -o ${meta.sample}_merged.bam --threads $task.cpus $bam
//     samtools index ${meta.sample}_merged.bam
//     """
// }
 
 
/*
 * Call peaks on merged samples
 */
 
 //  process call_peaks {
// 	tag "$meta.sample"
// 	conda "bioconda::macs3=3.01"
// 	container = "quay.io/biocontainers/macs3:3.0.1--py312he57d009_3"
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
//     macs3 callpeak \
//     -t $bam \
//     -n $meta.sample \
//     -f BAMPE \
//     --keep-dup all \
//     -g $params.macs_genome_size \
//     --call-summits
//     """
// }
 

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
 * combine peaks, blacklist, and exclude regions to get combined list of exclusion regions for generating nonpeaks
 */

  process combine_exclude_regions {
	container = "quay.io/biocontainers/bedtools:2.31.1--hf5e1c6e_2"
	
	input:
	path("peaks/*_peaks.txt")
	 
	
	output:
	path("combined_exclude_regions.bed")
	
	script:
    """
    cat peaks/*_peaks.txt | cut -f 1,2,3 > combined_exclude_regions.bed
    """
}

/*
 * generate non-peak regions for training bias model
 */

  process prep_nonpeaks {
	container = "kundajelab/chrombpnet:latest"
	
	input:
	tuple val(meta), path(bam), path(bam_index), path(peaks)
	path(fasta)
	path(chrom_sizes)
	path(chrom_splits)
	path(exclude_regions)
	 
	
	output:
	path("background_negatives.bed")
	
	script:
    """
    chrombpnet prep nonpeaks \
 	 -p $peaks \
     -g $fasta \
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
    publishDir "${params.results_dir}/bias_model/", mode: 'copy'
	container = "kundajelab/chrombpnet:latest"
	
	input:
	tuple val(meta), path(bam), path(bam_index), path(peaks)
	path(background_regions)
	path(fasta)
	path(chrom_sizes)
	path(chrom_splits)
	
	
	output:
	path("models/bias.h5"), emit: model
	tuple path("evaluation/overall_report.html"), path("evaluation/overall_report.pdf"), emit: reports
	
	script:
    """
    chrombpnet bias pipeline \
        -ibam $bam \
        -d "ATAC" \
        -g $fasta \
        -c $chrom_sizes \
        -p $peaks \
        -n $background_regions \
        -fl $chrom_splits \
        -b 0.5 \
        -o . \
        -fp bias
    """
 }
 
/*
 * train chrombpnet
 */

 process train_chrombpnet {
	tag "$meta.sample"
	publishDir "${params.results_dir}/${meta.sample}/", mode: 'copy'
	container = "kundajelab/chrombpnet:latest"
	
	input:
	tuple val(meta), path(bam), path(bam_index), path(peaks)
	path(background_regions)
	path(fasta)
	path(chrom_sizes)
	path(chrom_splits)
	path(bias_model)
	
	output:
	tuple val(meta), 
	path("models/*.h5")
	tuple path("evaluation/overall_report.html"), path("evaluation/overall_report.pdf"), emit: reports
	
	script:
    """
    chrombpnet pipeline \
        -ibam $bam \
        -d "ATAC" \
        -g $fasta \
        -c $chrom_sizes \
        -p $peaks \
        -n $background_regions \
        -fl $chrom_splits \
        -b $bias_model \
        -o .
    """
}

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

	
	prep_splits("${baseDir}/${params.chrom_sizes}")
	
	if (params.blacklist) {
	blacklist_ch = channel.fromPath(params.blacklist)
	}
	else {
	blacklist_ch = channel.empty()
	}
	
	if (params.background_exclude_regions) {
	exclude_ch = channel.fromPath(params.background_exclude_regions)
	}
	else {
	exclude_ch = channel.empty()
	}
	
	
	combine_exclude_ch = Channel.fromPath(params.samplesheet)
	| splitCsv(header:true)
	| map { row -> 
	peaks = [file(row.peaks, checkIfExists: true)] 
	}
	| mix(blacklist_ch, exclude_ch)
	| collect
	
	combine_exclude_regions(combine_exclude_ch)
	

	bias_ch = Channel.fromPath(params.samplesheet)
	| splitCsv(header:true)
	| filter { it["sample"] == params.bias_sample }
	 | map { bam_row ->
        bam_meta = bam_row.subMap('sample')
        [
        	bam_meta, 
        	file(bam_row.reads, checkIfExists: true),
            file(bam_row.read_index, checkIfExists: true),
            file(bam_row.peaks, checkIfExists: true)]
    }

	
	prep_nonpeaks(
	bias_ch, 
	"${baseDir}/${params.fasta}", 
	"${baseDir}/${params.chrom_sizes}", 
	prep_splits.out,
	combine_exclude_regions.out
	)

	train_bias(
	bias_ch,
	prep_nonpeaks.out,
	"${baseDir}/${params.fasta}", 
	"${baseDir}/${params.chrom_sizes}", 
	prep_splits.out
	)
	
	train_ch = Channel.fromPath(params.samplesheet)
	| splitCsv(header:true)
	 | map { bam_row ->
        bam_meta = bam_row.subMap('sample')
        [
        	bam_meta, 
        	file(bam_row.reads, checkIfExists: true),
            file(bam_row.read_index, checkIfExists: true),
            file(bam_row.peaks, checkIfExists: true)]
    }
	
		train_chrombpnet(
	train_ch,
	prep_nonpeaks.out,
	"${baseDir}/${params.fasta}", 
	"${baseDir}/${params.chrom_sizes}", 
	prep_splits.out,
	train_bias.out.model
	)
}
