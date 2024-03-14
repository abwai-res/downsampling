nextflow.enable.dsl=2

process FASTQC{
    
    publishDir params.out_dir, mode: 'copy'
    
    tag {"$sample_id"}
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "fastqc_${sample_id}"
        
    script:
    """
       mkdir fastqc_${sample_id}
       fastqc -o fastqc_${sample_id} -t 2 ${reads[0]}
       fastqc -o fastqc_${sample_id} -t 2 ${reads[1]}
    """
}

process seqtk{

    tag{"$sample_id"}
    
    publishDir params.out_dir2, mode: 'copy'

    input:
    tuple val(sample_id), path(reads)

    output:
    tuple val(sample_id), file("*.fastq.gz")
    
    script:
    """
    seqtk sample ${reads[0]} 0.25 > seqtk1_${sample_id}.fastq.gz
    seqtk sample ${reads[1]} 0.25 > seqtk2_${sample_id}.fastq.gz
    """
}

process FASTQC2{
    
    publishDir params.out_dir3, mode: 'copy'
    
    tag {"$sample_id"}
    
    input:
    tuple val(sample_id), path(reads)
    
    output:
    path "fastqc_seqtk_${sample_id}"
        
    script:
    """
        mkdir fastqc_seqtk_${sample_id}
        fastqc -o fastqc_seqtk_${sample_id} -t 2 ${reads[0]}
        fastqc -o fastqc_seqtk_${sample_id} -t 2 ${reads[1]}
    """
}

workflow {
    Channel
        .fromFilePairs(params.reads, checkIfExists: true)
        .set { read_pairs }
        
    read_pairs.view()
    
    fastqc_ch = FASTQC(read_pairs)
    fastqc_ch.view()
    
    seqtk_ch = seqtk(read_pairs)
    seqtk_ch.view()
    
    final_fastqc = FASTQC2(seqtk_ch)
}