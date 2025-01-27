nextflow.enable.dsl=2
process BWA_MEM {
    tag "$meta.pair_id"

    publishDir "${params.test_data}/bwamem/output", mode: 'copy'

    input:
    tuple val(meta), path(read1), path(read2), val(readGroups), val(fasta), val(sort_bam), val(threads), val(modules), val(addParem)

    output:
    tuple val(meta), path("*.bam"), emit: bam, optional: true
    tuple val(meta), path("*.bai"), emit: bai, optional: true
    tuple val(meta), path("*.cram"), emit: cram, optional: true
    tuple val(meta), path("*.csi"), emit: csi, optional: true
    tuple val(meta), path("*.crai"), emit: crai, optional: true

    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.pair_id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")    ? "sam" :
                    args2.contains("--output-fmt cram")   ? "cram":
                    sort_bam && args2.contains("-O cram") ? "cram":
                    !sort_bam && args2.contains("-C")     ? "cram":
                    "bam"
    def reference = fasta && extension=="cram" ? "--reference ${fasta}" : ""
    if (!fasta && extension=="cram") error "Fasta reference is required for CRAM output"
    def module_list = modules.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module}" }.join('\n')

    """
    ${module_load_cmds}
    echo "readgroup:" "$readGroups"

    bwa mem -M \\
        -R "$readGroups" \\
        $fasta \\
        $read1 \\
        $read2 \\
        | samtools $samtools_command $args2  -o ${prefix}.${extension} -

    samtools index ${prefix}.${extension}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    
    """

    stub:
    def args = task.ext.args ?: ''
    def args2 = task.ext.args2 ?: ''
    def prefix = task.ext.prefix ?: "${meta.pair_id}"
    def samtools_command = sort_bam ? 'sort' : 'view'
    def extension = args2.contains("--output-fmt sam")   ? "sam" :
                    args2.contains("--output-fmt cram")  ? "cram":
                    sort_bam && args2.contains("-O cram")? "cram":
                    !sort_bam && args2.contains("-C")    ? "cram":
                    "bam"
    """
    touch ${prefix}.${extension}
    touch ${prefix}.csi
    touch ${prefix}.crai

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        bwa: \$(echo \$(bwa 2>&1) | sed 's/^.*Version: //; s/Contact:.*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
