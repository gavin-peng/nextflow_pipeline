nextflow.enable.dsl=2

process RUN_VARSCAN {
    tag "${tumor_meta.library_name}"

    publishDir  "${params.test_data}/varscan/output", mode: "copy"

    input:
    val tumor_meta
    path tumor_bam
    path normal_bam
    val intervals
    val ref_fasta
    val ref_fai
    val ref_modules
    val varscan_module

    output:
    tuple val(tumor_meta), path("*.snp.vcf.gz"), emit: snp
    tuple val(tumor_meta), path("*.snp.vcf.gz.tbi"), emit: snp_tbi
    tuple val(tumor_meta), path("*.indel.vcf.gz"), emit: indel
    tuple val(tumor_meta), path("*.indel.vcf.gz.tbi"), emit: indel_tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${tumor_meta.library_name}"
    def interval_command = intervals ? "-l $intervals" : ""
    def min_coverage = params.varscan.min_coverage ?: 8
    def min_var_freq = params.varscan.min_var_freq ?: 0.01

    def avail_mem = 8
    if (!task.memory) {
        log.info '[VarScan] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    def module_list = ref_modules.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module.trim()}" }.join('\n')

    """
    ${module_load_cmds}
    module load ${varscan_module}
    module load tabix/0.2.6

    # Generate normal pileup
    samtools mpileup \\
        -f ${ref_fasta} \\
        ${interval_command} \\
        -q 1 \\
        -B \\
        ${normal_bam} > normal.pileup

    # Generate tumor pileup
    samtools mpileup \\
        -f ${ref_fasta} \\
        ${interval_command} \\
        -q 1 \\
        -B \\
        ${tumor_bam} > tumor.pileup

    # Run VarScan somatic
    varscan somatic \\
        normal.pileup \\
        tumor.pileup \\
        ${prefix} \\
        --min-coverage ${min_coverage} \\
        --min-var-freq ${min_var_freq} \\
        --output-vcf 1 \\
        ${args}

    # Compress and index outputs
    if [ -f "${prefix}.snp.vcf" ]; then
        bgzip -c ${prefix}.snp.vcf > ${prefix}.snp.vcf.gz
        tabix -p vcf ${prefix}.snp.vcf.gz
    fi

    if [ -f "${prefix}.indel.vcf" ]; then
        bgzip -c ${prefix}.indel.vcf > ${prefix}.indel.vcf.gz
        tabix -p vcf ${prefix}.indel.vcf.gz
    fi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(echo \$(varscan 2>&1) | sed 's/^.*VarScan v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${tumor_meta.library_name}"
    """
    touch ${prefix}.snp.vcf.gz
    touch ${prefix}.snp.vcf.gz.tbi
    touch ${prefix}.indel.vcf.gz
    touch ${prefix}.indel.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        varscan: \$(echo \$(varscan 2>&1) | sed 's/^.*VarScan v//; s/ .*\$//')
        samtools: \$(echo \$(samtools --version 2>&1) | sed 's/^.*samtools //; s/Using.*\$//')
    END_VERSIONS
    """
}
