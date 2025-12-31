nextflow.enable.dsl=2

process GATK4_HAPLOTYPECALLER {
    tag "${meta.library_name}"

    publishDir  "${params.test_data}/haplotypecaller/output", mode: "copy"

    input:
    val meta
    path bam
    path bai
    val intervals
    val ref_fasta
    val ref_fai
    val ref_dict
    val gatk

    output:
    tuple val(meta), path("*.g.vcf.gz"), emit: gvcf
    tuple val(meta), path("*.g.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${meta.library_name}"
    def interval_command = intervals ? "--intervals $intervals" : ""

    def avail_mem = 8
    if (!task.memory) {
        log.info '[GATK HaplotypeCaller] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
    } else {
        avail_mem = (task.memory.giga*0.8).intValue()
    }

    """
    module load python/3.8.12
    # Ensure Python3 is available and create a symlink if necessary
    if command -v python3 >/dev/null 2>&1; then
        ln -sf \$(which python3) ./python
        export PATH=./:\$PATH
    else
        echo "Error: Python3 is not available after loading the module" >&2
        exit 1
    fi

    module load $gatk

    gatk --java-options "-Xmx${avail_mem}G -XX:-UsePerfData" \\
        HaplotypeCaller \\
        --input $bam \\
        --output ${prefix}.g.vcf.gz \\
        --reference $ref_fasta \\
        --emit-ref-confidence GVCF \\
        $interval_command \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = "${meta.library_name}"
    """
    touch ${prefix}.g.vcf.gz
    touch ${prefix}.g.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
