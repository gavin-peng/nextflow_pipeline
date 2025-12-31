nextflow.enable.dsl=2

process GATK4_GENOTYPEGVCFS {
    tag "$donor"

    publishDir  "${params.test_data}/haplotypecaller/output", mode: "copy"

    input:
    val donor
    path gvcf
    path tbi
    val ref_fasta
    val ref_fai
    val ref_dict
    val gatk

    output:
    tuple val(donor), path("*.vcf.gz"), emit: vcf
    tuple val(donor), path("*.vcf.gz.tbi"), emit: tbi
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''

    def avail_mem = 8
    if (!task.memory) {
        log.info '[GATK GenotypeGVCFs] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
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
        GenotypeGVCFs \\
        --variant $gvcf \\
        --output ${donor}.vcf.gz \\
        --reference $ref_fasta \\
        --tmp-dir . \\
        $args

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """

    stub:
    """
    touch ${donor}.vcf.gz
    touch ${donor}.vcf.gz.tbi

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}
