nextflow.enable.dsl=2

process VCF2MAF {
    tag "$tumorName"
    publishDir  "${params.test_data}/vep/output", mode: "copy"

    input:
    val tumorName
    tuple val(meta), path(vcf)
    val genome
    val fasta
    val vcf2maf_modules
    val normalName
    val tumorOnly

    output:
    tuple val(tumorName), path("*.maf"), emit: maf
    path "versions.yml", emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = task.ext.prefix ?: "${tumorName}"
    def normal_id = tumorOnly || normalName == '' ? '' : "--normal-id ${normalName}"

    def module_list = vcf2maf_modules.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module.trim()}" }.join('\n')

    """
    ${module_load_cmds}

    # Decompress VCF if needed
    if [[ ${vcf} == *.gz ]]; then
        gunzip -c ${vcf} > input.vcf
    else
        cp ${vcf} input.vcf
    fi

    vcf2maf.pl \\
        --input-vcf input.vcf \\
        --output-maf ${prefix}.maf \\
        --tumor-id ${tumorName} \\
        ${normal_id} \\
        --ref-fasta ${fasta} \\
        --ncbi-build ${genome} \\
        --vep-data \${VEP_DATA:-/.mounts/labs/gsi/modulator/sw/data/vep-${genome.toLowerCase()}-cache-105/.vep} \\
        --vep-path \${VEP_PATH:-/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/} \\
        ${args}

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: \$(vcf2maf.pl --man | grep -m 1 'vcf2maf' | sed 's/^.*vcf2maf.pl v//' | sed 's/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${tumorName}"
    """
    touch ${prefix}.maf

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        vcf2maf: \$(vcf2maf.pl --man | grep -m 1 'vcf2maf' | sed 's/^.*vcf2maf.pl v//' | sed 's/ .*\$//')
    END_VERSIONS
    """
}
