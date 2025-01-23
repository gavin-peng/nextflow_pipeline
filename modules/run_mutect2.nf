nextflow.enable.dsl=2

process GATK4_MUTECT2 {
    tag "${tumor_meta.library_name}"

    publishDir  "${params.test_data}/mutect2/output", mode: "copy"

    input:
    val tumor_meta
    path input_bam
    val intervals
    val ref_fasta
    val ref_fai
    val ref_dict
    val germline_resource
    val germline_resource_tbi
    val panel_of_normals
    val panel_of_normals_tbi
    val gatk

    output:
    tuple val(tumor_meta), path("*.vcf.gz")     , emit: vcf
    tuple val(tumor_meta), path("*.tbi")        , emit: tbi
    tuple val(tumor_meta), path("*.stats")      , emit: stats
    tuple val(tumor_meta), path("*.f1r2.tar.gz"), optional:true, emit: f1r2
    path "versions.yml"                   , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def prefix = "${tumor_meta.library_name}_${params.mutect2.mutectTag}"
    def input_bams = input_bam.collect{ "--input $it"}.join(" ")
    def interval_command = intervals ? "--intervals $intervals" : ""
    def pon_command = panel_of_normals != 'NO_PON' ?  "--panel-of-normals $panel_of_normals" : ""
    def gr_command = germline_resource != 'no_gnomad' ? "--germline-resource $germline_resource" : ""

    def avail_mem = 8
    if (!task.memory) {
        log.info '[GATK Mutect2] Available memory not known - defaulting to 8GB. Specify process memory requirements to change this.'
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
        Mutect2 \\
        $input_bams \\
        --output ${prefix}.vcf.gz \\
        --reference $ref_fasta \\
        $pon_command \\
        $gr_command \\
        $interval_command \\
        --tmp-dir . \\
        $params.mutect2.mutect2ExtraArgs

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
    stub:
    def prefix = task.ext.prefix ?: "${tumor_meta.library_name}"
    """
    touch ${prefix}.vcf.gz
    touch ${prefix}.vcf.gz.tbi
    touch ${prefix}.vcf.gz.stats
    touch ${prefix}.f1r2.tar.gz

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        gatk4: \$(echo \$(gatk --version 2>&1) | sed 's/^.*(GATK) v//; s/ .*\$//')
    END_VERSIONS
    """
}