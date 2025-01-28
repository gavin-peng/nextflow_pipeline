nextflow.enable.dsl=2
process ENSEMBLVEP_VEP {
    tag "$meta.library_name"
    publishDir  "${params.test_data}/vep/output/${meta.project}/${meta.library_name}", mode: "copy"

    input:
    tuple val(meta), path(vcf), val(reference), val(target_bed)
    tuple val(fasta), val(cacheDir), val(genome), val(species), val(vep_modules), val(customTranscriptFile)

    output:
    tuple val(meta.library_name), path("*.vcf.gz")  , optional:true, emit: vcf
    tuple val(meta.library_name), path("*.tab.gz")  , optional:true, emit: tab
    tuple val(meta.library_name), path("*.json.gz") , optional:true, emit: json
    path "*.html"                      , optional:true, emit: report
    path "versions.yml"                , emit: versions

    when:
    task.ext.when == null || task.ext.when

    script:
    def args = task.ext.args ?: ''
    def file_extension = args.contains("--vcf") ? 'vcf' : args.contains("--json")? 'json' : args.contains("--tab")? 'tab' : 'vcf'
    def compress_cmd = args.contains("--compress_output") ? '' : '--compress_output bgzip'
    def prefix = task.ext.prefix ?: "${meta.library_name}"
    
    def module_list = vep_modules.split(',')
    def module_load_cmds = module_list.collect { module -> "module load ${module}" }.join('\n')
    def human_only_command_line = species == "homo_sapiens" ? "--polyphen b --af --af_1kg --af_esp --af_gnomad" : ""

    """
    ${module_load_cmds}

    vep \\
        --offline \\
        --dir $cacheDir \\
        -i $vcf \\
        -o ${prefix}.${file_extension}.gz \\
        $args \\
        $compress_cmd \\
        --fasta $fasta \\
        --assembly $genome \\
        --species $species \\
        --no_progress --sift b --ccds --uniprot --hgvs --symbol --numbers --domains --gene_phenotype --mane  \\
        --allele_number --no_escape --xref_refseq --failed 1 --flag_pick_allele --pick_order canonical,tsl,biotype,rank,ccds,length \\
        $human_only_command_line \\
        --pubmed --regulatory \\
        --fork $task.cpus

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """

    stub:
    def prefix = task.ext.prefix ?: "${meta.library_name}"
    """
    echo "" | gzip > ${prefix}.vcf.gz
    echo "" | gzip > ${prefix}.tab.gz
    echo "" | gzip > ${prefix}.json.gz
    touch ${prefix}_summary.html

    cat <<-END_VERSIONS > versions.yml
    "${task.process}":
        ensemblvep: \$( echo \$(vep --help 2>&1) | sed 's/^.*Versions:.*ensembl-vep : //;s/ .*\$//')
    END_VERSIONS
    """
}