nextflow.enable.dsl=2
include {GATK4_MUTECT2} from "../modules/run_mutect2"

workflow mutect2 {

    take:

    tumor_meta
    input_bam
    intervalFile
    pon
    ponIdx

    main:

    def GenomeResources = [
        "hg19": [
            "refDict" : '$HG19_ROOT/hg19_random.dict',
            "refFai" : '$HG19_ROOT/hg19_random.fa.fai',
            "refFasta" : '$HG19_ROOT/hg19_random.fa',
            "modules" : "hg19/p13 samtools/1.9 gatk/4.1.7.0",
            "gnomad": "",
            "gnomadIdx": ""
        ],
        "hg38": [
            "refDict" : '$HG38_ROOT/hg38_random.dict',
            "refFai" : '$HG38_ROOT/hg38_random.fa.fai',
            "refFasta" : '$HG38_ROOT/hg38_random.fa',
            "gnomad": '$HG38_GATK_GNOMAD_ROOT/af-only-gnomad.hg38.vcf.gz',
            "gnomadIdx": '$HG38_GATK_GNOMAD_ROOT/af-only-gnomad.hg38.vcf.gz.tbi',
            "modules" : "hg38/p12 samtools/1.9 hg38-gatk-gnomad/2.0 gatk/4.1.7.0"
        ],
        "mm10": [
            "refDict" : '$MM10_ROOT/mm10.dict',
            "refFai" : '$MM10_ROOT/mm10.fa.fai',
            "refFasta" : '$MM10_ROOT/mm10.fa',
            "modules" : "mm10/p6 samtools/1.9 gatk/4.1.7.0",
            "gnomad": "",
            "gnomadIdx": ""
        ]
    ]

    def reference  = tumor_meta
    .map { meta -> 
        def projectConfig = params.projects[meta.project]
        return projectConfig?.reference ?: 'hg38'
    }
    reference.view{"$it"}
    
    reference
    .map { ref ->
        tuple(
            GenomeResources[ref]['refDict'], 
            GenomeResources[ref]['refFasta'],
            GenomeResources[ref]['refFai'],
            GenomeResources[ref]['gnomad'],
            GenomeResources[ref]['gnomadIdx'],
            GenomeResources[ref]['modules']
        )
    }
    .set { mutect2_params } 

    mutect2_params.view{"$it"}

    GATK4_MUTECT2(
        tumor_meta,
        input_bam,
        intervalFile,
        mutect2_params,
        pon,
        ponIdx
    )

    emit: 
    vcf = GATK4_MUTECT2.out.vcf
    tbi = GATK4_MUTECT2.out.tbi
    stats = GATK4_MUTECT2.out.stats
    
}

