nextflow.enable.dsl=2
include {RUN_VARSCAN} from "../modules/run_varscan"

workflow varscan {

    take:

    tumor_meta
    tumor_bam
    normal_bam
    intervalFile
    reference
    varscan_module

    main:

    def GenomeResources = [
        hg19: [
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
                modules : "hg19/p13 samtools/1.9"
        ],
        hg38: [
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
                modules : "hg38/p12 samtools/1.9"
        ],
        mm10: [
                refFai : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.fa",
                modules : "mm10/p6 samtools/1.9"
            ]
    ]
    reference
    .map { ref ->
        return [
            refFasta: GenomeResources[ref]['refFasta'],
            refFai: GenomeResources[ref]['refFai'],
            modules: GenomeResources[ref]['modules']
        ]
    }
    .set { varscan_params }

    RUN_VARSCAN(
        tumor_meta,
        tumor_bam,
        normal_bam,
        intervalFile,
        varscan_params.refFasta,
        varscan_params.refFai,
        varscan_params.modules,
        varscan_module
    )

    emit:
    snp = RUN_VARSCAN.out.snp
    snp_tbi = RUN_VARSCAN.out.snp_tbi
    indel = RUN_VARSCAN.out.indel
    indel_tbi = RUN_VARSCAN.out.indel_tbi
}

