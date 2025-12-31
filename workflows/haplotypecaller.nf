nextflow.enable.dsl=2
include {GATK4_HAPLOTYPECALLER} from "../modules/run_haplotypecaller"
include {GATK4_COMBINEGVCFS} from "../modules/combine_gvcfs"
include {GATK4_GENOTYPEGVCFS} from "../modules/genotype_gvcfs"

workflow haplotypecaller {

    take:

    bams_with_meta  // channel: [meta, bam, bai]
    intervalFile
    reference
    gatk

    main:

    def GenomeResources = [
        hg19: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.dict",
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
                modules : "hg19/p13 samtools/1.9"
        ],
        hg38: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.dict",
                refFai : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa.fai",
                refFasta : "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
                modules : "hg38/p12 samtools/1.9"
        ],
        mm10: [
                refDict : "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm10.dict",
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
            refDict: GenomeResources[ref]['refDict'],
            modules: GenomeResources[ref]['modules']
        ]
    }
    .set { hc_params }

    // Step 1: Run HaplotypeCaller on each BAM
    GATK4_HAPLOTYPECALLER(
        bams_with_meta.map { meta, bam, _bai -> meta },
        bams_with_meta.map { _meta, bam, _bai -> bam },
        bams_with_meta.map { _meta, _bam, bai -> bai },
        intervalFile,
        hc_params.refFasta,
        hc_params.refFai,
        hc_params.refDict,
        gatk
    )

    // Step 2: Group GVCFs by donor and combine
    gvcfs_by_donor = GATK4_HAPLOTYPECALLER.out.gvcf
        .map { meta, gvcf -> [meta.donor, gvcf] }
        .groupTuple()

    GATK4_COMBINEGVCFS(
        gvcfs_by_donor.map { donor, _gvcfs -> donor },
        gvcfs_by_donor.map { _donor, gvcfs -> gvcfs },
        hc_params.refFasta,
        hc_params.refFai,
        hc_params.refDict,
        gatk
    )

    // Step 3: Genotype combined GVCFs
    GATK4_GENOTYPEGVCFS(
        GATK4_COMBINEGVCFS.out.gvcf.map { donor, _gvcf -> donor },
        GATK4_COMBINEGVCFS.out.gvcf.map { _donor, gvcf -> gvcf },
        GATK4_COMBINEGVCFS.out.tbi.map { _donor, tbi -> tbi },
        hc_params.refFasta,
        hc_params.refFai,
        hc_params.refDict,
        gatk
    )

    emit:
    gvcf = GATK4_HAPLOTYPECALLER.out.gvcf
    combined_gvcf = GATK4_COMBINEGVCFS.out.gvcf
    vcf = GATK4_GENOTYPEGVCFS.out.vcf
}

