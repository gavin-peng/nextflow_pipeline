nextflow.enable.dsl = 2

include { PICARD_MARKDUPLICATES } from "../modules/mark_duplicates"
include { RUN_DELLY } from "../modules/run_delly.nf"
include { MERGE_AND_ZIP as mergeAndZipALL; MERGE_AND_ZIP as mergeAndZipFiltered } from "../modules/merge_vcf.nf"

workflow delly {
    take:
        inputBams
        inputBais
        tumorName
        markdup
        reference
        picard_module

    main:
        def GenomeResources = [
            hg19: [
                rundelly_module: "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg19/p13 hg19-delly/1.0",
                rundelly_fasta: "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
                rundelly_exclude_list: "/.mounts/labs/gsi/modulator/sw/data/hg19-delly-1.0/human.hg19.excl.tsv"
            ],
            hg38: [
                rundelly_module: "delly/0.9.1 bcftools/1.9 tabix/0.2.6 hg38/p12 hg38-delly/1.0",
                rundelly_fasta: "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
                rundelly_exclude_list: "/.mounts/labs/gsi/modulator/sw/data/hg38-delly-1.0/human.hg38.excl.tsv"
            ]
        ]
        
        def markedBams = markdup ? PICARD_MARKDUPLICATES(inputBams, picard_module, tumorName).out.bam : inputBams
        def markedBais = markdup ? PICARD_MARKDUPLICATES(inputBams, picard_module, tumorName).out.bai : inputBais
        def bamCount = markedBams.map { it -> it.findAll { file -> file.name.endsWith('.bam') }.size() }
        def callType = bamCount.map { count -> count == 1 ? 'unmatched' : 'somatic' }
        def delly_modes = channel.from("DEL", "DUP", "INV", "INS", "BND")

        RUN_DELLY(
            delly_modes,
            markedBams,
            markedBais,
            callType,
            channel.value(tumorName),
            channel.value(GenomeResources[reference]['rundelly_module']),
            channel.value(GenomeResources[reference]['rundelly_fasta']),
            channel.value(GenomeResources[reference]['rundelly_exclude_list'])
        )

        mergeAndZipALL (
            RUN_DELLY.out.outVcf,
            RUN_DELLY.out.outTbi,
            tumorName,
            callType,
            params.delly.mergeAndZip_modules,
            '_all',
            0
        )

        mergeAndZipFiltered (
            RUN_DELLY.out.outVcf_filtered,
            RUN_DELLY.out.outTbi_filtered,
            tumorName,
            callType,
            params.delly.mergeAndZip_modules,
            '_filtered',
            0
        )

        emit:
        mergedIndex = mergeAndZipALL.out.dellyMergedTabixIndex
        mergedVcf   = mergeAndZipALL.out.dellyMergedVcf
        mergedFilteredIndex = mergeAndZipFiltered.out.dellyMergedTabixIndex
        mergedFilteredVcf   = mergeAndZipFiltered.out.dellyMergedVcf
        mergedFilteredPassIndex = mergeAndZipFiltered.out.dellyMergedPassTabixIndex
        mergedFilteredPassVcf = mergeAndZipFiltered.out.dellyMergedPassVcf
}