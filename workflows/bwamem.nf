nextflow.enable.dsl=2
include {BWA_MEM} from '../modules/run_bwamem'

workflow bwaMem {
    take:
    fastq_inputs
    sort_bam
    threads
    addParem

    main:

    def bwaMem_modules_by_genome = [
        hg19: "samtools/1.9 bwa/0.7.17 hg19-bwa-index/0.7.17",
        hg38: "samtools/1.9 bwa/0.7.17 hg38-bwa-index-with-alt/0.7.17",
        mm10: "samtools/1.9 bwa/0.7.17 mm10-bwa-index/0.7.17"
    ]
    def bwaMemRef_by_genome = [
        hg19: "/.mounts/labs/gsi/modulator/sw/data/hg19-bwa-index-0.7.17/hg19_random.fa",
        hg38: "/.mounts/labs/gsi/modulator/sw/data/hg38-bwa-index-with-alt-0.7.17/hg38_random.fa",
        mm10: "/.mounts/labs/gsi/modulator/sw/data/mm10-bwa-index-0.7.17/mm10.fa"
    ]


    process_inputs = fastq_inputs.combine(sort_bam)
                                 .combine(threads)
                                 .combine(addParem)
                                 .map { meta, read1, read2, rg, ref, sb, t, ap ->
                                     def bwaMemRef = bwaMemRef_by_genome[ref]
                                     def modules = bwaMem_modules_by_genome[ref]
                                     [meta, read1, read2, rg, bwaMemRef, sb, t, modules, ap]
                                 }

    //process_inputs.view { it -> "Debug - Combined input: $it" }

    BWA_MEM(process_inputs)
    
    emit:
    bam = BWA_MEM.out.bam
    bai = BWA_MEM.out.bai
}