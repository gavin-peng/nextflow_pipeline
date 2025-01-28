nextflow.enable.dsl=2
include {ENSEMBLVEP_VEP} from "../modules/run_vep"

workflow vep {

    take:
    vcf_input

    main:

    def GenomeResources = [
        hg19: [
            vep_modules: "vep/105.0 tabix/0.2.6 vep-hg19-cache/105 hg19/p13",
            vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 hg19/p13 vep-hg19-cache/105",
            vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-hg19-cache-105/.vep",
            referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/hg19-p13/hg19_random.fa",
            species: "homo_sapiens",
            ncbiBuild: "GRCh37",
            vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/",
            customTranscriptFile :"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.GRCh38.v1.3.ensembl_genomic.filtered.gtf.gz",
            customTranscriptENSTids:"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.filtered.ENST.ID.txt"
        ],
        hg38: [
            vep_modules: "vep/105.0 tabix/0.2.6 vep-hg38-cache/105 hg38/p12",
            vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 hg38/p12 vep-hg38-cache/105",
            vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-hg38-cache-105/.vep",
            referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/hg38-p12/hg38_random.fa",
            species: "homo_sapiens",
            ncbiBuild: "GRCh38",
            vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/",
            customTranscriptFile :"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.GRCh38.v1.3.ensembl_genomic.filtered.gtf.gz",
            customTranscriptENSTids:"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.filtered.ENST.ID.txt"
        ],
        mm39: [
            vep_modules: "vep/105.0 tabix/0.2.6 vep-mm39-cache/105 mm39/p6",
            vcf2maf_modules: "vcf2maf/1.6.21b tabix/0.2.6 mm39/p6 vep-mm39-cache/105",
            vepCacheDir: "/.mounts/labs/gsi/modulator/sw/data/vep-mm39-cache-105/.vep",
            referenceFasta: "/.mounts/labs/gsi/modulator/sw/data/mm39-p6/mm39.fa",
            species: "mus_musculus",
            ncbiBuild: "GRCm39",
            vepPath: "/.mounts/labs/gsi/modulator/sw/Ubuntu20.04/vep-105.0/bin/",
            customTranscriptFile :"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.GRCh38.v1.3.ensembl_genomic.filtered.gtf.gz",
            customTranscriptENSTids:"/.mounts/labs/gsi/src/variantEffectPredictor/mane/MANE.filtered.ENST.ID.txt"
        ]
    ]

    vcf_input
    .map { tuple ->
        def ref = tuple[2]
        [
            GenomeResources[ref]['referenceFasta'],
            GenomeResources[ref]['vepCacheDir'],
            GenomeResources[ref]['ncbiBuild'],
            GenomeResources[ref]['species'],
            GenomeResources[ref]['vep_modules'],
            GenomeResources[ref]['customTranscriptFile']
        ]
    }
    .set { vep_params }

    ENSEMBLVEP_VEP(
        vcf_input,
        vep_params
    )
}