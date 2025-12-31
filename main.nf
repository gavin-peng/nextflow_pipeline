nextflow.enable.dsl=2

include {BCL2FASTQ} from './workflows/bcl2fastq'
include {bwaMem} from './workflows/bwamem'
include {vep} from "./workflows/vep"
include {mutect2} from "./workflows/mutect2"
include {delly} from './workflows/delly'
include {varscan} from './workflows/varscan'
include {PICARD_MERGESAMFILES} from "./modules/merge_bams"

// Function to parse attributes string into a map
def parseAttributes(attr) {
    attr.split(';').collectEntries { 
        def parts = it.split('=', 2)
        [(parts[0]): parts.size() > 1 ? parts[1] : true]
    }
}

// Function to extract donor from Sample Name
def extractDonor(sampleName) {
    def parts = sampleName.split('_')
    if (parts.size() >= 2) {
        return parts[0] + '_' + parts[1]
    }
    return sampleName  // Return full sample name if it doesn't match expected format
}

workflow {

// Input data source data.tsv is a tsv file with 10 columns separated by tab. Columns are:
// Project    Sample Name    Sample Attributes    Sequencer Run Name    Sequencer Run Attributes    Lane Name       Lane Number     Lane Attributes    IUS Tag    File Path
// All columns except last "File Path" are considered metadata(which would be carried along the input/output channels), but some columns: Sample Attributes, Sequencer Run Attributes, Lane Attributes are themself multiple fields separated by ";".
// Create the fastq_inputs channel
file_inputs = Channel
    .fromPath('tests/data.tsv')
    .splitCsv(sep: '\t', header: true)
    .map { row ->
        // Initialize meta map with simple fields
        def meta = [
            project: row.'Project',
            sample_name: row.'Sample Name',
            sequencer_run_name: row.'Sequencer Run Name',
            lane_name: row.'Lane Name',
            lane_number: row.'Lane Number',
            ius_tag: row.'IUS Tag'
        ]

        // Parse and add nested attributes
        meta += parseAttributes(row.'Sample Attributes')
        meta += parseAttributes(row.'Sequencer Run Attributes')
        meta += parseAttributes(row.'Lane Attributes')

        // Extract donor from Sample Name
        meta.donor = extractDonor(row.'Sample Name')

        // Create library_name
        meta.library_name = [
            meta.donor,
            meta.geo_tissue_origin,
            meta.geo_tissue_type,
            meta.geo_library_source_template_type,
            meta.geo_group_id
        ].join('_')

        // Create tuple with meta and file
        [meta, file(row.'File Path')]
    }


// Create the fastq_inputs channel
fastq_inputs = file_inputs
    .map { meta, file -> 
        def key = meta.library_name
        def value = [meta: meta, file: file]
        return tuple(key, value)
    }
    .groupTuple()
    .flatMap { library_name, group ->
        def meta = group.first().meta
        
        // Group files by their base name (without R1/R2 suffix)
        def fileGroups = group.groupBy { item ->
            item.file.name.replaceFirst(/_R[12]\.fastq\.gz$/, '')
        }
        
        // For each group, find R1 and R2 files and emit as individual pairs
        return fileGroups.collect { _basename, files ->
            def r1 = files.find { it.file.name.endsWith('R1.fastq.gz') }?.file
            def r2 = files.find { it.file.name.endsWith('R2.fastq.gz') }?.file
            
            if (r1 && r2) {
                // Create a new meta map for each pair, including a unique identifier
                def pairMeta = meta.clone()
                pairMeta.pair_id = r1.name.replaceFirst(/_R1\.fastq\.gz$/, '')
                return tuple(pairMeta, r1, r2)
            } else {
                return null
            }
        }
    }
    .filter { it != null }

reference =  file_inputs
    .map { meta, file ->
        def projectConfig = params.projects[meta.project]
        def reference = projectConfig?.reference ?: 'default_reference'
        return reference
    }

bwaMem(
    fastq_inputs,
    channel.of(params.bwamem.readGroups),
    reference,
    channel.of(params.bwamem.sort_bam),
    channel.of(params.bwamem.threads),
    channel.of(params.bwamem.addParem)
)

// Merge bams of multiple lanes of same library
bams_to_merge = bwaMem.out.bam
    .map { meta, bam -> 
        def key = meta.library_name
        return tuple(key, tuple(meta, bam))
    }
    .groupTuple()
    .map { library_name, group ->
        def meta = group.first()[0]  // Get the meta from the first item
        def bams = group.collect { it[1] }  // Collect all BAM files
        return tuple(meta, bams)
    }

PICARD_MERGESAMFILES(bams_to_merge)
 
def tumor_bam = PICARD_MERGESAMFILES.out.bam
    .filter { meta, bam -> 
        meta.geo_tissue_type != 'R'
    }
    .map { meta, bam -> bam
}

def tumor_meta = PICARD_MERGESAMFILES.out.bam
    .filter { meta, bam -> 
        meta.geo_tissue_type != 'R'
    }
    .map { meta, bam -> meta
}

if (params.mutect2.tumor_only_mode) {
    mutect2_bams = tumor_bam
} else {
    def normal_bam = PICARD_MERGESAMFILES.out.bam
    .filter { meta, bam -> 
        meta.geo_tissue_type == 'R'
    }
    .map { meta, bam -> bam
    }

    mutect2_bams = tumor_bam.combine(normal_bam)
}

mutect2(
    tumor_meta,
    mutect2_bams,
    channel.fromPath(params.mutect2.intervals),
    params.mutect2.panel_of_normals ? channel.fromPath(params.mutect2.panel_of_normals) : channel.value('NO_PON'),
    params.mutect2.panel_of_normals_tbi ? channel.fromPath(params.mutect2.panel_of_normals_tbi) :  channel.value('NO_PON_TBI'),
    channel.value(params.mutect2.reference),
    channel.value(params.mutect2.gatk)
)

// VEP annotation downstream of mutect2
tumorName_vep = mutect2.out.vcf.map { meta, _vcf -> meta.library_name }
vcf_vep = mutect2.out.vcf
reference_vep = mutect2.out.vcf.map { meta, _vcf ->
    def projectConfig = params.projects[meta.project]
    return projectConfig?.reference ?: 'hg38'
}
normalName_vep = channel.value('')
vepTumorOnly_vep = channel.value(params.mutect2.tumor_only_mode)
targetBed_vep = channel.value('')
outputMaf_vep = channel.value(params.vep.output_maf)

vep(
    tumorName_vep,
    vcf_vep,
    reference_vep,
    normalName_vep,
    vepTumorOnly_vep,
    targetBed_vep,
    outputMaf_vep
)

// Delly structural variant calling
def tumor_bam_delly = PICARD_MERGESAMFILES.out.bam
    .filter { meta, bam ->
        meta.geo_tissue_type != 'R'
    }

def normal_bam_delly = PICARD_MERGESAMFILES.out.bam
    .filter { meta, bam ->
        meta.geo_tissue_type == 'R'
    }

if (params.delly.tumor_only_mode) {
    delly_bams = tumor_bam_delly.map { _meta, bam -> bam }
    delly_bais = tumor_bam_delly.map { _meta, bam -> file(bam.toString() + '.bai') }
    delly_meta = tumor_bam_delly.map { meta, _bam -> meta }
} else {
    delly_bams = tumor_bam_delly.combine(normal_bam_delly)
        .map { _meta_t, bam_t, _meta_n, bam_n -> [bam_t, bam_n] }
    delly_bais = tumor_bam_delly.combine(normal_bam_delly)
        .map { _meta_t, bam_t, _meta_n, bam_n -> [file(bam_t.toString() + '.bai'), file(bam_n.toString() + '.bai')] }
    delly_meta = tumor_bam_delly.map { meta, _bam -> meta }
}

reference_delly = delly_meta.map { meta ->
    def projectConfig = params.projects[meta.project]
    return projectConfig?.reference ?: 'hg38'
}

tumorName_delly = delly_meta.map { meta -> meta.library_name }

delly(
    delly_bams,
    delly_bais,
    tumorName_delly,
    channel.value(params.delly.markdup),
    reference_delly,
    channel.value(params.delly.picard_module)
)

// VarScan somatic variant calling (tumor-normal paired)
def tumor_bam_varscan = PICARD_MERGESAMFILES.out.bam
    .filter { meta, _bam ->
        meta.geo_tissue_type != 'R'
    }

def normal_bam_varscan = PICARD_MERGESAMFILES.out.bam
    .filter { meta, _bam ->
        meta.geo_tissue_type == 'R'
    }

// Combine tumor and normal BAMs
varscan_bams = tumor_bam_varscan.combine(normal_bam_varscan)
    .map { meta_t, bam_t, _meta_n, bam_n ->
        [meta_t, bam_t, bam_n]
    }

// Extract metadata and reference
varscan_meta = varscan_bams.map { meta, _t_bam, _n_bam -> meta }
reference_varscan = varscan_meta.map { meta ->
    def projectConfig = params.projects[meta.project]
    return projectConfig?.reference ?: 'hg38'
}

// Extract just the BAMs
varscan_tumor_bam = varscan_bams.map { _meta, t_bam, _n_bam -> t_bam }
varscan_normal_bam = varscan_bams.map { _meta, _t_bam, n_bam -> n_bam }

varscan(
    varscan_meta,
    varscan_tumor_bam,
    varscan_normal_bam,
    channel.value(params.varscan.intervals ?: ''),
    reference_varscan,
    channel.value(params.varscan.varscan_module)
)

}