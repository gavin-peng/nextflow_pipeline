nextflow.enable.dsl=2

include {BCL2FASTQ} from './workflows/bcl2fastq'
include {bwaMem} from './workflows/bwamem'
include {vep} from "./workflows/vep"
include {mutect2} from "./workflows/mutect2"
include {delly} from './workflows/delly'
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
// Project Sample Name     Sample Attributes       Sequencer Run Name      Sequencer Run Attributes        Lane Name       Lane Number     Lane Attributes IUS Tag File Path
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
        return fileGroups.collect { _, files ->
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

bwaMem(
    fastq_inputs,
    channel.of(params.bwamem.readGroups),
    channel.of(params.bwamem.reference),
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
   

/*   
    tumor_name = params.vep.tumorName
    reference = params.vep.reference
    normal_name = params.vep.normalName
    vep_tumor_only = params.vep.onlyTumor
    target_bed = params.vep.targetBed

    vep_vcf = channel.fromPath("${params.test_data}/vep/input/*${tumor_name}*.vcf.gz")

    vep(
        channel.value(tumor_name),
        vep_vcf,
        channel.value(reference),
        channel.value(normal_name),
        channel.value(vep_tumor_only),
        channel.value(target_bed)
    )

    tumor_bam_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.tumorName}*.bam")
    tumor_bam_index_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.tumorName}*.bam.bai")

    if (!params.delly.tumor_only_mode) {
        normal_bam_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.normalName}*.bam")
        normal_bam_index_files = channel.fromPath("${params.test_data}/delly/input/*${params.delly.normalName}*.bam.bai")
    }

    if (params.delly.tumor_only_mode) {
        delly_bams = tumor_bam_files
        delly_indexes = tumor_bam_index_files
    } else {
        delly_bams = tumor_bam_files.combine(normal_bam_files)
        delly_indexes = tumor_bam_index_files.combine(normal_bam_index_files)
    }
    
    delly (
        delly_bams,
        delly_indexes,
        params.delly.tumorName,
        params.delly.markdup,
        params.delly.reference,
        params.delly.picard_module
    )
    */
}