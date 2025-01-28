nextflow.enable.dsl=2

/**
 * Creates a read group string in the form: @RG:Z:run_name-IUS_tag_lane_number
 *
 * @param runName    The run name (sequencing run identifier)
 * @param IUSTag     The IUS tag (library or sample identifier)
 * @param laneNumber The lane number 
 * @return           A formatted read group string
 */
def createReadGroup(String runName, String IUSTag, String laneNumber) {
    return "@RG\\tID:${runName}-${IUSTag}_${laneNumber}\\tSM:${runName}-${IUSTag}_${laneNumber}\\tPL:ILLUMINA\\tLB:Lib1"
}

def getSampleName(String donor, String tissue_origin, String tissue_type, String library_type, String group_id) {
    return "${donor}_${tissue_origin}_${tissue_type}_${library_type}_${group_id}"
}

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