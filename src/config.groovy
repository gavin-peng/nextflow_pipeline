/**
 * Creates a read group string in the form: @RG:Z:run_name-IUS_tag_lane_number
 *
 * @param runName    The run name (sequencing run identifier)
 * @param IUSTag     The IUS tag (library or sample identifier)
 * @param laneNumber The lane number 
 * @return           A formatted read group string
 */
def createReadGroup(String runName, String IUSTag, int laneNumber) {
    return "@RG:Z:${runName}-${IUSTag}_${laneNumber}"
}
