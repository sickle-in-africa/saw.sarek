#!/usr/bin/env nextflow
/*
 *  GET RAW READS QUALITY REPORTS
 *  =============================
 *
 *  
 ***************************************/
nextflow.enable.dsl=2

include {
    GetSoftwareVersions;
} from "${params.modulesDir}/inputs.nf"

include {
    initializeInputChannelsForRawReadQualityReporting
} from "${params.modulesDir}/inputs.nf"

include {
    branchReadGroupsIntoBamOrFastqChannels;
} from "${params.modulesDir}/preprocess.nf"

include {
    GetFastqcQualityReport;
    GetUnmappedBamQualityReport;
    SaveCohortRawReadsQualityReport;
} from "${params.modulesDir}/qualityReports.nf"

workflow {

    (readGroupsFromInput,\
     referenceSequenceFasta,
     bwaIndexTuple,
     referenceSequenceIndex,
     __genderMap__,
     __statusMap__,
     ch_multiqc_config)\
        = initializeInputChannelsForRawReadQualityReporting()

    softwareVersions = GetSoftwareVersions().collect()

    (readGroupsAsBam,\
     readGroupsAsFastq)\
        = branchReadGroupsIntoBamOrFastqChannels(\
            readGroupsFromInput)

    fastqQualityReports\
        = GetFastqcQualityReport(
            readGroupsAsFastq)

    unmappedBamQualityReports\
        = GetUnmappedBamQualityReport(
            readGroupsAsBam)

    rawReadGroupQualityReports\
        = fastqQualityReports
            .mix(unmappedBamQualityReports)
            .collect()

    SaveCohortRawReadsQualityReport(
            ch_multiqc_config,
            softwareVersions,
            rawReadGroupQualityReports)
}