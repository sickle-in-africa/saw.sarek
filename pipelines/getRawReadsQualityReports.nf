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
} from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRawReadQualityReporting
} from "${params.modulesDir}/inputs.nf"

include {
    branchReadGroupsIntoBamOrFastqChannels;
} from "${params.modulesDir}/preprocess.nf"

include {
    GetFastqcQualityReport;
    GetUnmappedBamQualityReport;
    GetCohortRawReadsQualityReport;
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

    ch_software_versions_yaml\
        = GetSoftwareVersions()
        .collect()  

    SaveCohortRawReadsQualityReport(
            ch_multiqc_config,
            ch_software_versions_yaml,
            rawReadGroupQualityReports.collect())
}