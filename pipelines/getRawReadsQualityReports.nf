#!/usr/bin/env nextflow
/*
 *  GET RAW READS QUALITY REPORTS
 *  =============================
 *
 *  This pipeline takes in raw read groups and performs quality
 *  analysis. Specifically we use the tool fastqc to generate quality
 *  reports for each input read group, and then we use the tool
 *  multiqc to merge the reports from each individual read group into
 *  a composite quality report for the whole cohort.
 *
 *  Note that that by "read group" we mean a single unmapped bam file 
 *  or fastq file pair containing a list of paired-end reads. This may
 *  be the set of all reads for a single sample, or it might be one of
 *  of many files corresponding to multiple sequencing runs of the 
 *  same sample. 
 *  
 ********************************************************************/
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