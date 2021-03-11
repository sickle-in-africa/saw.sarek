#!/usr/bin/env nextflow
/*
 *  PREPARING THE WORKSPACE FOR BASE RECALIBRATION
 *  ==============================================
 *
 *      This script prepares the user's workspace for recalibrating the aligned
 *      sample read groups (bam files) by generating the necessary base quality
 *      score recalibration reports with gatk. To speed up the process, we first
 *      create intervals and then the reports are generated in parallel across 
 *      each interval. After generating the reports, we merge them across 
 *      intervals for each sample. 
 *
 *      Note on intervals and interval lists
 *      ------------------------------------
 *      Here we define an *interval list* as a single file containing 
 *      intervals of a reference genome sequence (in position coordinates, for 
 *      example GRCh37). A reference interval list is a file containing a set
 *      of intervals that partition the callable parts of the reference. 
 *      This script checks to see if a reference interval list was supplied 
 *      by the user (e.g. in an igenomes repo) and if not it builds one from
 *      the reference fasta index file.
 *
 *      For parallelization, the reference interval list is split up into many 
 *      smaller interval lists (actually one interval per list/file), and then
 *      reports are build on these interval lists in parallel. 
 *
 *****************************************************************************/

nextflow.enable.dsl=2

include { groupByPatientSample } from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    BuildReferenceIntervalList;
    SplitIntervalList;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    GetBaseRecalibrationReport;
    MergeBaseRecalibrationReportsForSample;
    writeTsvFilesForRecalibrationReports
} from "${params.modulesDir}/recalibration.nf"


workflow {

    (sampleReadGroups,
     referenceSequenceFasta,
     referenceSequenceDictionary,
     referenceSequenceIndex,
     dbsnp,
     dbsnpIndex,
     referenceIntervalListFromInput,
     knownIndels,
     knownIndelsIndex,
     __genderMap__,
     __statusMap__)
        = initializeInputChannelsForRecalibration()

    //  set up intervals for recalibrating bases in parallel  //

    referenceIntervalListFromIndex\
        = BuildReferenceIntervalList(\
            referenceSequenceIndex,\
            referenceIntervalListFromInput)

    referenceIntervalList\
        = referenceIntervalListFromInput.mix(referenceIntervalListFromIndex)

    intervalLists\
        = SplitIntervalList(referenceIntervalList)\
        .flatten()

    intervalsWithDurations\
        = addDurationToInterval(intervals)

    sampleReadGroupAndIntervalListPairs\
        = sampleReadGroups.combine(intervalListsWithDurations)

    //  create recalibration reports  //

    baseRecalibrationReports\
        = GetBaseRecalibrationReport(\
            sampleReadGroupAndIntervalListPairs,\
            dbsnp,\
            dbsnpIndex,\
            referenceSequenceFasta,\
            referenceSequenceDictionary,\
            referenceSequenceIndex,\
            knownIndels,\
            knownIndelsIndex)

    sampleGroupsOfBaseRecalibrationReports\
        = groupByPatientSample(baseRecalibrationReports)

    //   merge recalibration reports  //

    (sampleBaseRecalibrationReports,\
     recalReportTSV)\
        = MergeBaseRecalibrationReportsForSample(\
            sampleGroupsOfBaseRecalibrationReports)

    //  write details of outputs to tsv files  //

    writeTsvFilesForRecalibrationReports(\
        recalReportTSV,\
        __genderMap__,\
        __statusMap__)

    /**/

}