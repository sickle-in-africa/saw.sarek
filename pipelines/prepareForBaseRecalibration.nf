#!/usr/bin/env nextflow
/*
 *  PREPARING THE WORKSPACE FOR BASE RECALIBRATION
 *  ==============================================
 *
 *      This script prepares the user's workspace for recalibrating the aligned
 *      sample read groups (bam files) by generating the necessary base quality
 *      score recalibration reports with gatk. To speed up the process, we first
 *      create intervals and then the reportss are generated in parallel across 
 *      each interval. After generating the reports, we merge them across 
 *      intervals for each sample. 
 *
 *****************************************************************************/

nextflow.enable.dsl=2

include { groupByPatientSample } from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    GetIntervalsPlan;
    GetIntervals;
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
     intervalsPlanFromInput,
     knownIndels,
     knownIndelsIndex,
     __genderMap__,
     __statusMap__)
        = initializeInputChannelsForRecalibration()

    //  set up intervals for recalibrating bases in parallel  //

    intervalsPlan\
        = GetIntervalsPlan(\
            referenceSequenceIndex,\
            intervalsPlanFromInput)

    intervals = GetIntervals(intervalsPlan).flatten()

    intervalsWithDurations\
        = addDurationToInterval(intervals)

    sampleReadGroupAndIntervalPairs\
        = sampleReadGroups.combine(intervalsWithDurations)

    //  create recalibration reports  //

    baseRecalibrationReports\
        = GetBaseRecalibrationReport(\
            sampleReadGroupAndIntervalPairs,\
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