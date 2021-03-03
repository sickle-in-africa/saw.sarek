#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    groupByPatientSample
} from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    GetIntervalsPlan;
    GetIntervals;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    RecalibrateBasesInReadGroup;
    MergeRecalibratedReadGroupsForSample;
    writeTsvFilesForRecalibratedBams
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
     __statusMap__)\
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

    //  recalibrate aligned sample reads in parallel over intervals  //

    (readGroupsRecalibrated)\
        = RecalibrateBasesInReadGroup(\
            sampleReadGroupAndIntervalPairs,\
            referenceSequenceDictionary,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    //  group and merge recalibrated bam files for each sample  //

    sampleGroupsOfReadGroupsRecalibrated\
        = groupByPatientSample(\
            readGroupsRecalibrated)

    (sampleReadGroupsRecalibrated,\
     bam_recalibrated_qc,\
     tsv_bam_recalibrated)\
        = MergeRecalibratedReadGroupsForSample(\
            sampleGroupsOfReadGroupsRecalibrated)

    //  write tsv files for the next step  //

    writeTsvFilesForRecalibratedBams(\
        tsv_bam_recalibrated,\
        __genderMap__,\
        __statusMap__)

    /**/

}