#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include {
    groupByPatientSample
} from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    GetBwaIndexes;
    GetGatkDictionary;
    GetReferenceSequenceIndex;
} from "${params.modulesDir}/indices.nf"

include {
    GetIntervalsPlan;
    GetIntervals;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    RecalibrateBasesInReadGroup;
    MergeRecalibratedReadGroupsForSample;
    IndexRecalibratedSampleReadGoup;
    writeTsvFilesForRecalibratedBams
} from "${params.modulesDir}/recalibration.nf"


workflow {

    (sampleReadGroups,\
     igenomesReferenceSequenceFasta,
     igenomesReferenceSequenceDictionary,
     igenomesReferenceSequenceIndex,
     dbsnp,
     intervalsPlanFromInput,
     knownIndels,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForRecalibration()

    // get reference indexes as channels

    referenceSequenceDictionary\
        = GetGatkDictionary(\
            igenomesReferenceSequenceFasta,\
            igenomesReferenceSequenceDictionary)

    referenceSequenceIndex\
        = GetReferenceSequenceIndex(\
            igenomesReferenceSequenceFasta,\
            igenomesReferenceSequenceIndex)

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
            igenomesReferenceSequenceFasta,\
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