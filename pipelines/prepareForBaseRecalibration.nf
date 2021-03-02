#!/usr/bin/env nextflow

nextflow.enable.dsl=2

include { isChannelActive } from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    GetBwaIndexes;
    GetGatkDictionary;
    GetReferenceSequenceIndex;
    GetDbsnpIndex;
    GetKnownIndelsIndex
} from "${params.modulesDir}/indices.nf"

include {
    GetIntervalsPlan;
    GetIntervals;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    BaseRecalibrator;
    GatherBQSRReports;
    writeTsvFilesForRecalibrationTables
} from "${params.modulesDir}/recalibration.nf"


workflow {

    (inputSample,\
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

    dbsnpIndex = GetDbsnpIndex(dbsnp)

    knownIndelsIndex = GetKnownIndelsIndex(knownIndels)

    //  set up intervals for recalibrating bases in parallel  //

    intervalsPlan\
        = GetIntervalsPlan(\
            referenceSequenceIndex,\
            intervalsPlanFromInput)

    intervals = GetIntervals(intervalsPlan).flatten()

    intervalsWithDurations\
        = addDurationToInterval(intervals)

    bamBaseRecalibrator\
        = inputSample.combine(intervalsWithDurations)

    //  create recalibration tables  //

    (tableGatherBQSRReports,\
     recalTableTSVnoInt)\
        = BaseRecalibrator(\
            bamBaseRecalibrator,\
            dbsnp,\
            dbsnpIndex,\
            igenomesReferenceSequenceFasta,\
            referenceSequenceDictionary,\
            referenceSequenceIndex,\
            knownIndels,\
            knownIndelsIndex)

    tableGatherBQSRReports = tableGatherBQSRReports.groupTuple(by:[0, 1])

    //   merge recalibration tables  //

    (recalTable,\
     baseRecalibratorReport,\
     recalTableTSV)\
        = GatherBQSRReports(\
            tableGatherBQSRReports)


    //  write details of outputs to tsv files  //

    recalTableSampleTSV = recalTableTSV.mix(recalTableTSVnoInt)

    writeTsvFilesForRecalibrationTables(\
        recalTableTSV,\
        recalTableSampleTSV,\
        __genderMap__,\
        __statusMap__)

    /**/

}