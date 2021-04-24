#!/usr/bin/env nextflow
/*
 *  RECALIBRATE VARIANT QUALITY SCORES FOR SAMPLE VARIANT SETS
 *  ==========================================================
 *
 *
 ***************************************************************/

nextflow.enable.dsl=2

include {
    initializeInputChannelsForVariantRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    GetIndelRecalibrationReport;
    RecalibrateIndelQualityScores;
    GetVariantRecalibrationReport;
    RecalibrateVariantQualityScores
} from "${params.modulesDir}/variantRecalibration.nf"

workflow {

    (variantSetsFromInput,
     referenceSequenceFasta,
     referenceSequenceDictionary,
     referenceSequenceIndex,
     dbsnp,
     dbsnpIndex,
     hapmap,
     hapmapIndex,
     onekgSnps,
     onekgSnpsIndex,
     onekgIndels,
     onekgIndelsIndex,
     onekgOmni,
     onekgOmniIndex,
     axiomExomePlus,
     axiomExomePlusIndex)\
        = initializeInputChannelsForVariantRecalibration()


    indelRecalibrationReports\
        = GetIndelRecalibrationReport(
            variantSetsFromInput,
            onekgIndels,
            onekgIndelsIndex,
            axiomExomePlus,
            axiomExomePlusIndex,
            dbsnp,
            dbsnpIndex)

    RecalibrateIndelQualityScores(
        variantSetsFromInput,
        indelRecalibrationReports)

    /*

    variantRecalibrationReports\
        = GetVariantRecalibrationReport(
            variantSetsFromInput,
            referenceSequenceFasta,
            referenceSequenceDictionary,
            referenceSequenceIndex,
            dbsnp,
            dbsnpIndex,
            hapmap,
            hapmapIndex,
            onekgSnps,
            onekgSnpsIndex,
            onekgIndels,
            onekgIndelsIndex,
            onekgOmni,
            onekgOmniIndex)


    RecalibrateVariantQualityScores(
        variantSetsFromInput,
        variantRecalibrationReports)

    /**/

}