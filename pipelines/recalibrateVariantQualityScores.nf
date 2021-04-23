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
     hapmapIndex
     onekgSnps,
     onekgSnpsIndex,
     onekgIndels,
     onekgIndelsIndex,
     onekgOmni,
     onekgOmniIndex)\
        = initializeInputChannelsForVariantRecalibration()


    variantRecalibrationReports\
        = GetVariantRecalibrationReport(
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

}