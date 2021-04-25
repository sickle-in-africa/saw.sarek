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
    PrintMessage;
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

/*

    RecalibrateIndelQualityScores(
        variantSetsFromInput,
        indelRecalibrationReports)

    /*

    /**/

}

/*
process PrintMessage {

    input:
        tuple val(variantCaller), val(idSample), path(vcf)

    output:
        stdout

    script:
        """
        echo ${vcf}
        """
}
*/
