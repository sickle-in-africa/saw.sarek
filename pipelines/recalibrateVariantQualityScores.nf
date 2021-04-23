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

workflow {

    variantSetsFromInput\
        = initializeInputChannelsForVariantRecalibration()


    variantSetsFromInput.view()
}