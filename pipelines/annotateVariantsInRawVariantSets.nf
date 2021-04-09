#!/usr/bin/env nextflow
/*
 *  ANNOTATING VARIANTS IN VARIANT SETS
 *  ===================================
 *
 *
 ***************************************/
nextflow.enable.dsl = 2

include {
    initializeInputChannelsForAnnotation
} from "${params.modulesDir}/inputs.nf"

include {
    AnnotateVariantsWithSnpeff;
    CompressVariantSetFromSnpeff
} from "${params.modulesDir}/annotation.nf"


workflow {

    (variantSetsFromInput,
     snpeff_config,
     snpeff_cache,
     snpeff_db)\
        = initializeInputChannelsForAnnotation()

    (snpeffReport,
     annotatedVariantSetsFromSnpeff)\
        = AnnotateVariantsWithSnpeff(
            variantSetsFromInput,
            snpeff_config,
            snpeff_cache,
            snpeff_db)

    compressVCFsnpEffOut\
        = CompressVariantSetFromSnpeff(
            annotatedVariantSetsFromSnpeff)
}

