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
    CompressVariantSetFromSnpeff;
    AnnotateVariantsWithVep;
    MergeVariantSetsFromVepAndSnpeff;
    CompressVariantSetFromVep
} from "${params.modulesDir}/annotation.nf"


workflow {

    (variantSetsFromInput,
     snpeff_config,
     snpeff_cache,
     snpeff_db,
     vep_cache,
     vep_cache_version,
     cadd_cache,
     cadd_indels,
     cadd_indels_tbi,
     cadd_wg_snvs,
     cadd_wg_snvs_tbi)\
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

    (vepVCF,
     vepReport)\
        = AnnotateVariantsWithVep(
            variantSetsFromInput,
            vep_cache,
            vep_cache_version,
            cadd_indels,
            cadd_indels_tbi,
            cadd_wg_snvs,
            cadd_wg_snvs_tbi)

    (vepVCFmerge, vepReportMerge)\
        = MergeVariantSetsFromVepAndSnpeff(
            compressVCFsnpEffOut,
            vep_cache,
            vep_cache_version,
            cadd_indels,
            cadd_indels_tbi,
            cadd_wg_snvs,
            cadd_wg_snvs_tbi)

    vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

    compressVCFOutVEP\
        = CompressVariantSetFromVep(
            vcfCompressVCFvep)

    /**/

}

