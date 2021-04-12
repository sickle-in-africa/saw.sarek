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

    (vepVCF,
     vepReport)\
        = AnnotateVariantsWithVep(
            variantSetsFromInput,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_indels,
            ch_cadd_indels_tbi,
            ch_cadd_wg_snvs,
            ch_cadd_wg_snvs_tbi)

    (vepVCFmerge, vepReportMerge)\
        = MergeVariantSetsFromVepAndSnpeff(
            compressVCFsnpEffOut,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_indels,
            ch_cadd_indels_tbi,
            ch_cadd_wg_snvs,
            ch_cadd_wg_snvs_tbi)

    vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

    compressVCFOutVEP\
        = CompressVariantSetFromVep(
            vcfCompressVCFvep)

}

