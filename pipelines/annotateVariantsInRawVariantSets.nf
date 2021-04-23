#!/usr/bin/env nextflow
/*
 *  ANNOTATING VARIANTS IN VARIANT SETS
 *  ===================================
 *
 *  In this pipeline we annotate variants from the variant calllers.
 *  This step is applied after variant quality score recalibration
 *  to avoid wasting time annotating variants that then do not pass
 *  basic quality filters. 
 *
 *  We use the following annotation tools:
 *      + SnpEff
 *      + VEP (Variant Effect Predictor)
 *
 *  We produce three different annotated vcf files as output:
 *      + a SnpEff annotated file;
 *      + a VEP annotated file;
 *      + a VEP annotated file that uses input from the SnpEff 
 *          annotated file as input.
 *
 ********************************************************************/
nextflow.enable.dsl = 2

include {
    initializeInputChannelsForAnnotation
} from "${params.modulesDir}/inputs.nf"

include {
    AnnotateVariantsWithSnpeff;
    CompressVariantSetFromSnpeff;
    AnnotateVariantsWithVep;
    AnnotateVariantsWithVepAndSnpeff;
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
        = AnnotateVariantsWithVepAndSnpeff(
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

