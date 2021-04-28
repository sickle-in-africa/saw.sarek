#!/usr/bin/env nextflow
/*
 *  FILTER VARIANTS IN SAMPLE VARIANT SETS
 *  ======================================
 *
 *  This pipeline is to be used after variant calling if you have
 *  called variants in the per-sample mode (i.e. you have not joint-
 *  called the variants across the cohort). 
 *
 *  We filter the variant sets from each snp & indel variant caller:
 *      + Gatk
 *      + Freebayes
 *      + Strelka
 *  according to the reccomendations/best practices given in each of
 *  the tools' own documentation. 
 *
 *  For joint-called variants the filtering process is more
 *  more complicated, and we provide alternative pipelines for this
 *  task.
 *
 ********************************************************************/
nextflow.enable.dsl=2

include {
    initializeInputChannelsForSampleVariantSetFiltering;
} from "${params.modulesDir}/inputs.nf"

include {
    FilterVariantsFromFreebayes;
    branchIntoGatkStrelkaOrFreebayesChannels;
} from "${params.modulesDir}/variantFiltering.nf"

workflow {

    (variantSets,
     referenceSequenceFasta,
     referenceSequenceDictionary,
     referenceSequenceIndex)\
        = initializeInputChannelsForSampleVariantSetFiltering()

    // branch into GATK, Freebayes, and Strelka channels
    (variantSetsFromGatk,
     variantSetsFromStrelka,
     variantSetsFromFreebayes)\
        = branchIntoGatkStrelkaOrFreebayesChannels(
            variantSets)

    // filter gatk

    // filter strelka

    // filter freebayes
    FilterVariantsFromFreebayes(variantSetsFromFreebayes)


}