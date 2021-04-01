#!/usr/bin/env nextflow
/*
 *  VARIANT CALLING
 *  ===============
 *
 *      This script takes in *aligned sample read groups* (single files of 
 *      reads, such that one read file corresponds to one and only one patient 
 *      sample) in the bam format, creates intervals (one interval per 
 *      interval/bed file) and then calls variants in the sample read groups 
 *      for a range of variant callers. 
 *      
 *      The callers are either for snvs and small indels:
 *          + gatk 4 Haplotype caller
 *          + strelka
 *          + freebayes
 *      or for structural variants:
 *          + manta
 *          + tiddit
 *
 *      Calling variants on intervals is helpful for speeding up the process, 
 *      by parralellising over each interval, however it is only possible for 
 *      the gatk and freebayes callers. Once variant calling is done, all the 
 *      variant sets corresponding to each sample are merged across intervals 
 *      (where intervals were used) and the sample variant sets are output. 
 *
 *      Note on intervals and interval lists
 *      ------------------------------------
 *      Here we define an *interval list* as a single file containing 
 *      intervals of a reference genome sequence (in position coordinates, for 
 *      example GRCh37). A reference interval list is a file containing a set
 *      of intervals that partition the callable parts of the reference. 
 *      This script checks to see if a reference interval list was supplied 
 *      by the user (e.g. in an igenomes repo) and if not it builds one from
 *      the reference fasta index file.
 *
 *      For parallelization, the reference interval list is split up into many 
 *      smaller interval lists (actually one interval per list/file), and then
 *      variants are called on these interval lists in parallel. 
 *
 ******************************************************************************/

nextflow.enable.dsl=2

include {
    groupByPatientSample
} from "${params.modulesDir}/sarek.nf"

include {
    GetSoftwareVersions;
    initializeInputChannelsForCalling
} from "${params.modulesDir}/inputs.nf"

include {
    BuildReferenceIntervalList;
    SplitIntervalList;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    CallVariantsWithGatk;
    GenotypeVariantsFromGatk;
    CallVariantsWithStrelka;
    CallVariantsWithManta;
    CallVariantsWithTiddit;
    CallVariantsWithFreebayes;
    MergeVariantSetsForSample;
    branchIntoGenotypingOrNoGenotypingChannels;
    removeIntervalList;
    writeInputsForNextStep
} from "${params.modulesDir}/calling.nf"


workflow {

    (sampleReadGroups,
     referenceSequenceFasta,
     referenceSequenceDictionary,
     referenceSequenceIndex,
     dbsnp,
     dbsnpIndex,
     referenceIntervalListFromInput,
     targetIntervalFromInput,
     __genderMap__,
     __statusMap__)
        = initializeInputChannelsForCalling()

    softwareVersions = GetSoftwareVersions()

   //  set up intervals for calling variants in parallel  //

    referenceIntervalListFromIndex\
        = BuildReferenceIntervalList(
            referenceSequenceIndex,
            referenceIntervalListFromInput.ifEmpty('empty'))

    referenceIntervalList\
        = referenceIntervalListFromInput.mix(referenceIntervalListFromIndex)

    intervalLists\
        = SplitIntervalList(referenceIntervalList)\
        .flatten()

    intervalListsWithDurations\
        = addDurationToInterval(intervalLists)

    sampleReadGroupAndIntervalListPairs\
        = sampleReadGroups.combine(intervalListsWithDurations)

    //  call variants  //

    variantSetAndIntervalPairsFromGatk\
        = CallVariantsWithGatk(
            sampleReadGroupAndIntervalListPairs,
            dbsnp,
            dbsnpIndex,
            referenceSequenceDictionary,
            referenceSequenceFasta,
            referenceSequenceIndex)

    (forGenotyping,\
     noGenotyping)\
        = branchIntoGenotypingOrNoGenotypingChannels(\
            variantSetAndIntervalPairsFromGatk)

    (genotyped)\
        = GenotypeVariantsFromGatk(\
            forGenotyping,\
            dbsnp,\
            dbsnpIndex,\
            referenceSequenceDictionary,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    variantSetAndIntervalListPairsFromGatk\
        = genotyped.mix(noGenotyping)

    variantSetsFromGatk\
        = removeIntervalList(variantSetAndIntervalListPairsFromGatk)

    (sampleVariantSetsFromStrelka)\
        = CallVariantsWithStrelka(\
            sampleReadGroups,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    //(sampleVariantSetsFromManta)\
    //    = CallVariantsWithManta(\
    //        sampleReadGroups,\
    //        referenceSequenceFasta,\
    //        referenceSequenceIndex)

    (sampleVariantSetsFromTiddit,\
     tidditOut)\
        = CallVariantsWithTiddit(\
            sampleReadGroups,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    (variantSetsFromFreebayes)\
        = CallVariantsWithFreebayes(\
            sampleReadGroupAndIntervalListPairs,\
            referenceSequenceFasta,\
            softwareVersions)

    //  concatenate variant sets across intervals  //
    //  -- only valid for freebayes and gatk  //

    sampleGroupsOfVariantSetsFromGatk\
        = groupByPatientSample(variantSetsFromGatk)
    sampleGroupsOfVariantSetsFromFreeBayes\
        = groupByPatientSample(variantSetsFromFreebayes)

    sampleGroupsOfVariantSets\
        = sampleGroupsOfVariantSetsFromGatk.mix(
            sampleGroupsOfVariantSetsFromFreeBayes)

    sampleVariantSetsFromGatkAndFreebayes\
        = MergeVariantSetsForSample(
            sampleGroupsOfVariantSets,
            referenceSequenceIndex)

    sampleVariantSets\
        = sampleVariantSetsFromGatkAndFreebayes.mix(
            sampleVariantSetsFromStrelka,
            //sampleVariantSetsFromManta,
            sampleVariantSetsFromTiddit)

    writeInputsForNextStep(
        sampleVariantSets,
        __genderMap__,
        __statusMap__)

    /**/

}
