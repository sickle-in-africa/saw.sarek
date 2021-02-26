#!/usr/bin/env nextflow
/*
 *  VARIANT CALLING
 *  ===============
 *
 *      This script takes in *aligned sample read groups* (single files of reads
 *      such that one read file corresponds to one and only one patient sample)
 *      in the bam format, creates intervals (one interval per interval/bed 
 *      file) and then calls variants in the sample read groups for a range of 
 *      variant callers. 
 *      
 *      The callers are either for snvs and small indels:
 *          + gatk 4 Haplotype caller
 *          + strelka
 *          + freebayes
 *      or for structural variants:
 *          + manta
 *          + tiddit
 *
 *      Calling variants on intervals is helpful for speeding up the process, by
 *      parralellising over each interval, however it is only possible for 
 *      the gatk and freebayes callers. Once variant calling is done, all the 
 *      variant sets corresponding to each sample are merged across intervals 
 *      (where intervals were used) and the sample variant sets are output.  
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
    GetBwaIndexes;
    GetGatkDictionary;
    GetReferenceSequenceIndex;
    GetDbsnpIndex;
    GetKnownIndelsIndex
} from "${params.modulesDir}/indices.nf"

include {
    GetIntervalsPlan;
    GetIntervals;
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
    removeIntervals;
} from "${params.modulesDir}/calling.nf"


workflow {

    (sampleReadGroups,\
     igenomesReferenceSequenceFasta,
     igenomesReferenceSequenceDictionary,
     igenomesReferenceSequenceIndex,
     dbsnp,
     intervalsPlanFromInput,
     targetIntervalFromInput,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForCalling()

    softwareVersions = GetSoftwareVersions()

    //  get reference indexes as channels  //
    
    referenceSequenceDictionary\
        = GetGatkDictionary(\
            igenomesReferenceSequenceFasta,\
            igenomesReferenceSequenceDictionary)

    referenceSequenceIndex\
        = GetReferenceSequenceIndex(\
            igenomesReferenceSequenceFasta,\
            igenomesReferenceSequenceIndex)

    dbsnpIndex = GetDbsnpIndex(dbsnp)

   //  set up intervals for calling variants in parallel  //

    intervalsPlan\
        = GetIntervalsPlan(\
            referenceSequenceIndex,\
            intervalsPlanFromInput)

    intervals = GetIntervals(intervalsPlan).flatten()

    intervalsWithDurations\
        = addDurationToInterval(intervals)

    sampleReadGroupAndIntervalPairs\
        = sampleReadGroups.combine(intervalsWithDurations)

    //  call variants  //

    variantSetAndIntervalPairsFromGatk\
        = CallVariantsWithGatk(\
            sampleReadGroupAndIntervalPairs,\
            dbsnp,
            dbsnpIndex,\
            referenceSequenceDictionary,\
            igenomesReferenceSequenceFasta,\
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
            igenomesReferenceSequenceFasta,\
            referenceSequenceIndex)

    variantSetAndIntervalPairsFromGatk\
        = genotyped.mix(noGenotyping)

    variantSetsFromGatk\
        = removeIntervals(variantSetAndIntervalPairsFromGatk)    

    (sampleVariantSetsFromStrelka)\
        = CallVariantsWithStrelka(\
            sampleReadGroups,\
            igenomesReferenceSequenceFasta,\
            referenceSequenceIndex,\
            targetIntervalFromInput)

    (sampleVariantSetsFromManta)\
        = CallVariantsWithManta(\
            sampleReadGroups,\
            igenomesReferenceSequenceFasta,\
            referenceSequenceIndex,\
            targetIntervalFromInput)

    (sampleVariantSetsFromTiddit,\
     tidditOut)\
        = CallVariantsWithTiddit(\
            sampleReadGroups,\
            igenomesReferenceSequenceFasta,\
            referenceSequenceIndex)

    (variantSetsFromFreebayes)\
        = CallVariantsWithFreebayes(\
            sampleReadGroupAndIntervalPairs,\
            igenomesReferenceSequenceFasta,\
            softwareVersions)

    //  concatenate variant sets across intervals  //
    //  -- only valid for freebayes and gatk  //

    sampleGroupsOfVariantSetsFromGatk\
        = groupByPatientSample(variantSetsFromGatk)
    sampleGroupsOfVariantSetsFromFreeBayes\
        = groupByPatientSample(variantSetsFromFreebayes)

    sampleGroupsOfVariantSets\
        = sampleGroupsOfVariantSetsFromGatk.mix(\
            sampleGroupsOfVariantSetsFromFreeBayes)

    sampleVariantSetsFromGatkAndFreebayes\
        = MergeVariantSetsForSample(\
            sampleGroupsOfVariantSets,\
            referenceSequenceIndex,\
            targetIntervalFromInput)

    sampleVariantSets\
        = sampleVariantSetsFromGatkAndFreebayes.mix(\
            sampleVariantSetsFromStrelka,\
            sampleVariantSetsFromManta,
            sampleVariantSetsFromTiddit)

    /**/

}