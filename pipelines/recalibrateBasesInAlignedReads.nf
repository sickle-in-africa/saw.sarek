#!/usr/bin/env nextflow
/*
 *  RECALIBRATE BASE CALLS IN ALIGNED READ GROUPS
 *  =============================================
 *
 *      This script uses the base quality score recalibration reports generated
 *      in the previous step to recalibrate the base call quality scores in the
 *      aligned sample read groups (.bam files). To speed up the process, we 
 *      first create intervals and then the sample read groups are recalibrated
 *      in parallel across each interval. After generating the reports, we merge
 *      them across intervals for each sample. 
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
 *      recalibration is done on these interval lists in parallel. 
 *
 *****************************************************************************/

nextflow.enable.dsl=2

include {
    groupByPatientSample
} from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForRecalibration
} from "${params.modulesDir}/inputs.nf"

include {
    BuildReferenceIntervalList;
    SplitIntervalList;
    addDurationToInterval
} from "${params.modulesDir}/intervals.nf"

include {
    RecalibrateBasesInReadGroup;
    MergeRecalibratedReadGroupsForSample;
    writeTsvFilesForRecalibratedBams
} from "${params.modulesDir}/baseRecalibration.nf"


workflow {

    (sampleReadGroups,
     referenceSequenceFasta,
     referenceSequenceDictionary,
     referenceSequenceIndex,
     dbsnp,
     dbsnpIndex,
     referenceIntervalListFromInput,
     knownIndels,
     knownIndelsIndex,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForRecalibration()

    //  set up intervals for recalibrating bases in parallel  //

    referenceIntervalListFromIndex\
        = BuildReferenceIntervalList(\
            referenceSequenceIndex,\
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

    //  recalibrate aligned sample reads in parallel over intervals  //

    readGroupsRecalibrated\
        = RecalibrateBasesInReadGroup(\
            sampleReadGroupAndIntervalListPairs,\
            referenceSequenceDictionary,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    //  group and merge recalibrated bam files for each sample  //

    sampleGroupsOfReadGroupsRecalibrated\
        = groupByPatientSample(\
            readGroupsRecalibrated)

    (sampleReadGroupsRecalibrated,\
     bam_recalibrated_qc,\
     tsv_bam_recalibrated)\
        = MergeRecalibratedReadGroupsForSample(\
            sampleGroupsOfReadGroupsRecalibrated)

    //  write tsv files for the next step  //

    writeTsvFilesForRecalibratedBams(\
        tsv_bam_recalibrated,\
        __genderMap__,\
        __statusMap__)

    /**/

}
