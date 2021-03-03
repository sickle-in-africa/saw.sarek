#!/usr/bin/env nextflow
/*
 *  READS PREPROCESSING AND ALIGNMENT
 *  =================================
 *
 *      This script takes in unaligned read groups,
 *      preprocesses them, aligns them to a reference, 
 *      and then marks duplicates. The output of this step
 *      can then be fed into the base recalibration step.
 *
 *      Here we define a *read group* as a single file that
 *      contains reads from a single patient sample. Note that
 *      multiple read groups may correspond to the same patient
 *      sample, if for example the same sample was sequenced
 *      over several runs. The input read groups can be in
 *      fastq pairs or (unmapped) bam format.  
 *
 *      There is an option to specify splitting fastq read groups
 *      across multiple files, which can be passed to the aligner
 *      in parallel to speed up the alignment. After alignment,
 *      the read groups are grouped by patient sample, (into groups
 *      of read groups) and these are merged across runs into
 *      *sample read groups*, where each read group (file of reads)
 *      corresponds to one and only one patient sample.
 *
 *****************************************************************/

nextflow.enable.dsl=2

include {
    groupByPatientSample
} from "${params.modulesDir}/sarek.nf"

include {
    initializeInputChannelsForMapping
} from "${params.modulesDir}/inputs.nf"

include {
    FastQCFQ;
    FastQCBAM;
    TrimReads;
    UMIFastqToBAM;
    UMIMapBamFile;
    GroupReadsByUmi;
    CallMolecularConsensusReads;
    branchReadGroupsIntoBamOrFastqChannels;
    splitReadGroups;
    branchReadGroupsIntoPreProcessingChannels
} from "${params.modulesDir}/preprocess.nf"

include {
    AlignReadsToReferenceSequence;
    MergeReadGroupsForSample;
    GetIndexOfAlignedSampleReadGroup;
    MarkDuplicatesInSampleReadGroup;
    branchIntoSingleOrMultipleGroupChannels;
    writeTsvFilesForBams;
    writeTsvFilesForBamsWithDuplicatesMarked
} from "${params.modulesDir}/alignment.nf"


workflow {

    (readGroupsFromInput,\
     referenceSequenceFasta,
     bwaIndexTuple,
     referenceSequenceIndex,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForMapping()

    //  preprocess reads  //

    (readGroupsAsBam,\
     readGroupsAsFastq)\
        = branchReadGroupsIntoBamOrFastqChannels(\
            readGroupsFromInput)

    readGroupsAsFastqSplit = splitReadGroups(readGroupsAsFastq)

    readGroups = readGroupsAsBam.mix(readGroupsAsFastqSplit)

    (readGroupsForTrimming,\
     readGroupsForUmiProcessing,\
     readGroupsNoPreProcessing)\
        = branchReadGroupsIntoPreProcessingChannels(\
            readGroups)

    (trimGaloreReport,\
     readGroupsTrimmed)
        = TrimReads(\
            readGroupsForTrimming)

    umi_converted_bams_ch\
        = UMIFastqToBAM(\
            readGroupsForUmiProcessing)

    umi_aligned_bams_ch\
        = UMIMapBamFile(\
            umi_converted_bams_ch,\
            bwaIndexTuple,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    (umi_histogram_ch,\
     umi_grouped_bams_ch)\
        = GroupReadsByUmi(\
            umi_aligned_bams_ch)

    readGroupsUmiProcessed\
        = CallMolecularConsensusReads(\
            umi_grouped_bams_ch)

    //  map reads to reference with bwa  //

    readGroupsForAligning\
        = readGroupsTrimmed.mix(\
            readGroupsUmiProcessed,\
            readGroupsNoPreProcessing)

    readGroupsAligned\
        = AlignReadsToReferenceSequence(\
            readGroupsForAligning,\
            bwaIndexTuple,\
            referenceSequenceFasta,\
            referenceSequenceIndex)

    sampleGroupsOfAlignedReadGroups\
        = groupByPatientSample(readGroupsAligned)

    (singleGroups,\
     multipleGroups)\
        = branchIntoSingleOrMultipleGroupChannels(\
        sampleGroupsOfAlignedReadGroups)

    sampleReadGroupsMerged\
        = MergeReadGroupsForSample(multipleGroups)

    sampleReadGroupsAligned = singleGroups.mix(sampleReadGroupsMerged)

    (sampleReadGroupsAlignedIndex,\
     tsv_bam_indexed)\
        = GetIndexOfAlignedSampleReadGroup(sampleReadGroupsAligned)

    //  mark duplicate reads  //

    (sampleReadGroupsMarked,\
     tsv_bam_duplicates_marked,\
     duplicates_marked_report)\
        = MarkDuplicatesInSampleReadGroup(sampleReadGroupsAligned)

    //  write tsv files for input to the next step  //

    writeTsvFilesForBams(\
        tsv_bam_indexed,\
        __genderMap__,\
        __statusMap__)

    writeTsvFilesForBamsWithDuplicatesMarked(\
        tsv_bam_duplicates_marked,\
        __genderMap__,\
        __statusMap__)

    /*

    inputBamFastQC = stripSecondInputFile(readGroupsAsBam)
    fastQCFQReport = FastQCFQ(readGroups)
    fastQCBAMReport = FastQCBAM(inputBamFastQC)
    fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

    /**/

}