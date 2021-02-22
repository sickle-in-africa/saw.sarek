#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    initializeInputChannelsForMapping
} from "${params.modulesDir}/inputs.nf"

include {
    GetBwaIndexes;
    GetGatkDictionary;
    GetSamtoolsFastaIndex;
    GetDbsnpIndex;
    GetKnownIndelsIndex
} from "${params.modulesDir}/indices.nf"

include {
    FastQCFQ;
    FastQCBAM;
    TrimGalore;
    UMIFastqToBAM;
    UMIMapBamFile;
    GroupReadsByUmi;
    CallMolecularConsensusReads;
    splitInputsIntoBamAndFastaPairs;
    stripSecondInputFile;
    splitFastqFiles;
} from "${params.modulesDir}/preprocess.nf"

include {
    MapReads;
    MergeBamMapped;
    IndexBamFile;
    MarkDuplicates;
    selectPairReadsChannelForMapping;
    splitMappedBamsIntoSingleAndMultipeLanes;
    writeTsvFilesForBams;
    writeTsvFilesForBamsWithDuplicatesMarked
} from "${params.modulesDir}/alignment.nf"

workflow {

    (inputSample,\
     _fasta_,
     _bwa_,
     _dict_,
     _fastaFai_,
     _dbsnp_,
     _knownIndels_,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForMapping()

    //  get reference indexes as channels  //

    ch_bwa = GetBwaIndexes(_fasta_, _bwa_)
    
    ch_dict = GetGatkDictionary(_fasta_, _dict_)

    ch_fai = GetSamtoolsFastaIndex(_fasta_, _fastaFai_)

    ch_dbsnp_tbi = GetDbsnpIndex(_dbsnp_)

    ch_known_indels_tbi = GetKnownIndelsIndex(_knownIndels_)


    //  preprocess reads  //

    (inputBam,\
     inputPairReads)\
        = splitInputsIntoBamAndFastaPairs(inputSample)

    inputBamFastQC = stripSecondInputFile(inputBam)

    inputPairReadsSplit\
        = splitFastqFiles(\
            inputPairReads) \
        | mix(inputBam)

    fastQCFQReport = FastQCFQ(inputPairReadsSplit)

    fastQCBAMReport = FastQCBAM(inputBamFastQC)

    fastQCReport = fastQCFQReport.mix(fastQCBAMReport)

    (trimGaloreReport,\
     outputPairReadsTrimGalore)
        = TrimGalore(\
            inputPairReadsSplit)

    //  UMIs processing (optional, default: off)  //

    umi_converted_bams_ch\
        = UMIFastqToBAM(\
            inputPairReadsSplit)

    umi_aligned_bams_ch\
        = UMIMapBamFile(\
            umi_converted_bams_ch,\
            ch_bwa,\
            _fasta_,\
            ch_fai)

    (umi_histogram_ch,\
     umi_grouped_bams_ch)\
        = GroupReadsByUmi(\
            umi_aligned_bams_ch)

    consensus_bam_ch\
        = CallMolecularConsensusReads(\
            umi_grouped_bams_ch)

    //  map reads to reference with bwa  //

    preprocessedPairReads\
        = selectPairReadsChannelForMapping(\
            consensus_bam_ch,\
            outputPairReadsTrimGalore,\
            inputPairReadsSplit) 

    (bamMapped,\
     bamMappedBamQC)\
        = MapReads(\
            preprocessedPairReads,\
            ch_bwa,\
            _fasta_,\
            ch_fai)

    (singleBam,\
     multipleBam)\
        = splitMappedBamsIntoSingleAndMultipeLanes(bamMapped)

    bam_mapped_merged\
        = MergeBamMapped(multipleBam)
            .mix(singleBam)

    (bam_mapped_merged_indexed,\
     tsv_bam_indexed)\
        = IndexBamFile(bam_mapped_merged)

    (bam_duplicates_marked,\
     tsv_bam_duplicates_marked,\
     duplicates_marked_report)\
        = MarkDuplicates(bam_mapped_merged)

    writeTsvFilesForBams(\
        tsv_bam_indexed,\
        __genderMap__,\
        __statusMap__)

    writeTsvFilesForBamsWithDuplicatesMarked(\
        tsv_bam_duplicates_marked,\
        __genderMap__,\
        __statusMap__)

    /**/

}