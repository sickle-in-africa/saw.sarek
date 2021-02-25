#!/usr/bin/env nextflow
nextflow.enable.dsl=2

include {
    initializeInputChannelsForCalling
} from "${params.modulesDir}/inputs.nf"

include {
    GetBwaIndexes;
    GetGatkDictionary;
    GetSamtoolsFastaIndex;
    GetDbsnpIndex;
    GetKnownIndelsIndex;
    GetIntervalsList
} from "${params.modulesDir}/indices.nf"

include {
    CreateIntervalBeds;
    addIntervalDurationsToIntervalsChannel;
} from "${params.modulesDir}/preprocess.nf"

include {
    GetSoftwareVersions;
    HaplotypeCaller;
    GenotypeGVCFs;
    StrelkaSingle;
    MantaSingle;
    TIDDIT;
    FreebayesSingle;
    ConcatVCF;
    groupVcfChannelsAcrossIntervalsAndMix
} from "${params.modulesDir}/calling.nf"


workflow {

    (inputSample,\
     _fasta_,
     _dict_,
     _fastaFai_,
     _dbsnp_,
     _intervalsList_,
     _targetBed_,
     __genderMap__,
     __statusMap__)\
        = initializeInputChannelsForCalling()

    //  get reference indexes as channels  //
    
    ch_dict = GetGatkDictionary(_fasta_, _dict_)

    ch_fai = GetSamtoolsFastaIndex(_fasta_, _fastaFai_)

    ch_dbsnp_tbi = GetDbsnpIndex(_dbsnp_)

    ch_intervals = GetIntervalsList(ch_fai, _intervalsList_)

    ch_software_versions_yaml = GetSoftwareVersions()

   //  preprocess reads  //

    bedIntervals = CreateIntervalBeds(ch_intervals).flatten()

    bedIntervalsWithDurations\
        = addIntervalDurationsToIntervalsChannel(bedIntervals)

    //  call variants  //

    inputSampleWithIntervals = inputSample.combine(bedIntervalsWithDurations)

    (gvcfHaplotypeCaller,\
     gvcfGenotypeGVCFs)\
        = HaplotypeCaller(\
            inputSampleWithIntervals,\
            _dbsnp_,
            ch_dbsnp_tbi,\
            ch_dict,\
            _fasta_,\
            ch_fai)

    (vcfGenotypeGVCFs)\
        = GenotypeGVCFs(\
            gvcfGenotypeGVCFs,\
            _dbsnp_,\
            ch_dbsnp_tbi,\
            ch_dict,\
            _fasta_,\
            ch_fai)

    (vcfStrelkaSingle)\
        = StrelkaSingle(\
            inputSample,\
            _fasta_,\
            ch_fai,\
            _targetBed_)

    (vcfMantaSingle)\
        = MantaSingle(\
            inputSample,\
            _fasta_,\
            ch_fai,\
            _targetBed_)

    (vcfTIDDIT,\
     tidditOut)\
        = TIDDIT(\
            inputSample,\
            _fasta_,\
            ch_fai)

    (vcfFreebayesSingle)\
        = FreebayesSingle(\
            inputSampleWithIntervals,\
            _fasta_,\
            ch_software_versions_yaml)

    //  concatenate variant sets across intervals  //

    groupedVcfs\
        = groupVcfChannelsAcrossIntervalsAndMix(\
            gvcfHaplotypeCaller,\
            vcfGenotypeGVCFs,\
            vcfFreebayesSingle)

    (vcfConcatenated)\
        = ConcatVCF(\
            groupedVcfs,\
            ch_fai,\
            _targetBed_)

    /**/

}