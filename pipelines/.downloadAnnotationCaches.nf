#!/usr/bin/env nextflow
/*
 *  DOWNLOAD ANNOATION CACHES
 *  =========================
 *
 *  This script is for administrators only. Before trying to use this 
 *  pipeline yourself, first check that the annotation caches are not
 *  already saved in the public data space you are accessing. 
 *
 *  Caches are databases used by the annotation tools:
 *      + SnpEff
 *      + VEP
 *      + cadd
 *  They are downloaded before runtime and saved locally so we can save
 *  on downloading data. 
 *
 **********************************************************************/
nextflow.enable.dsl = 2

include {
    getInactiveValueChannel
} from "${params.modulesDir}/sarek.nf"

include {
    DownloadAnnotationCacheForVep
} from "${params.modulesDir}/annotation.nf"

workflow {

    params.vep_cache = params.genome ? params.genomes[params.genome].vep_cache : null
    params.vep_cache_version = params.genome ? params.genomes[params.genome].vep_cache_version : null

    ch_vep_cache = params.vep_cache ? Channel.value(params.vep_cache) : getInactiveValueChannel()
    ch_vep_cache_version = params.vep_cache_version ? Channel.value(params.vep_cache_version) : getInactiveValueChannel()
    ch_species = Channel.value(params.species)

    vep_cache_out\
        = DownloadAnnotationCacheForVep(
            ch_vep_cache_version,
            ch_species)
}