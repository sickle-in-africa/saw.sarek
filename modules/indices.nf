include {
    isChannelActive;
    getInactiveChannelFlag
} from "${params.modulesDir}/sarek.nf"


process GetBwaIndexes {
    // if the BWA indexes are available on igenomes, 
    // pull them, otherwise build them locally

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/BWAIndex/${it}" : null }

    input:
        path(fasta)
        path(_bwa_)

    output:
        path("${fasta}.*"), includeInputs: true

    script:
        aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
        if ( isChannelActive(_bwa_) == true )
            """
            : # do nothing. Indexes will get automatically pulled from iGenomes
            """
        else
            """
            ${aligner} index ${fasta}
            """
}

process GetGatkDictionary {
    // if the gatk dictionary is available on igenomes,
    // pull it, otherwise build it locally
    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(fasta)
        path(_dict_)

    output:
        path("${fasta.baseName}.dict"), includeInputs: true

    script:
        if ( isChannelActive(_dict_) == true )
            """
            : # do nothing. Dictionary will get automatically pulled from iGenomes
            """
        else
            """
            gatk --java-options "-Xmx${task.memory.toGiga()}g" \
                CreateSequenceDictionary \
                --REFERENCE ${fasta} \
                --OUTPUT ${fasta.baseName}.dict
            """

}

process GetReferenceSequenceIndex {
    // if the samtools fasta index (fai) is available on iGenomes,
    // pull it, otherwise build it locally

    tag "${fasta}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(fasta)
        path(_fastaFai_)

    output:
        path("${fasta}.fai"), includeInputs: true

    script:
        if ( isChannelActive(_fastaFai_) == true )
            """
            : do nothing. Fasta index will be pulled from iGenomes
            """
        else
            """
            samtools faidx ${fasta}
            """
}

process GetDbsnpIndex {
    // if the user has specified an input dbsnp file, 
    // build the index, otherwise, emit an inactive channel

    tag "${dbsnp}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        path(_dbsnp_)

    output:
        path("${outputPath}")

    script:
        outputPath = "${_dbsnp_}.tbi"
        if ( isChannelActive(_dbsnp_) == true )
            """
            tabix -p vcf ${_dbsnp_}
            """
        else
            outputPath = getInactiveChannelFlag()
            """
            touch ${outputPath}
            """
}

process GetKnownIndelsIndex {
    tag "${knownIndels}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {params.save_reference ? "reference_genome/${it}" : null }

    input:
        each path(_knownIndels_)

    output:
        path("${outputPath}")

    script:
        outputPath = "${_knownIndels_}.tbi"
        if ( isChannelActive(_knownIndels_) == true )
            """
            tabix -p vcf ${_knownIndels_}
            """
        else
            outputPath = getInactiveChannelFlag()
            """
            touch ${outputPath}
            """
}