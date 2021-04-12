include {
    reduceVCF;
    isChannelActive
} from "${params.modulesDir}/sarek.nf" 

process AnnotateVariantsWithSnpeff {
    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${it}"
    }

    input:
        tuple val(variantCaller), val(idSample), file(vcf)
        file(snpEff_config)
        file(dataDir)
        val snpeffDb

    output:
        tuple file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv")
        tuple val(variantCaller), val(idSample), file("${reducedVCF}_snpEff.ann.vcf")

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        -dataDir \${PWD}/${dataDir} \
        -canon \
        -v ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    """
}

process CompressVariantSetFromSnpeff {
    tag "${idSample} - ${vcf}"

    publishDir "${params.outdir}/Annotation/${idSample}/snpEff", mode: params.publish_dir_mode

    input:
        tuple val(variantCaller), val(idSample), file(vcf)

    output:
        tuple val(variantCaller), val(idSample), file("*.vcf.gz"), file("*.vcf.gz.tbi")

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

process AnnotateVariantsWithVep {
    //label 'VEP'
    label 'cpus_4'

    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        tuple val(variantCaller), val(idSample), file(vcf)
        file(dataDir)
        val cache_version
        file(cadd_InDels)
        file(cadd_InDels_tbi)
        file(cadd_WG_SNVs)
        file(cadd_WG_SNVs_tbi)

    output:
        tuple val(variantCaller), val(idSample), file("${reducedVCF}_VEP.ann.vcf")
        file("${reducedVCF}_VEP.summary.html")

    script:
    species = params.genomes[params.genome].species
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = " \${PWD}/${dataDir}"
    cadd = (isChannelActive(cadd_InDels) && isChannelActive(cadd_WG_SNVs)) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf

    rm -rf ${reducedVCF}
    """
}

process MergeVariantSetsFromVepAndSnpeff {
    //label 'VEP'
    label 'cpus_4'

    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        tuple val(variantCaller), val(idSample), file(vcf)
        file(dataDir)
        val cache_version
        file(cadd_InDels)
        file(cadd_InDels_tbi)
        file(cadd_WG_SNVs)
        file(cadd_WG_SNVs_tbi)

    output:
        tuple val(variantCaller), val(idSample), file("${reducedVCF}_VEP.ann.vcf")
        file("${reducedVCF}_VEP.summary.html")

    script:
    species = params.genomes[params.genome].species
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = " \${PWD}/${dataDir}"
    cadd = (isChannelActive(cadd_InDels) && isChannelActive(cadd_WG_SNVs)) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${species} \
        ${cadd} \
        ${genesplicer} \
        --cache \
        --cache_version ${cache_version} \
        --dir_cache ${dir_cache} \
        --everything \
        --filter_common \
        --fork ${task.cpus} \
        --format vcf \
        --per_gene \
        --stats_file ${reducedVCF}_VEP.summary.html \
        --total_length \
        --vcf

    rm -rf ${reducedVCF}
    """
}

process CompressVariantSetFromVep {
    tag "${idSample} - ${vcf}"

    publishDir "${params.outdir}/Annotation/${idSample}/VEP", mode: params.publish_dir_mode

    input:
        tuple val(variantCaller), val(idSample), file(vcf)

    output:
        tuple val(variantCaller), val(idSample), file("*.vcf.gz"), file("*.vcf.gz.tbi")

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

