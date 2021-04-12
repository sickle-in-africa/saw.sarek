include {
    reduceVCF
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
