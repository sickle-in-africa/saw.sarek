process FilterVariantsFromFreebayes {

    publishDir "${params.outdir}/Filtered/${idSample}/Freebayes", mode: 'copy'

    input:
        tuple val(variantCaller), val(idSample), path(vcf)

    script:
        """
        vcffilter -f "QUAL > 20" ${vcf} > ${vcf}
        """

}

def branchIntoGatkStrelkaOrFreebayesChannels(variantSets) {
    result = variantSets
        .branch {
            gatk: it[0] == 'HaplotypeCaller'
            strelka: it[0] == 'Strelka'
            freebayes: it[0] == 'FreeBayes'
        }

    return [
        result.gatk,
        result.strelka,
        result.freebayes]
}
