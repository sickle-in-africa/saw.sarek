process FilterVariantsFromFreebayes {

    publishDir "${params.outdir}/Filtered/${idSample}/Freebayes", mode: 'copy'

    input:
        tuple val(variantCaller), val(idSample), path(vcf)

    script:
        """
        bcftools view -i'QUAL>20 && DP>10' > ${vcf}
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
