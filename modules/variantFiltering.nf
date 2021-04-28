include {
    reduceVCF
} from "${params.modulesDir}/sarek.nf"

process FilterVariantsFromFreebayes {

    publishDir "${params.outdir}/Filtered/${idSample}/Freebayes", mode: 'copy'

    input:
        tuple val(variantCaller), val(idSample), path(vcf)

    output:
        path (outputVcf + '.gz')

    script:
        outputVcf = reduceVCF(vcf) + '.filtered.vcf'
        """
        bcftools view -i'QUAL>20 && DP>10' ${vcf} > ${outputVcf}
        gzip ${outputVcf}
        """
}

process FilterVariantsFromStrelka {

    publishDir "${params.outdir}/Filtered/${idSample}/Strelka", mode: 'copy'

    input:
    tuple val(variantCaller), val(idSample), path(vcf)

    script:
    outputVcf = reduceVCF(vcf) + '.filtered.vcf'
        """
        bcftools view -i'FILTER="PASS"' ${vcf} > ${outputVcf}
        gzip ${outputVcf}
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
