process GetVariantRecalibrationReport {
    label 'cpus_1'
    //label 'withGatkContainer'

    tag "${variantCaller}-${idSample}"

    input:
        tuple val(variantCaller), val(idSample), file(vcf)
        file(fasta)
        file(dict)
        file(fastaFai)
        file(dbsnp)
        file(dbsnpIndex)
        path hapmap
        path hapmap_index
        path onekgSnps
        path onekgSnpsIndex
        path onekgIndels
        path onekgIndelsIndex
        path onekgOmni
        path onekgOmniIndex


    output:
        path "${variantCaller}.${idSample}.recal"

    script:
        """
        gatk VariantRecalibrator \
            -R ${fasta} \
            -V ${vcf} \
            --resource:hapmap,known=false,training=true,truth=true,prior=15.0 ${hapmap} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekgSnps} \
            --resource:1000G,known=false,training=true,truth=false,prior=10.0 ${onekgIndels} \
            --resource:omni,known=false,training=true,truth=false,prior=12.0 ${onekgOmni} \
            --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 ${dbsnp} \
            -an QD \
            -an MQ \
            -an MQRankSum \
            -an ReadPosRankSum \
            -an FS \
            -an SOR \
            -an InbreedingCoeff \
            -mode BOTH \
            -O ${variantCaller}.${idSample}.recal \
            --tranches-file ${variantCaller}.${idSample}.tranches \
            --rscript-file ${variantCaller}.${idSample}.plots.R
        """

}

process RecalibrateVariantQualityScores {
    label 'cpus_1'
    label 'withGatkContainer'

    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/recalibratedVariants/${idSample}/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple val(variantCaller), val(idSample), file(vcf)
        path recalibrationReport

    script:
        """
        gatk ApplyVQSR \
            -V ${vcf} \
            --recal-file ${recalibrationReport} \
            -O ${variantCaller}_${idSample}.vcf.gz
        """
}