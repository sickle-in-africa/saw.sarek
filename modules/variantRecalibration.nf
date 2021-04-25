process PrintMessage {

    input:
    tuple val(variantCaller), val(idSample), path(vcf)

    output: 
        stdout

    script:
        """
        echo "${vcf}"
        """
}



process GetIndelRecalibrationReport {
    label 'cpus_1'
    tag "${variantCaller}-${idSample}"

    input:
        tuple val(variantCaller), val(idSample), path(vcf)
        path onekgIndels
        path onekgIndelsIndex
        path axiomExomePlus
        path axiomExomePlusIndex
        path dbsnp
        path dbsnpIndex

    output:
        tuple path("${variantCaller}_${idSample}_indels.recal"), path("${variantCaller}_${idSample}_indels.tranches")

    script:
        """
        gatk IndexFeatureFile \
            -I ${vcf}

        gatk --java-options "-Xmx${task.memory.toGiga()}g -Xms${task.memory.toGiga()}g" VariantRecalibrator \
            -V ${vcf} \
            --trust-all-polymorphic \
            -tranche 100.0 -tranche 99.95 -tranche 99.9 -tranche 99.5 -tranche 99.0 -tranche 97.0 -tranche 96.0 -tranche 95.0 -tranche 94.0 -tranche 93.5 -tranche 93.0 -tranche 92.0 -tranche 91.0 -tranche 90.0 \
            -an FS -an ReadPosRankSum -an MQRankSum -an QD -an SOR -an DP \
            -mode INDEL \
            --max-gaussians 4 \
            -resource:mills,known=false,training=true,truth=true,prior=12 ${onekgIndels} \
            -resource:axiomPoly,known=false,training=true,truth=false,prior=10 ${axiomExomePlus} \
            -resource:dbsnp,known=true,training=false,truth=false,prior=2 ${dbsnp} \
            -O ${variantCaller}_${idSample}_indels.recal \
            --tranches-file ${variantCaller}_${idSample}_indels.tranches
        """

}

process RecalibrateIndelQualityScores {
    label 'process_low'
    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/recalibratedVariants/${idSample}/${variantCaller}", mode: params.publish_dir_mode

    input:
        tuple val(variantCaller), val(idSample), path(vcf), path(vcfIndex)
        tuple path(recalTable), path(tranches)

    script:
        """
        gatk --java-options "-Xmx5g -Xms5g" \
            ApplyVQSR \
            -V ${vcf} \
            --recal-file ${recalTable} \
            --tranches-file ${tranches} \
            --truth-sensitivity-filter-level 99.7 \
            --create-output-variant-index true \
            -mode INDEL \
            -O ${idSample}_indel.recalibrated.vcf.gz
        """

}






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
        gatk IndexFeatureFile -I ${vcf}

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
