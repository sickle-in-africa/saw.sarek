include {
    hasExtension;
    isChannelActive
} from "${params.modulesDir}/sarek.nf"


process CallVariantsWithGatk {
    label 'memory_singleCPU_task_sq'
    label 'cpus_2'
    label 'withGatkContainer'

    tag "${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(intervalBed)
        path(dbsnp)
        path(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
        tuple val(idPatient), val(idSample), val("HaplotypeCallerGVCF"), file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf")

    //when: 'haplotypecaller' in tools

    script:
    javaOptions = (params.profile == 'chpc') ? "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10" : ""
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = isChannelActive(dbsnp) ? "--D ${dbsnp}" : ""
    """
    gatk --java-options "${javaOptions}" \
        HaplotypeCaller \
        -R ${fasta} \
        -I ${bam} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -O ${intervalBed.baseName}_${idSample}.g.vcf \
        -ERC GVCF
    """
}

process GenotypeVariantsFromGatk {
    label 'withGatkContainer'

    tag "${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), val(variantCaller), file(intervalBed), file(gvcf)
        file(dbsnp)
        file(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
    tuple val(idPatient), val(idSample), val("HaplotypeCaller"), file(intervalBed), file("${intervalBed.baseName}_${idSample}.vcf")

    //when: 'haplotypecaller' in tools

    script:
    // Using -L is important for speed and we have to index the interval files also
    intervalsOptions = params.no_intervals ? "" : "-L ${intervalBed}"
    dbsnpOptions = isChannelActive(dbsnp) ? "--D ${dbsnp}" : ""
    """
    gatk --java-options -Xmx${task.memory.toGiga()}g \
        IndexFeatureFile \
        -I ${gvcf}

    gatk --java-options -Xmx${task.memory.toGiga()}g \
        GenotypeGVCFs \
        -R ${fasta} \
        ${intervalsOptions} \
        ${dbsnpOptions} \
        -V ${gvcf} \
        -O ${intervalBed.baseName}_${idSample}.vcf
    """
}

process CallVariantsWithStrelka {
    label 'cpus_max'
    label 'memory_max'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Strelka", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        tuple val(idPatient), val(idSample), val("Strelka"), file("*.vcf.gz"), file("*.vcf.gz.tbi")

    //when: 'strelka' in tools

    script:
    beforeScript = isChannelActive(targetBED) ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = isChannelActive(targetBED) ? "--exome --callRegions call_targets.bed.gz" : ""
    """
    ${beforeScript}
    configureStrelkaGermlineWorkflow.py \
        --bam ${bam} \
        --referenceFasta ${fasta} \
        ${options} \
        --runDir Strelka

    python Strelka/runWorkflow.py -m local -j ${task.cpus}

    mv Strelka/results/variants/genome.*.vcf.gz \
        Strelka_${idSample}_genome.vcf.gz
    mv Strelka/results/variants/genome.*.vcf.gz.tbi \
        Strelka_${idSample}_genome.vcf.gz.tbi
    mv Strelka/results/variants/variants.vcf.gz \
        Strelka_${idSample}_variants.vcf.gz
    mv Strelka/results/variants/variants.vcf.gz.tbi \
        Strelka_${idSample}_variants.vcf.gz.tbi
    """
}


process CallVariantsWithManta {
    label 'cpus_max'
    label 'memory_max'

    tag "${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/Manta", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai)
        file(fasta)
        file(fastaFai)
        file(targetBED)

    output:
        tuple val(idPatient), val(idSample), val("Manta"), file("*.vcf.gz"), file("*.vcf.gz.tbi")

    //when: 'manta' in tools

    script:
    beforeScript = params.target_bed ? "bgzip --threads ${task.cpus} -c ${targetBED} > call_targets.bed.gz ; tabix call_targets.bed.gz" : ""
    options = params.target_bed ? "--exome --callRegions call_targets.bed.gz" : ""
    inputbam = "--bam" 
    vcftype = "diploid"
    """
    ${beforeScript}
    configManta.py \
        ${inputbam} ${bam} \
        --reference ${fasta} \
        ${options} \
        --runDir Manta

    python Manta/runWorkflow.py -m local -j ${task.cpus}

    mv Manta/results/variants/candidateSmallIndels.vcf.gz \
        Manta_${idSample}.candidateSmallIndels.vcf.gz
    mv Manta/results/variants/candidateSmallIndels.vcf.gz.tbi \
        Manta_${idSample}.candidateSmallIndels.vcf.gz.tbi
    mv Manta/results/variants/candidateSV.vcf.gz \
        Manta_${idSample}.candidateSV.vcf.gz
    mv Manta/results/variants/candidateSV.vcf.gz.tbi \
        Manta_${idSample}.candidateSV.vcf.gz.tbi
    mv Manta/results/variants/${vcftype}SV.vcf.gz \
        Manta_${idSample}.${vcftype}SV.vcf.gz
    mv Manta/results/variants/${vcftype}SV.vcf.gz.tbi \
        Manta_${idSample}.${vcftype}SV.vcf.gz.tbi
    """
}

process CallVariantsWithTiddit {
    tag "${idSample}"

    publishDir params.outdir, mode: params.publish_dir_mode,
        saveAs: {
            if (it == "TIDDIT_${idSample}.vcf") "VariantCalling/${idSample}/TIDDIT/${it}"
            else "Reports/${idSample}/TIDDIT/${it}"
        }

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai)
        file(fasta)
        file(fastaFai)

    output:
        tuple val(idPatient), val(idSample), val("TIDDIT"), file("*.vcf.gz"), file("*.tbi")
        tuple file("TIDDIT_${idSample}.old.vcf"), file("TIDDIT_${idSample}.ploidy.tab"), file("TIDDIT_${idSample}.signals.tab"), file("TIDDIT_${idSample}.wig"), file("TIDDIT_${idSample}.gc.wig")

    //when: 'tiddit' in tools

    script:
    """
    tiddit --sv -o TIDDIT_${idSample} --bam ${bam} --ref ${fasta}

    mv TIDDIT_${idSample}.vcf TIDDIT_${idSample}.old.vcf

    grep -E "#|PASS" TIDDIT_${idSample}.old.vcf > TIDDIT_${idSample}.vcf

    bgzip --threads ${task.cpus} -c TIDDIT_${idSample}.vcf > TIDDIT_${idSample}.vcf.gz

    tabix TIDDIT_${idSample}.vcf.gz
    """
}

process CallVariantsWithFreebayes {
    tag "${idSample}-${intervalBed.baseName}"

    label 'cpus_1'
    
    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(intervalBed)
        file(fasta)
        file(fastaFai)
    
    output:
        tuple val(idPatient), val(idSample), val("FreeBayes"), file("${intervalBed.baseName}_${idSample}.vcf")
    
    //when: 'freebayes' in tools

    script:
    intervalsOptions = "-t ${intervalBed}"
    """
    freebayes \
        -f ${fasta} \
        --min-alternate-fraction 0.1 \
        --min-mapping-quality 1 \
        ${intervalsOptions} \
        ${bam} > ${intervalBed.baseName}_${idSample}.vcf
    """
}

process MergeVariantSetsForSample {
    label 'concat_vcf'
    label 'cpus_8'

    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

    input:
        tuple val(idPatient), val(idSample), val(variantCaller), file(vcf)
        file(fastaFai)
        file(targetBED)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        tuple val(idPatient), val(idSample), val(variantCaller), file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi")

    //when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    variantCaller = variantCaller[0]
    outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    intervalsOptions = params.no_intervals ? "-n" : ""
    """
    ${params.sarekDir}/bin/concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
    """
}

def branchIntoGenotypingOrNoGenotypingChannels(variantSets) {
    result = variantSets.branch {
        genotype: params.genotypeGatkVariants
        noGenotype: !params.genotypeGatkVariants
    }
    return [result.genotype, result.noGenotype]
}

def removeIntervalList(variantSetAndIntervalPairs) {
    return variantSetAndIntervalPairs.map {
        idPatient, idSample, variantCaller, intervalBed, variantSet\
        -> [idPatient, idSample, variantCaller, variantSet]
    }
}