include {
    isChannelActive
} from "${params.modulesDir}/sarek.nf"

process GetSoftwareVersions {
    publishDir "${params.outdir}/pipeline_info", mode: params.publish_dir_mode,
        saveAs: {it.indexOf(".csv") > 0 ? it : null}

    output:
        file 'software_versions_mqc.yaml'
        //file "software_versions.csv"

    when: !('versions' in skipQC)

    script:
    aligner = params.aligner == "bwa-mem2" ? "bwa-mem2" : "bwa"
    """
    alleleCounter --version &> v_allelecount.txt 2>&1 || true
    bcftools --version &> v_bcftools.txt 2>&1 || true
    ${aligner} version &> v_bwa.txt 2>&1 || true
    cnvkit.py version &> v_cnvkit.txt 2>&1 || true
    configManta.py --version &> v_manta.txt 2>&1 || true
    configureStrelkaGermlineWorkflow.py --version &> v_strelka.txt 2>&1 || true
    echo "${workflow.manifest.version}" &> v_pipeline.txt 2>&1 || true
    echo "${workflow.nextflow.version}" &> v_nextflow.txt 2>&1 || true
    snpEff -version &> v_snpeff.txt 2>&1 || true
    fastqc --version &> v_fastqc.txt 2>&1 || true
    freebayes --version &> v_freebayes.txt 2>&1 || true
    freec &> v_controlfreec.txt 2>&1 || true
    gatk ApplyBQSR --help &> v_gatk.txt 2>&1 || true
    msisensor &> v_msisensor.txt 2>&1 || true
    multiqc --version &> v_multiqc.txt 2>&1 || true
    qualimap --version &> v_qualimap.txt 2>&1 || true
    R --version &> v_r.txt 2>&1 || true
    R -e "library(ASCAT); help(package='ASCAT')" &> v_ascat.txt 2>&1 || true
    samtools --version &> v_samtools.txt 2>&1 || true
    tiddit &> v_tiddit.txt 2>&1 || true
    trim_galore -v &> v_trim_galore.txt 2>&1 || true
    vcftools --version &> v_vcftools.txt 2>&1 || true
    vep --help &> v_vep.txt 2>&1 || true

    ${params.sarekDir}/bin/scrape_software_versions.py &> software_versions_mqc.yaml
    """
}

process HaplotypeCaller {
    label 'memory_singleCPU_task_sq'
    label 'cpus_2'

    tag "${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(intervalBed)
        path(dbsnp)
        path(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
        tuple val("HaplotypeCallerGVCF"), val(idPatient), val(idSample), file("${intervalBed.baseName}_${idSample}.g.vcf")
        tuple val(idPatient), val(idSample), file(intervalBed), file("${intervalBed.baseName}_${idSample}.g.vcf")

    //when: 'haplotypecaller' in tools

    script:
    //javaOptions = "-Xmx${task.memory.toGiga()}g -Xms6000m -XX:GCTimeLimit=50 -XX:GCHeapFreeLimit=10"
    javaOptions = ""
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

process GenotypeGVCFs {
    tag "${idSample}-${intervalBed.baseName}"

    input:
        tuple val(idPatient), val(idSample), file(intervalBed), file(gvcf)
        file(dbsnp)
        file(dbsnpIndex)
        file(dict)
        file(fasta)
        file(fastaFai)

    output:
    tuple val("HaplotypeCaller"), val(idPatient), val(idSample), file("${intervalBed.baseName}_${idSample}.vcf")

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

process StrelkaSingle {
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
        tuple val("Strelka"), val(idPatient), val(idSample), file("*.vcf.gz"), file("*.vcf.gz.tbi")

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


process MantaSingle {
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
        tuple val("Manta"), val(idPatient), val(idSample), file("*.vcf.gz"), file("*.vcf.gz.tbi")

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

process TIDDIT {
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
        tuple val("TIDDIT"), val(idPatient), val(idSample), file("*.vcf.gz"), file("*.tbi")
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

process FreebayesSingle {
    tag "${idSample}-${intervalBed.baseName}"

    label 'cpus_1'
    
    input:
        tuple val(idPatient), val(idSample), file(bam), file(bai), file(intervalBed)
        file(fasta)
        file(fastaFai)
    
    output:
        tuple val("FreeBayes"), val(idPatient), val(idSample), file("${intervalBed.baseName}_${idSample}.vcf")
    
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

process ConcatVCF {
    label 'concat_vcf'
    label 'cpus_8'

    tag "${variantCaller}-${idSample}"

    publishDir "${params.outdir}/VariantCalling/${idSample}/${"$variantCaller"}", mode: params.publish_dir_mode

    input:
        tuple val(variantCaller), val(idPatient), val(idSample), file(vcf)
        file(fastaFai)
        file(targetBED)

    output:
    // we have this funny *_* pattern to avoid copying the raw calls to publishdir
        tuple val(variantCaller), val(idPatient), val(idSample), file("*_*.vcf.gz"), file("*_*.vcf.gz.tbi")

    //when: ('haplotypecaller' in tools || 'mutect2' in tools || 'freebayes' in tools)

    script:
    if (variantCaller == 'HaplotypeCallerGVCF')
        outputFile = "HaplotypeCaller_${idSample}.g.vcf"
    else
        outputFile = "${variantCaller}_${idSample}.vcf"
    options = params.target_bed ? "-t ${targetBED}" : ""
    intervalsOptions = params.no_intervals ? "-n" : ""
    """
    ${params.sarekDir}/bin/concatenateVCFs.sh -i ${fastaFai} -c ${task.cpus} -o ${outputFile} ${options} ${intervalsOptions}
    """
}

def groupVcfChannelsAcrossIntervalsAndMix(\
    gvcfHaplotypeCaller,\
    vcfGenotypeGVCFs,\
    vcfFreebayesSingle) {

    gvcfHaplotypeCallerGrouped = gvcfHaplotypeCaller.groupTuple(by:[0, 1, 2])
    vcfGenotypeGVCFsGrouped = vcfGenotypeGVCFs.groupTuple(by:[0, 1, 2])
    vcfFreebayesSingleGrouped = vcfFreebayesSingle.groupTuple(by: [0,1,2])

    return Channel
            .empty()
            .mix(vcfFreebayesSingleGrouped, vcfGenotypeGVCFsGrouped, gvcfHaplotypeCallerGrouped)

}