#!/usr/bin/env nextflow
/*
 *  RAW (UNFILTERED) VARIANT ANNOTATION
 *  ===================================
 *
 *
 *
 *******************************************/

nextflow.enable.dsl=2

include {
    initializeInputChannelsForAnnotation
} from "${params.modulesDir}/inputs.nf"

include {
    Snpeff;
    CompressVCFsnpEff;
    VEP;
    VEPmerge;
    CompressVCFvep
} from "${params.modulesDir}/variantAnnotation.nf"


workflow {

    (vcfSnpeff,
     vcfVep,
     ch_snpeff_cache,
     ch_snpeff_db,
     ch_vep_cache_version,
     ch_vep_cache,
     ch_cadd_indels,
     ch_cadd_indels_tbi,
     ch_cadd_wg_snvs,
     ch_cadd_wg_snvs_tbi)\
        = initializeInputChannelsForAnnotation()

    (snpeffReport, 
     snpeffVCF)\
        = Snpeff(
            vcfSnpeff,
            ch_snpeff_cache,
            ch_snpeff_db)

    (compressVCFsnpEffOut)\
        = CompressVCFsnpEff(snpeffVCF)

    (vepVCF,
     vepReport)\
        = VEP(
            vcfVep,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_indels,
            ch_cadd_indels_tbi,
            ch_cadd_wg_snvs,
            ch_cadd_wg_snvs_tbi)

    (vepVCFmerge,
     vepReportMerge)\
        = VEPmerge(
            compressVCFsnpEffOut,
            ch_vep_cache,
            ch_vep_cache_version,
            ch_cadd_indels,
            ch_cadd_indels_tbi,
            ch_cadd_wg_snvs,
            ch_cadd_wg_snvs_tbi)

    vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

    (compressVCFOutVEP) = CompressVCFvep(vcfCompressVCFvep)

}





if (step == 'annotate') {
    vcfToAnnotate = Channel.create()
    vcfNoAnnotate = Channel.create()

    if (tsvPath == []) {
    // Sarek, by default, annotates all available vcfs that it can find in the VariantCalling directory
    // Excluding vcfs from FreeBayes, and g.vcf from HaplotypeCaller
    // Basically it's: results/VariantCalling/*/{HaplotypeCaller,Manta,Mutect2,SentieonDNAseq,SentieonDNAscope,SentieonTNscope,Strelka,TIDDIT}/*.vcf.gz
    // Without *SmallIndels.vcf.gz from Manta, and *.genome.vcf.gz from Strelka
    // The small snippet `vcf.minus(vcf.fileName)[-2]` catches idSample
    // This field is used to output final annotated VCFs in the correct directory
      Channel.empty().mix(
        Channel.fromPath("${params.outdir}/VariantCalling/*/HaplotypeCaller/*.vcf.gz")
          .flatten().map{vcf -> ['HaplotypeCaller', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Manta/*[!candidate]SV.vcf.gz")
          .flatten().map{vcf -> ['Manta', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Mutect2/*.vcf.gz")
          .flatten().map{vcf -> ['Mutect2', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAseq/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAseq', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonDNAscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonDNAscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/SentieonTNscope/*.vcf.gz")
          .flatten().map{vcf -> ['SentieonTNscope', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/Strelka/*{somatic,variant}*.vcf.gz")
          .flatten().map{vcf -> ['Strelka', vcf.minus(vcf.fileName)[-2].toString(), vcf]},
        Channel.fromPath("${params.outdir}/VariantCalling/*/TIDDIT/*.vcf.gz")
          .flatten().map{vcf -> ['TIDDIT', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
      ).choice(vcfToAnnotate, vcfNoAnnotate) {
        annotate_tools == [] || (annotate_tools != [] && it[0] in annotate_tools) ? 0 : 1
      }
    } else if (annotate_tools == []) {
    // Annotate user-submitted VCFs
    // If user-submitted, Sarek assume that the idSample should be assumed automatically
      vcfToAnnotate = Channel.fromPath(tsvPath)
        .map{vcf -> ['userspecified', vcf.minus(vcf.fileName)[-2].toString(), vcf]}
    } else exit 1, "specify only tools or files to annotate, not both"

    vcfNoAnnotate.close()
    vcfAnnotation = vcfAnnotation.mix(vcfToAnnotate)
}

// as now have the list of VCFs to annotate, the first step is to annotate with allele frequencies, if there are any

(vcfSnpeff, vcfVep) = vcfAnnotation.into(2)

vcfVep = vcfVep.map {
  variantCaller, idSample, vcf ->
  [variantCaller, idSample, vcf, null]
}

// STEP SNPEFF

process Snpeff {
    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_snpEff.ann.vcf") null
        else "Reports/${idSample}/snpEff/${it}"
    }

    input:
        set variantCaller, idSample, file(vcf) from vcfSnpeff
        file(dataDir) from ch_snpeff_cache
        val snpeffDb from ch_snpeff_db

    output:
        set file("${reducedVCF}_snpEff.genes.txt"), file("${reducedVCF}_snpEff.html"), file("${reducedVCF}_snpEff.csv") into snpeffReport
        set variantCaller, idSample, file("${reducedVCF}_snpEff.ann.vcf") into snpeffVCF

    when: 'snpeff' in tools || 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    cache = (params.snpeff_cache && params.annotation_cache) ? "-dataDir \${PWD}/${dataDir}" : ""
    """
    snpEff -Xmx${task.memory.toGiga()}g \
        ${snpeffDb} \
        -csvStats ${reducedVCF}_snpEff.csv \
        -nodownload \
        ${cache} \
        -canon \
        -v \
        ${vcf} \
        > ${reducedVCF}_snpEff.ann.vcf

    mv snpEff_summary.html ${reducedVCF}_snpEff.html
    """
}

snpeffReport = snpeffReport.dump(tag:'snpEff report')

// STEP COMPRESS AND INDEX VCF.1 - SNPEFF

process CompressVCFsnpEff {
    tag "${idSample} - ${vcf}"

    publishDir "${params.outdir}/Annotation/${idSample}/snpEff", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from snpeffVCF

    output:
        set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into (compressVCFsnpEffOut)

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}

compressVCFsnpEffOut = compressVCFsnpEffOut.dump(tag:'VCF')

// STEP VEP.1

process VEP {
    label 'VEP'
    label 'cpus_4'

    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        set variantCaller, idSample, file(vcf), file(idx) from vcfVep
        file(dataDir) from ch_vep_cache
        val cache_version from ch_vep_cache_version
        file(cadd_InDels) from ch_cadd_indels
        file(cadd_InDels_tbi) from ch_cadd_indels_tbi
        file(cadd_WG_SNVs) from ch_cadd_wg_snvs
        file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

    output:
        set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCF
        file("${reducedVCF}_VEP.summary.html") into vepReport

    when: 'vep' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome

    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
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

vepReport = vepReport.dump(tag:'VEP')

// STEP VEP.2 - VEP AFTER SNPEFF

process VEPmerge {
    label 'VEP'
    label 'cpus_4'

    tag "${idSample} - ${variantCaller} - ${vcf}"

    publishDir params.outdir, mode: params.publish_dir_mode, saveAs: {
        if (it == "${reducedVCF}_VEP.summary.html") "Reports/${idSample}/VEP/${it}"
        else null
    }

    input:
        set variantCaller, idSample, file(vcf), file(idx) from compressVCFsnpEffOut
        file(dataDir) from ch_vep_cache
        val cache_version from ch_vep_cache_version
        file(cadd_InDels) from ch_cadd_indels
        file(cadd_InDels_tbi) from ch_cadd_indels_tbi
        file(cadd_WG_SNVs) from ch_cadd_wg_snvs
        file(cadd_WG_SNVs_tbi) from ch_cadd_wg_snvs_tbi

    output:
        set variantCaller, idSample, file("${reducedVCF}_VEP.ann.vcf") into vepVCFmerge
        file("${reducedVCF}_VEP.summary.html") into vepReportMerge

    when: 'merge' in tools

    script:
    reducedVCF = reduceVCF(vcf.fileName)
    genome = params.genome == 'smallGRCh37' ? 'GRCh37' : params.genome
    dir_cache = (params.vep_cache && params.annotation_cache) ? " \${PWD}/${dataDir}" : "/.vep"
    cadd = (params.cadd_cache && params.cadd_wg_snvs && params.cadd_indels) ? "--plugin CADD,whole_genome_SNVs.tsv.gz,InDels.tsv.gz" : ""
    genesplicer = params.genesplicer ? "--plugin GeneSplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/bin/genesplicer,/opt/conda/envs/nf-core-sarek-${workflow.manifest.version}/share/genesplicer-1.0-1/human,context=200,tmpdir=\$PWD/${reducedVCF}" : "--offline"
    """
    mkdir ${reducedVCF}

    vep \
        -i ${vcf} \
        -o ${reducedVCF}_VEP.ann.vcf \
        --assembly ${genome} \
        --species ${params.species} \
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

vepReportMerge = vepReportMerge.dump(tag:'VEP')

vcfCompressVCFvep = vepVCF.mix(vepVCFmerge)

// STEP COMPRESS AND INDEX VCF.2 - VEP

process CompressVCFvep {
    tag "${idSample} - ${vcf}"

    publishDir "${params.outdir}/Annotation/${idSample}/VEP", mode: params.publish_dir_mode

    input:
        set variantCaller, idSample, file(vcf) from vcfCompressVCFvep

    output:
        set variantCaller, idSample, file("*.vcf.gz"), file("*.vcf.gz.tbi") into compressVCFOutVEP

    script:
    """
    bgzip < ${vcf} > ${vcf}.gz
    tabix ${vcf}.gz
    """
}