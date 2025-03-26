/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT MODULES / SUBWORKFLOWS / FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
include { BWA_INDEX              } from '../modules/nf-core/bwa/index/main'
include { BWA_MEM                } from '../modules/nf-core/bwa/mem/main'

include { SAMTOOLS_FIXMATE } from '../modules/nf-core/samtools/fixmate/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT1 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_SORT as SAMTOOLS_SORT2 } from '../modules/nf-core/samtools/sort/main'
include { SAMTOOLS_MARKDUP } from '../modules/nf-core/samtools/markdup/main'
include { SAMTOOLS_INDEX } from '../modules/nf-core/samtools/index/main'

include { BEDTOOLS_GENOMECOV } from '../modules/nf-core/bedtools/genomecov/main'
include { BEDTOOLS_MERGE } from '../modules/nf-core/bedtools/merge/main'

include { GAWK } from '../modules/nf-core/gawk/main'
// include { GUNZIP as GZIP } from '../modules/nf-core/gunzip/main'
include { GUNZIP as GZIP } from '../modules/local/gunzip/main'

//include { BCFTOOLS_MPILEUP       } from '../modules/nf-core/bcftools/mpileup/main'
include { BCFTOOLS_MPILEUP       } from '../modules/local/bcftools/mpileup/main'
// include { BCFTOOLS_CALL          } from '../modules/nf-core/bcftools/call/main'
// include { BCFTOOLS_VIEW } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_SITES       } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_FILTER as BCFTOOLS_FILTER_SAMPLES     } from '../modules/nf-core/bcftools/filter/main'
include { BCFTOOLS_VIEW          } from '../modules/nf-core/bcftools/view/main'
include { BCFTOOLS_STATS         } from '../modules/local/bcftools/stats/main'
include { BCFTOOLS_CONSENSUS     } from '../modules/nf-core/bcftools/consensus/main'

include { paramsSummaryMap       } from 'plugin/nf-schema'
// include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_multilocus_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/


workflow MULTILOCUS {

    take:
    ch_samplesheet  // channel: samplesheet read in from --input

    main:
    // ch_versions = Channel.empty()

    // Create input channels
    //ch_input_fasta = Channel.empty()
    //ch_input_reads = Channel.empty()
    //ch_samplesheet.map { it -> println it}
    ch_ref = Channel.value([params.refname, params.ref])

    //
    // MODULE: Run bwa index
    //
    BWA_INDEX (
       ch_ref
    )
    // BWA_INDEX.out.index.view()

    //ch_id = ch_samplesheet.map { it -> it[0].id }
    //ch_sex = ch_samplesheet.map { it -> it[1] }
    //ch_pop = ch_samplesheet.map { it -> it[2] }
    ch_reads = ch_samplesheet.map { it -> [it[0], it[3]] }

    //ch_id.view()
    //ch_sex.view()
    //ch_pop.view()
    //ch_reads.view()

    //
    // MODULE: Read mapping
    //
    BWA_MEM (
        ch_reads,
        BWA_INDEX.out.index,
        ch_ref,
        false
    )
    // BWA_MEM.out.bam.view()

    // sort by read names
    SAMTOOLS_SORT1 (
        BWA_MEM.out.bam,
        ch_ref
    )

    // fill in mate coordinates and insert size
    SAMTOOLS_FIXMATE (
        SAMTOOLS_SORT1.out.bam
    )

    // sort by reference chromosomes and coordinates
    SAMTOOLS_SORT2 (
        SAMTOOLS_FIXMATE.out.bam,
        ch_ref
    )

    // remove all duplicates
    SAMTOOLS_MARKDUP (
        SAMTOOLS_SORT2.out.bam,
        ch_ref
    )

    SAMTOOLS_INDEX (
        SAMTOOLS_MARKDUP.out.bam
    )

    // create depth files
    ch_genomecov_input = SAMTOOLS_MARKDUP.out.bam.map { meta, bam -> [ meta, bam, 1 ] }
    // ch_genomecov_input.view()

    BEDTOOLS_GENOMECOV (
        ch_genomecov_input,
        [],
        'bed',
        false
    )

    // get intervals with depth >= minDP
    GAWK (
        BEDTOOLS_GENOMECOV.out.genomecov,
        [],
        false
    )

    // clean up intervals
    BEDTOOLS_MERGE (
        GAWK.out.output
    )

    GZIP (
        BEDTOOLS_MERGE.out.bed
    )


    // variant calling per population
    // SAMTOOLS_MARKDUP.out.bam.view()

    // ch_samplesheet.groupTuple(by: 2).view()

    // SAMTOOLS_MARKDUP.out.bam.view()

    // add info from ch_samplesheet to SAMTOOLS_MARKDUP.out.bam by sample id (by: 0)
    // then group by pop (by: 2)
    ch_pop = ch_samplesheet.map { it -> [it[0], it[2]] }
    // ch_pop.view()

    // SAMTOOLS_MARKDUP.out.bam.combine(ch_pop, by: 0).view()
    // SAMTOOLS_MARKDUP.out.bam.combine(ch_pop, by: 0).groupTuple(by: 2).view()
    // ch_bam_pop = SAMTOOLS_MARKDUP.out.bam.combine(ch_pop, by: 0).groupTuple(by: 2).map { it -> [ it[0], it[1], [] ] }

    // use pop as id instead of sample name
    ch_bam_pop = SAMTOOLS_MARKDUP.out.bam.combine(ch_pop, by: 0).groupTuple(by: 2).map { it -> [ [id: it[2][0]], it[1], [] ] }

    // ch_bam_pop.view()

    // ch_bam_pop.map{ it -> it[1] }.view()

    BCFTOOLS_MPILEUP (
        ch_bam_pop,
        ch_ref,
        false
    )

    // mark site filter
    ch_vcf_pop = BCFTOOLS_MPILEUP.out.vcf.combine(BCFTOOLS_MPILEUP.out.tbi, by: 0)
    // ch_vcf_pop.view()

    BCFTOOLS_FILTER_SITES (
        ch_vcf_pop
    )

    // mark sample filter
    ch_vcf_pop_mark_sites = BCFTOOLS_FILTER_SITES.out.vcf.combine(BCFTOOLS_FILTER_SITES.out.tbi, by: 0)

    BCFTOOLS_FILTER_SAMPLES (
        ch_vcf_pop_mark_sites
    )

    // remove sites that do not pass site-filter
    ch_vcf_pop_mark_samples = BCFTOOLS_FILTER_SAMPLES.out.vcf.combine(BCFTOOLS_FILTER_SAMPLES.out.tbi, by: 0)
    // ch_vcf_pop_mark_samples.view()

    BCFTOOLS_VIEW (
        ch_vcf_pop_mark_samples,
        [], [], []
    )
    // TODO: site filtering requires a site mask later in bcftools consensus


    ch_vcf_pop_filt = BCFTOOLS_VIEW.out.vcf.combine(BCFTOOLS_VIEW.out.tbi, by: 0)
    // ch_vcf_pop_filt.view()

    BCFTOOLS_STATS (
        ch_vcf_pop_filt
    )

    






    // FASTQC (
    //     ch_samplesheet
    // )
    // ch_multiqc_files = ch_multiqc_files.mix(FASTQC.out.zip.collect{it[1]})
    // ch_versions = ch_versions.mix(FASTQC.out.versions.first())

    // //
    // // Collate and save software versions
    // //
    // softwareVersionsToYAML(ch_versions)
    //     .collectFile(
    //         storeDir: "${params.outdir}/pipeline_info",
    //         name: 'nf_core_'  +  'multilocus_software_'  + 'mqc_'  + 'versions.yml',
    //         sort: true,
    //         newLine: true
    //     ).set { ch_collated_versions }


    // //
    // // MODULE: MultiQC
    // //
    // ch_multiqc_config        = Channel.fromPath(
    //     "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    // ch_multiqc_custom_config = params.multiqc_config ?
    //     Channel.fromPath(params.multiqc_config, checkIfExists: true) :
    //     Channel.empty()
    // ch_multiqc_logo          = params.multiqc_logo ?
    //     Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
    //     Channel.empty()

    // summary_params      = paramsSummaryMap(
    //     workflow, parameters_schema: "nextflow_schema.json")
    // ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    // ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
    //     file(params.multiqc_methods_description, checkIfExists: true) :
    //     file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    // ch_methods_description                = Channel.value(
    //     methodsDescriptionText(ch_multiqc_custom_methods_description))

    // ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    // ch_multiqc_files = ch_multiqc_files.mix(
    //     ch_methods_description.collectFile(
    //         name: 'methods_description_mqc.yaml',
    //         sort: true
    //     )
    // )

    // MULTIQC (
    //     ch_multiqc_files.collect(),
    //     ch_multiqc_config.toList(),
    //     ch_multiqc_custom_config.toList(),
    //     ch_multiqc_logo.toList(),
    //     [],
    //     []
    // )

    // emit:multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    // versions       = ch_versions                 // channel: [ path(versions.yml) ]

}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
