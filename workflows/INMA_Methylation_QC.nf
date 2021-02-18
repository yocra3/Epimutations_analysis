/*
 * QC and normalization of methylation data
 */

date = java.time.LocalDate.now()

if (params.inputTable) {
      Channel
            .fromPath(params.inputTable)
            .splitCsv(sep: "\t")
            .map { row -> [ row[0], file(row[1], checkIfExists: true), path(row[2], checkIfExists: true) ,
             file(row[3], checkIfExists: true) }
            .ifEmpty { exit 1, "params.inputTable was empty - no input files supplied" }
            .into { ch_input_QC }
} else { exit 1, "--inputTable should be text file with the cohort id, the path to IDAT folder, the path to a Rdata with the phenotypes." }

ch_geno = Channel.fromPath( [file("${params.genoPath}.bed", checkIfExists: true),
                            file("${params.genoPath}.bim", checkIfExists: true),
                            file("${params.genoPath}.fam", checkIfExists: true) })

workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

process getMethylationSNPs {

  output:
  file("snp-names.txt") into ch_methy_SNPs

  script:
  """
  Rscript 'writeLines(meffil::meffil.snp.names(), con="snp-names.txt")'
  """
}

process convertGenotypes {

  input:
  set file("geno.bed"), file("geno.bim"), file("geno.fam") from ch_geno
  file(selSNPs) from ch_methy_SNPs.collect()

  output:
  file("${genoID}.raw") into ch_geno_raw

  script:
  """
  plink --bfile geno --extract $selSNPs --recode A --out ${genoID}
  """
}

// Run QC and normalization pipeline
process runQC_normalization {

  publishDir "${params.outdir}/Methylation_QC/${date}/${cohort}", mode: 'copy', pattern: '*.Rdata'
  publishDir "${params.outdir}/Methylation_QC/${date}/${cohort}", mode: 'copy', pattern: '*.txt'
  publishDir "${params.outdir}/Methylation_QC/${date}/${cohort}", mode: 'copy', pattern: '*.pdf'
  publishDir "${params.reports}/Methylation_QC/${date}/${cohort}", mode: 'copy', pattern: '*.html'

  input:
  set val(cohort), file(config), path("idats/"), file("phenotypes.Rdata") from ch_input_QC
  file("genotypes.raw") from ch_geno_raw.collect()

  output:
  set val(cohort), file '*GenomicRatioSet.Rdata' into ch_raw_gset
  set val(cohort), file '*.norm.obj.pc.Rdata' into ch_raw_gset
  file("${cohort}*.Rdata")
  file("${cohort}*.txt")
  file("${cohort}*.html")
  file("${cohort}*.pdf")

  """
  run_methylation_QC.R $config $cohort
  """
}
/*
// Create GenomicRatioSets with filtered probes
process prepareGenomicRatioSets {

 	publishDir "${params.outdir}/${date}", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.outdir}/${params.version}", mode: 'copy'
  }

  input:
  file 'gset.Rdata' from gset
  file(methyAnnot) from methyAnnot
  val logText from "$workflowInfo"

  output:
  file 'gset.allProbes.Rdata' into gset_all
  file 'gset.allCpGs.Rdata' into gset_cpgs
  file 'gset.filterAnnotatedProbes.Rdata' into gset_filt
  file 'gset.autosomic.Rdata' into gset_aut

  file 'log.txt' into logCh2

  """
  Rscript $baseDir/scripts/createGenomicRatioSets.R gset.Rdata $methyAnnot
  echo "$logText" > log.txt
  """
}*/
