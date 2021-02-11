/*
 * QC and normalization of methylation data
 */

date = java.time.LocalDate.now()

params.data = path("data/")
params.configMEDALL = path("configs/QC_Medall.R")
params.configEsteller = path("configs/QC_Estller.R")

params.genosPath = ""
params.methyAnnot = ""

phenoPath = file(params.phenoPath)
genosPath = file(params.genosPath)
methyAnnot = file(params.methyAnnot)

workflowInfo = """
Project : $workflow.projectDir
Nextflow version: $nextflow.version - build: $nextflow.build
Date: $date
Git info: $workflow.repository - $workflow.revision [$workflow.commitId]
Cmd line: $workflow.commandLine
"""

// Run QC and normalization pipeline
process runQC_normalization {

  publishDir "${params.qcdir}/${date}", mode: 'copy'

  if ( params.version != null ){
    publishDir "${params.qcdir}/${params.version}", mode: 'copy'
  }

  input:
  path inFold from Channel.fromPath("${params.inFold}")
  val sampleSheet from "${params.sampleSheet}"
  file(pheno) from phenoPath
  file(genos) from genosPath
  val logText from "$workflowInfo"
  val cores from "${params.cores}"

  output:
  file 'gset.Rdata' into gset
  file 'norm.beta.Rdata' into norm_betas
  file 'norm.obj.pc.Rdata' into norm_obj
  file 'norm.summary.Rdata' into norm_summary
  file 'pc.fit.pdf' into pc_fit
  file 'qc.objects.clean.Rdata' into qc_clean
  file 'qc.objects.Rdata' into qc_obj
  file 'qcsummary.clean.Rdata' into qc_sum_clean
  file 'qcsummary.Rdata' into qc_sum
  file 'removed.probes.txt' into badprobes
  file 'removed.samples.txt' into badsamples
  file 'log.txt' into logCh1

  """
  Rscript $baseDir/scripts/methylation_QC.R $inFold $sampleSheet $phenoPath $cores $genosPath
  echo "$logText" > log.txt
  """
}

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
}
