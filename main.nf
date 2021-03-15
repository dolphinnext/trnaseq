$HOSTNAME = ""
params.outdir = 'results'  


if (!params.reads){params.reads = ""} 
if (!params.mate){params.mate = ""} 

Channel
	.fromFilePairs( params.reads , size: params.mate == "single" ? 1 : params.mate == "pair" ? 2 : params.mate == "triple" ? 3 : params.mate == "quadruple" ? 4 : -1 )
	.ifEmpty { error "Cannot find any reads matching: ${params.reads}" }
	.set{g_1_reads_g_0}

Channel.value(params.mate).set{g_2_mate_g_0}


process mimseq {

input:
 set val(name), file(reads) from g_1_reads_g_0
 val mate from g_2_mate_g_0


script:
"""
mimseq --species Hsap --cluster --cluster-id 0.95 --snp-tolerance --cca-analysis --threads 1 --min-cov 2000 --max-mismatches 0.1 --control-condition HEK293T -n hg19_test --out-dir hg19_HEK239vsK562 --max-multi 1 --remap --remap-mismatches 0.075 sampleData_HEKvsK562.txt
"""
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
