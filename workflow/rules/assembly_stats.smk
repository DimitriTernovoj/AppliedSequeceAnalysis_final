rule quast:
	input:
		"results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta" if config["assembly_type"] == "de_novo" else "results/reference_based_assembly/ivar/{Sample}/assembly.fa"  
	output:
		"results/assembly_statistics/{Sample}/report.txt"
	log:
		"logs/assembly_statistics/{Sample}.log"
	threads: 8
	conda:
		"../envs/assembly.yaml"
	params:
		ref = "-r " +  config["reference_genome"] if config["reference_genome"] != "" else ""
	shell:
		"quast.py -l {wildcards.Sample} -o results/assembly_statistics/{wildcards.Sample}/ {params.ref} -t {threads} {input} 1> {log} 2>&1" 

rule multiqc_assembly:
    input:
        expand("results/assembly_statistics/{Sample}/report.txt", Sample=list(samples.index))
    output:
        "results/multiqc/multiqc_assembly.html"
    log:
        "logs/multiqc/multiqc_assembly.log"
    conda:
        "../envs/qc.yaml"
    params:
        outdir = "results/multiqc",
        input = "results/assembly_statistics/",
    shell:
        "multiqc -n multiqc_assembly.html -o {params.outdir} {params.input} 1> {log} 2>&1"

