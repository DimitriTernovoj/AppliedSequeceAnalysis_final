rule prokka:
	input:
		"results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta" if config["assembly_type"] == "de_novo" else "results/reference_based_assembly/ivar/{Sample}/assembly.fa"
	output:
		"results/annotation/{Sample}/{Sample}.gff"
	log:
		"logs/annotation/{Sample}.log"
	threads: workflow.cores
        conda:
                "../envs/annotation.yaml"
	params:
		genus = " --genus " + config["prokka"]["genus"] + " -usegenus " if config["prokka"]["genus"] != "" else "",
		species = " --species " + config["prokka"]["species"] if config["prokka"]["species"] != "" else "",
		extra = config["prokka"]["extra"]
	shell:
		"prokka --force --outdir results/annotation/{wildcards.Sample} --prefix {wildcards.Sample} {params.genus} {params.species} {params.extra} --cpus {threads} {input} 1> {log} 2>&1"

rule outgroup_prokka:
	input:
		config["phylogeny"]["outgroup"]
	output:
		expand("results/annotation/outgroup/{outgroup_sample_name}.gff", outgroup_sample_name = config["phylogeny"]["outgroup_sample_name"])
	log:
		"logs/annotation/outgroup.log"
	threads: workflow.cores
	conda:
		"../envs/annotation.yaml"
	params:
                outgroup_sample_name = config["phylogeny"]["outgroup_sample_name"],
		genus = " --genus " + config["prokka"]["genus"] + " -usegenus " if config["prokka"]["genus"] != "" else "",
                species = " --species " + config["prokka"]["species"] if config["prokka"]["species"] != "" else "",
		extra = config["prokka"]["extra"]
	shell:
		"prokka --force --outdir results/annotation/outgroup/ --prefix {params.outgroup_sample_name} {params.genus} {params.species} {params.extra} --cpus {threads} {input} 1> {log} 2>&1"
