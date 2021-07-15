rule mlst:
	input:
		expand("results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta", Sample=list(samples.index)) if config["assembly_type"] == "de_novo" else expand("results/reference_based_assembly/ivar/{Sample}/assembly.fa", Sample=list(samples.index))
	output:
		"results/mlst/mlst.tsv"
	log:
		"logs/mlst/mlst.log"
	conda:
		"../envs/mlst.yaml"
	params:
		scheme = "--scheme " + config["mlst"]["scheme"] if config["mlst"]["scheme"] != "" else "",
		extra = config["mlst"]["extra"]
	shell:
		"mlst {params.scheme} {params.extra} {input} > {output} 2> {log}"

