rule RGI:
	input:
		"results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta" if config["assembly_type"] == "de_novo" else "results/reference_based_assembly/ivar/{Sample}/assembly.fa"
	output:
		"results/antibiotic_resistance_genes/{Sample}.txt"
	log:
		"logs/antibiotic_resistance_genes/{Sample}.log"
	threads: workflow.cores
	conda:
		"../envs/antibiotic_resistance_genes.yaml"
	params:
		output = "results/antibiotic_resistance_genes/{Sample}",
		extra = config["antibiotic_resistance_genes"]["extra"],
		card_database = "rgi load --card_json " + config["antibiotic_resistance_genes"]["card_database"] if config["antibiotic_resistance_genes"]["card_database"] != "" else ""
	shell:
		"""
		{params.card_database}
		rgi main -i {input} -o {params.output} {params.extra} --n {threads} --input_type contig --clean  1> {log} 2>&1
		"""

rule RGI_heatmap:
	input:
		expand("results/antibiotic_resistance_genes/{Sample}.txt", Sample=list(samples.index))
	output:
		"results/antibiotic_resistance_genes/heatmap.png"
	conda:
		"../envs/antibiotic_resistance_genes.yaml"
	log:
		"logs/antibiotic_resistance_genes/heatmap.log"
	params:
		input_path = "results/antibiotic_resistance_genes/",
		output = "results/antibiotic_resistance_genes/heatmap"
	shell:
		"""
		rgi heatmap -i {params.input_path} -o {params.output} 1> {log} 2>&1
		mv results/antibiotic_resistance_genes/heatmap*.png results/antibiotic_resistance_genes/heatmap.png
		"""
