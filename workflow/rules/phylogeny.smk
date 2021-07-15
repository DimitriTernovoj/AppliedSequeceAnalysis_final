rule raxml_ng:
	input:
		"results/core_genomes/roary/core_gene_alignment.aln"
	output:
		"results/phylogeny/phylogeny.raxml.bestTree"
	log:
		"logs/phylogeny/raxml-ng.log"
	params:
		outgroup = "--outgroup " + config["phylogeny"]["outgroup_sample_name"] if config["phylogeny"]["outgroup_sample_name"] != "" else "",
		extra = config["phylogeny"]["extra"],
		bs_trees = config["phylogeny"]["bs_trees"]
	conda:
		"../envs/phylogeny.yaml"
	shell:
		"raxml-ng --msa {input} {params.outgroup} {params.extra} --msa-format FASTA --all --bs-trees {params.bs_trees} --model GTR+G --seed 15 --prefix results/phylogeny/phylogeny 1> {log} 2>&1"

rule visualization:
	input:
		"results/phylogeny/phylogeny.raxml.bestTree"
	output:
		"results/phylogeny/newick_tree.png"
	conda:
		"../envs/phylogeny.yaml"
	script:
		"../scripts/newick_tree.py"
