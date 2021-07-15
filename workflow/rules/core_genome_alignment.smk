#def roary_input(wildcards):
#	input = expand("results/annotation/{Sample}/{Sample}_annotated.gff", Sample=list(samples.index)) if 

def roary_input():
	sample_list = expand(["results/annotation/{Sample}/{Sample}.gff"], Sample=list(samples.index)) if config["phylogeny"]["outgroup"] == "" else expand(["results/annotation/{Sample}/{Sample}.gff","results/annotation/outgroup/{outgroup_sample_name}.gff"], Sample=list(samples.index), outgroup_sample_name = config["phylogeny"]["outgroup_sample_name"])
	if config["phylogeny"]["blacklist"] != "":
		blacklist = config["phylogeny"]["blacklist"].split(",")
		for i in blacklist:
			sample_list.remove(f"results/annotation/{i}/{i}.gff".format())

	return sample_list

rule roary:
	input:
		roary_input()
	output:
		"results/core_genomes/roary/core_gene_alignment.aln"
	log:
		"logs/core_genomes/roary.log"
	threads: 8
	conda:
		"../envs/core_genome_alignment.yaml"
	params:
		alignment_type = "--mafft" if config["roary"]["alignment_type"] == "fast" else "",
		extra = config["roary"]["extra"]
	shell:
		"""
		roary -p {threads} -f results/core_genomes/core_genome_inference/work_dir -e {params.alignment_type} {params.extra} {input} 1> {log} 2>&1
		mv results/core_genomes/core_genome_inference/work_dir/core_gene_alignment.aln {output}
		"""


#rule roary:
#        input:
#                expand("results/annotation/{Sample}/{Sample}_annotated.gff", Sample=list(samples.index)) if config["phylogeny"]["outgroup"] == "" else expand(["results/annotation/{Sample}/{Sample}_annota$
#        output:
#                "results/core_genomes/roary/core_gene_alignment.aln"
#        log:
#                "logs/core_genomes/roary.log"
#        threads: 8
#        conda:
#                "../envs/core_genome_alignment.yaml"
#        shell:
#                """
#                roary -p {threads} -f results/core_genomes/core_genome_inference/work_dir -e --mafft {input} 1> {log} 2>&1
#                mv results/core_genomes/core_genome_inference/work_dir/core_gene_alignment.aln {output}
#                """


#enc2xs -C
