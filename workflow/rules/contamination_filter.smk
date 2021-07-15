rule kraken2:
	input:
		fq1 ="results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz",
		fq2 ="results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz"
	output:
		"results/kraken2/report_{Sample}.tsv",
		"results/kraken2/kraken_rep_{Sample}.txt"
	log:
		"logs/contamination_filter/kraken2/{Sample}.log"
	conda:
		"../envs/contamination_filter.yaml"
	params:
		database = config["contamination_filter"]["kraken2_database"]
	threads: 8
	shell:
		"kraken2 --db {params.database} --threads {threads} --gzip-compressed --paired {input.fq1} {input.fq2} --report {output[0]} --output {output[1]} 1> {log} 2>&1"

rule contamination_info:
	input:
		expand("results/kraken2/report_{Sample}.tsv", Sample=list(samples.index))
	output:
		"results/multiqc/multiqc_contamination.html"
	log:
		"logs/multiqc/multiqc_contamination.log"	
	conda:
		"../envs/qc.yaml"
	params:
		outdir = "results/multiqc",
		input = "results/kraken2"
	shell:
		"multiqc -n multiqc_contamination.html -o {params.outdir} {params.input} 1> {log} 2>&1"

rule kraken_tools:
	input:
		"results/kraken2/kraken_rep_{Sample}.txt",
		"results/kraken2/report_{Sample}.tsv",
                "results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz",
                "results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz"
	output:
		"results/decontaminated_reads/{Sample}_1.fq", #"results/decontaminated_reads/{Sample}_1.fastq.gz",
		"results/decontaminated_reads/{Sample}_2.fq"
	log:
		"logs/contamination_filter/krakentools/{Sample}.log"
	conda:
		"../envs/contamination_filter.yaml"
	params:
		taxonomy_list = config["contamination_filter"]["krakentools"]["taxonomy_IDs"],
		filter_type = "--exclude" if config["contamination_filter"]["krakentools"]["filter_type"] == "negative" else "",
		parents = "--include-parents" if config["contamination_filter"]["krakentools"]["include_parents"] == "yes" else "",
		children = "--include-children" if config["contamination_filter"]["krakentools"]["include_children"] == "yes" else ""
	shell:
		"extract_kraken_reads.py -k {input[0]} -r {input[1]} -s1 {input[2]} -s2 {input[3]} -o {output[0]} -o2 {output[1]} -t {params.taxonomy_list} {params.filter_type} {params.parents} {params.children} --fastq-output 1> {log} 2>&1"

#rule kraken2_temp:
#        input:
#                "results/decontaminated_reads/{Sample}_1.fq",
#                "results/decontaminated_reads/{Sample}_2.fq"
#        output:
#                "results/temp/report_{Sample}.tsv"
#        log:
#                "logs/contamination_filter/kraken2/{Sample}.log"
#        conda:
#                "../envs/contamination_filter.yaml"
#        params:
#                database = config["contamination_filter"]["kraken2_database"]
#        threads: 8
#        shell:
#                "kraken2 --db {params.database} --threads {threads} --paired {input[0]} {input[1]} --report {output}"

#rule contamination_info_temp:
#        input:   
#                expand("results/temp/report_{Sample}.tsv", Sample=list(samples.index))
#        output:
#                "results/temp/multiqc/multiqc_contamination.html"
#        #log:
#        #        "logs/multiqc/multiqc_contamination.log"
#        conda:
#                "../envs/qc.yaml"
#        params:
#                outdir = "results/temp/multiqc",
#                input = "results/temp"
#        shell:
#                "multiqc -n multiqc_contamination.html -o {params.outdir} {params.input}"

