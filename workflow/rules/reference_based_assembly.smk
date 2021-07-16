rule bwa_index:
	input:
		expand(["results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz", "results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz"], Sample=list(samples.index))
	output:
		"results/reference_based_assembly/ref_index/ref_index.bwt"
	conda:
		"../envs/reference_based_assembly.yaml"
	log:
		"logs/reference_based_assembly/bwa_index.log"
	params:
		ref = config["reference_genome"]
	shell:
		"""
		cd results/reference_based_assembly/ref_index
		bwa index -p ref_index {params.ref} 1> ../../../{log} 2>&1
		cd ../../../
		"""
		
rule bwa_mem:
	input:
		fastq1="results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz",
		fastq2="results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz",
		index="results/reference_based_assembly/ref_index/ref_index.bwt"
	output:
		"results/reference_based_assembly/{Sample}/bwa_assembly.sam"
	conda:
		"../envs/reference_based_assembly.yaml"
	log:
		"logs/reference_based_assembly/bwa_mem/{Sample}.log"
	params:
		ref = config["reference_genome"],
		ref_index = "results/reference_based_assembly/ref_index/ref_index"
	threads: workflow.cores
	shell:
		"bwa mem -t {threads} {params.ref_index} {input.fastq1} {input.fastq2} -o {output} 1> {log} 2>&1"
		
rule samtools_samtobam:
	input:
		"results/reference_based_assembly/{Sample}/bwa_assembly.sam"
	output:
		"results/reference_based_assembly/{Sample}/bwa_assembly.bam"
	conda:
		"../envs/reference_based_assembly.yaml"
	log:
		"logs/reference_based_assembly/samtools/{Sample}.log"
	shell:
		"samtools view -b {input} -o {output} > {log}"

rule samtools_sorting:
	input:
		"results/reference_based_assembly/{Sample}/bwa_assembly.bam"
	output:
		"results/reference_based_assembly/{Sample}/bwa_assembly_sorted.bam"
	conda:
		"../envs/reference_based_assembly.yaml"
	log:
		"logs/reference_based_assembly/samtools/{Sample}_sorted.log"
	threads: workflow.cores * 0.5
	shell:
		"samtools sort -@ {threads} -o {output} -O bam {input} 1> {log} 2>&1"

rule samtools_indexing:
        input:
                "results/reference_based_assembly/{Sample}/bwa_assembly_sorted.bam"
        output:
                "results/reference_based_assembly/{Sample}/bwa_assembly_sorted.bai"
        conda:
                "../envs/reference_based_assembly.yaml"
        log:
                "logs/reference_based_assembly/samtools/{Sample}_indexed.log"
        threads: workflow.cores * 0.5
        shell:
                "samtools index -@ {threads} {input} {output} 1> {log} 2>&1"

rule ivar:
        input:  
                "results/reference_based_assembly/{Sample}/bwa_assembly_sorted.bam",
                "results/reference_based_assembly/{Sample}/bwa_assembly_sorted.bai"
        output: 
                "results/reference_based_assembly/ivar/{Sample}/assembly.fa"
        conda: 
               "../envs/reference_based_assembly.yaml"
        log:
                "logs/reference_based_assembly/ivar/{Sample}_ivar.log"
        shell:  
                "(samtools mpileup -A -d 0 -Q 0 {input[0]} | ivar consensus -p {output} -i {wildcards.Sample}) 1> {log} 2>&1"

