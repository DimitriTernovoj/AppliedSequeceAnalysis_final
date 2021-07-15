def fastqc_input():
	fastq_list = []
	for Sample in samples.index:
		for i in ["fq1","fq2"]:
			fastq_list.append(samples.at[Sample,i])
	
	for i in fastq_list:
		fastq = lambda wildcards: "/".join((i.split(".")[0]).split("/")[0:-1]) + "/" +  wildcards.Sample + wildcards.Orientation + "." + ".".join(i.split(".")[1:]) if wildcards.Sample in samples.index and wildcards.Orientation in ["_1","_2"] else ""
		return fastq


rule first_fastqc:
    input:
        fastqc_input()
    output:
        html = "results/quality_control/first_fastqc/{Sample}/{Sample}{Orientation}.html",
        zip = "results/quality_control/first_fastqc/{Sample}/{Sample}{Orientation}_fastqc.zip",
    log:
        "logs/fastqc/{Sample}{Orientation}.log"
    threads:
        1
    wrapper:
        "0.77.0/bio/fastqc"


rule trimmomatic:
    input:
        fastq1=lambda wildcards: samples.at[wildcards.Sample,'fq1'] if wildcards.Sample in samples.index else '',
        fastq2=lambda wildcards: samples.at[wildcards.Sample,'fq2'] if wildcards.Sample in samples.index else '',
    output:
        "results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz",
        "results/trimmed_reads/{Sample}_1_trimmed_pe_unpaired.fastq.gz",
        "results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz",
        "results/trimmed_reads/{Sample}_2_trimmed_pe_unpaired.fastq.gz",
    log:
        "logs/trimmomatic/{Sample}.log"
    conda:
        "../envs/qc.yaml"
    params:
        adapter_trim = "ILLUMINACLIP:" + config["trimmomatic"]["adapters"] + ":2:25:8" if config["trimmomatic"]["adapters"] != "" else "",
	sliding_window = "SLIDINGWINDOW:" + config["trimmomatic"]["sliding_window"]["window_size"] + ":" + config["trimmomatic"]["sliding_window"]["requiredQuality"] if config["trimmomatic"]["sliding_window"]["window_size"] != "" and config["trimmomatic"]["sliding_window"]["requiredQuality"] != "" else "SLIDINGWINDOW:4:15",
	minlen = "MINLEN:" + config["trimmomatic"]["minlen"] if config["trimmomatic"]["minlen"] != "" else "MINLEN:40"
    shell:
        "trimmomatic PE -trimlog {log} {input.fastq1} {input.fastq2} {output} {params.adapter_trim} {params.sliding_window} {params.minlen} 1> {log} 2>&1"


rule second_fastqc:
    input:
        fastq1="results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz",
        fastq2="results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz",
    output:
        "results/quality_control/second_fastqc/{Sample}/{Sample}_1_trimmed_pe_fastqc.zip",
        "results/quality_control/second_fastqc/{Sample}/{Sample}_2_trimmed_pe_fastqc.zip"
    log:
        "logs/fastqc_trimmed/{Sample}.log"
    threads:
        4
    conda:
        "../envs/qc.yaml"
    params:
        outdir = "results/quality_control/second_fastqc/{Sample}"
    shell:
        "fastqc -t {threads} {input.fastq1} {input.fastq2} -o {params.outdir} 1> {log} 2>&1"

#rule second_fastqc:
#    input:
#        fastq = "results/trimmed_reads/{Sample}{Orientation}_trimmed.fastq.gz"
#    output:
#        html = "results/quality_control/second_fastqc/{Sample}/{Sample}{Orientation}_trimmed.html",
#        zip = "results/quality_control/second_fastqc/{Sample}/{Sample}{Orientation}_trimmed_fastqc.zip",
#    log:
#        "logs/fastqc/{Sample}{Orientation}_trimmed.log"
#    threads:
#        1
#    wrapper:
#        "0.77.0/bio/fastqc"


rule multiqc:
    input:
        expand(["results/quality_control/second_fastqc/{Sample}/{Sample}{Orientation}_trimmed_pe_fastqc.zip", "results/quality_control/first_fastqc/{Sample}/{Sample}{Orientation}.html"], Sample=list(samples.index), Orientation=["_1","_2"])
    output:
        "results/multiqc/multiqc_qc.html"
    log:
        "logs/multiqc/multiqc.log"
    conda:
        "../envs/qc.yaml"
    params:
        outdir = "results/multiqc",
        input = "results/quality_control/",
    shell:
        "multiqc -n multiqc_qc.html -o {params.outdir} {params.input} 1> {log} 2>&1"


#rule multiqc:
#    input:
#        expand(["results/quality_control/second_fastqc/{Sample}/{Sample}{Orientation}_trimmed.html", "results/quality_control/first_fastqc/{Sample}/{Sample}{Orientation}.html"], Sample=list(samples.index$
#    output:
#        "results/multiqc/multiqc_qc.html"
#    log:
#        "logs/multiqc/multiqc_qc.log"
#    wrapper:
#        "0.77.0/bio/multiqc"

#rule multiqc_qc:
#    input:
#        expand(["results/quality_control/second_fastqc/{Sample}_trimmed/{Sample}{Orientation}_trimmed.html"], Sample=list(samples.index), Orientation = ["_1","_2"])
#    output:
#        "results/multiqc/multiqc_qc.html"
#    log:
#        "logs/multiqc/mutltiqc_qc.log"
#    wrapper:
#        "0.77.0/bio/multiqc"

