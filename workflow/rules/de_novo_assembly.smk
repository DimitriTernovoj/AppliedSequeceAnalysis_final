rule bbnorm:
    input:
        "results/trimmed_reads/{Sample}_1_trimmed_pe.fastq.gz" if config["contamination_filter"]["execution"] == "no" else "results/decontaminated_reads/{Sample}_1.fq",
        "results/trimmed_reads/{Sample}_2_trimmed_pe.fastq.gz" if config["contamination_filter"]["execution"] == "no" else "results/decontaminated_reads/{Sample}_2.fq"
    output:
        "results/de_novo_assembly/bbnorm/normalized_{Sample}_1.fastq.gz",
        "results/de_novo_assembly/bbnorm/normalized_{Sample}_2.fastq.gz"
    log:
        "logs/de_novo_assembly/bbnorm/normalized_{Sample}.log"
    conda:
        "../envs/assembly.yaml"
    threads: workflow.cores * 0.5
    shell:
        "bbnorm.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} target=100 min=5 t={threads} 1>{log} 2>&1"

rule spades:
    input:
        fastq1="results/de_novo_assembly/bbnorm/normalized_{Sample}_1.fastq.gz",
        fastq2="results/de_novo_assembly/bbnorm/normalized_{Sample}_2.fastq.gz",
    output:
        "results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta"
    log:
        "logs/de_novo_assembly/Spades/{Sample}.log"
    threads: workflow.cores
    conda:
        "../envs/assembly.yaml"
    shell:
        "spades.py -1 {input.fastq1} -2 {input.fastq2} -t {threads} --isolate -o results/de_novo_assembly/assembly/{wildcards.Sample} > {log}"

