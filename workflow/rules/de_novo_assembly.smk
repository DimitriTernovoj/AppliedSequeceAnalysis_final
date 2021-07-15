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
    threads: 8
    shell:
        "bbnorm.sh in1={input[0]} in2={input[1]} out1={output[0]} out2={output[1]} target=100 min=5 t={threads} 1>{log} 2>&1"

rule spades:
    input:
        fastq1="results/de_novo_assembly/bbnorm/normalized_{Sample}_1.fastq.gz",
        fastq2="results/de_novo_assembly/bbnorm/normalized_{Sample}_2.fastq.gz",
        #fastqc="results/second_fastqc/{Sample}/{Sample}_1_trimmed_pe_fastqc.zip",
    output:
        "results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta"
    log:
        "logs/de_novo_assembly/Spades/{Sample}.log"
    conda:
        "../envs/assembly.yaml"
    shell:
        "spades.py -1 {input.fastq1} -2 {input.fastq2} --isolate -o results/de_novo_assembly/assembly/{wildcards.Sample} > {log}"

#rule abacas:
#    input:
#        "results/de_novo_assembly/assembly/{Sample}/scaffolds.fasta",
#    output:
#        "results/de_novo_assembly/abacas/{Sample}/final_assembly.fasta"
#    log:
#        "logs/de_novo_assembly/abacas/{Sample}.log"
#    conda:
#        "../envs/assembly.yaml"
#    params:
#         ref = config["reference_genome"]
#    shell:
#        """
#        cd results/de_novo_assembly/abacas/{wildcards.Sample}/
#        abacas.pl -r {params.ref} -q ../../../../{input} -p nucmer -m -o final_assembly 1> "../../../../logs/de_novo_assembly/abacas/{wildcards.Sample}.log" 2>&1
#        cd ../../../../
#        """
