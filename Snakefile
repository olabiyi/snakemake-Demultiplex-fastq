from os import path, makedirs, rename, remove

# Show the rulegraph
#  snakemake -s Snakefile --rulegraph | dot -Tpng > rulegraph.png

configfile: "config.yaml"

onsuccess:
    print("Workflow completed without any error")

onerror:
    print("An error occurred")

# Should an excel file with barcodes and sample names be coverted to TSV ?
# Please note that it is assumed that the 5th and 2nd columns of sheet one
#  in the Excel file are the sample names and barcodes, respectively
isExcel=config['isExcel']  # True or False
EXCEL_FILE=config['EXCEL_FILE']

SAMPLES=config['SAMPLES']



rule all:
    input: 
        "02.Join_fastq/indices.fastq",
        "07.Count_Seqs/seqs_stat.txt"
      

rule Join_fastq:
    input:
        forward_index = "01.raw_data/index1.fastq",
        reverse_index = "01.raw_data/index2.fastq"
    output: "02.Join_fastq/indices.fastq"
    params:
        program=config['programs_path']['usearch']
    threads: 1
    shell:
        """
        {params.program} \
                 -fastq_join {input.forward_index} \
                 -reverse {input.reverse_index} \
                 -join_padgap "" \
                 -threads {threads} \
                 -fastqout {output}
        """

# Convert barcodes file from TSV to FASTA
rule Parse_barcodes:
    input: EXCEL_FILE if isExcel else "01.raw_data/sample2barcode.tsv"
    output: "03.Parse_barcodes/bar.fasta"
    params:
        isExcel=isExcel,
        program=config['programs_path']['xlsx2csv']
    threads: 1
    shell:
        """
        ISEXCEL={params.isExcel}
        if [ ${{ISEXCEL}} == True ]; then

             {params.program} -d "tab" -s 1 {input} | \
              awk 'BEGIN{{FS=OFS="\t"}} NR>1{{print($5,$2)}}' | \
              awk '{{print ">"$1"\\n"$2}}' > {output}

        else

            awk '{{print ">"$1"\\n"$2}}' {input} > {output}

        fi
        """
        

rule Reformat_barcodes:
    input: rules.Parse_barcodes.output
    output: "04.Reformat_barcodes/bar.fasta"
    threads: 1
    script:
        "./reformat_barcode.py"
        

# Demultiplex the reads per sample using the reformated barcodes
rule Demultiplex:
    input: 
        barcodes=rules.Reformat_barcodes.output,
        indices=rules.Join_fastq.output,
        forward="01.raw_data/read1.fastq",
        rev="01.raw_data/read2.fastq"
    output:
        forward="05.Demultiplex/demux_R1.fastq",
        rev="05.Demultiplex/demux_R2.fastq"
    params:
        program=config['programs_path']['usearch']
    threads: 1
    shell:
        """
        {params.program} \
                -fastx_demux {input.forward} \
                -reverse {input.rev} \
                -index {input.indices} \
                -barcodes {input.barcodes} \
                -fastqout {output.forward} \
                -output2 {output.rev}
        """



# Separate the demultiplexed fastq into folders containing 
# the forward and reverse reads for each sample

rule Split_forward:
    input: rules.Demultiplex.output.forward
    output: expand("06.Split/{sample}/{sample}_R1.fastq.gz", sample=SAMPLES)
    params:
        program=config['programs_path']['parallel'],
        samples=SAMPLES
    threads: 10
    shell:
        """
        function split_samples(){{
             
             local OUTDIR=$1
             local DEMUX_FILE=$2
             local SAMPLE=$4
             local DIRECTION=$3

              [ -d ${{OUTDIR}}/${{SAMPLE}} ] || mkdir -p ${{OUTDIR}}/${{SAMPLE}}
             grep -A 3 \
                 --no-group-separator "sample=${{SAMPLE}};" ${{DEMUX_FILE}} \
                 > "${{OUTDIR}}/${{SAMPLE}}/${{SAMPLE}}_${{DIRECTION}}.fastq" \
                 && gzip "${{OUTDIR}}/${{SAMPLE}}/${{SAMPLE}}_${{DIRECTION}}.fastq"
        }}

        export -f split_samples
        {params.program} -j 10 "split_samples 06.Split/ {input}  R1 {{}}" \
              ::: {params.samples}

        """

rule Split_reverse:
    input: rules.Demultiplex.output.rev
    output: expand("06.Split/{sample}/{sample}_R2.fastq.gz", sample=SAMPLES)
    params:
        program=config['programs_path']['parallel'],
        samples=SAMPLES
    threads: 10
    shell:
        """
        function split_samples(){{

             local OUTDIR=$1
             local DEMUX_FILE=$2
             local SAMPLE=$4
             local DIRECTION=$3

              [ -d ${{OUTDIR}}/${{SAMPLE}} ] || mkdir -p ${{OUTDIR}}/${{SAMPLE}}
             grep -A 3 \
                 --no-group-separator "sample=${{SAMPLE}};" ${{DEMUX_FILE}} \
                 > "${{OUTDIR}}/${{SAMPLE}}/${{SAMPLE}}_${{DIRECTION}}.fastq" \
                 && gzip "${{OUTDIR}}/${{SAMPLE}}/${{SAMPLE}}_${{DIRECTION}}.fastq"
        }}

        export -f split_samples
        {params.program} -j 10 "split_samples 06.Split/ {input}  R2 {{}}" \
                  ::: {params.samples}

        """



rule Count_Seqs:
    input: expand(["06.Split/{sample}/{sample}_R1.fastq.gz", "06.Split/{sample}/{sample}_R2.fastq.gz"], sample=SAMPLES)
    output: "07.Count_Seqs/seqs_stat.txt"
    params:
        program=config['programs_path']['seqkit']
    shell:
        """
        # Get the stats on the sequences using seqkit
        {params.program} stats {input} > temp.txt
           
         # Sort the sequence statistics
         (sed -n '1p' temp.txt; awk 'NR>1{{print}}' temp.txt | \
           sort -V -k1,1) > {output} \
           && rm temp.txt
        """

