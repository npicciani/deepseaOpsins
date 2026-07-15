# Modular replacement for workflow/scripts/pia.py — opsin identification
# pipeline (adapted from Speiser et al. 2010). Each function in the original
# script becomes its own rule below; the biopython/ete3 logic that doesn't
# translate cleanly to shell lives in small scripts under workflow/scripts/.
#
# Reference resources are fixed across all transcriptomes (not wildcards).

BAIT = "resources/baits0206.fasta"
REF_ALIGNMENT = "resources/all_0512_gt1_rs50_s65.aln"
REF_TREE = "resources/rep17.treefile"

# Reference taxa used to root / delimit the opsin clade in the pre-built tree
OUTGROUP_LEAF = "196_Hpoly_JEL142_contig01071"
CLADE_ANCHOR_1 = "Patiria_miniata_ops5_chaopsin"
CLADE_ANCHOR_2 = "Cavia_porcellus_H0W479_RGR_retinochrome"


rule make_blastdb:
    input:
        pep="resources/{transcriptome}.pep",
    output:
        multiext(
            "results/opsins/{transcriptome}/{transcriptome}.blastdb",
            ".phr",
            ".pin",
            ".psq",
        ),
    params:
        dbname=lambda wc, output: output[0][: -len(".phr")],
    log:
        "logs/opsins/{transcriptome}_make_blastdb.log",
    shell:
        "makeblastdb -in {input.pep} -dbtype prot -parse_seqids "
        "-out {params.dbname} &> {log}"


rule blastp_opsins:
    input:
        db=rules.make_blastdb.output,
        baits=BAIT,
    output:
        "results/opsins/{transcriptome}/{transcriptome}.opsins_blast_results",
    params:
        dbname=lambda wc, input: input.db[0][: -len(".phr")],
    threads: 8
    log:
        "logs/opsins/{transcriptome}_blastp.log",
    shell:
        "blastp -db {params.dbname} -query {input.baits} -out {output} "
        "-evalue 1e-10 -outfmt 7 &> {log}"


rule extract_blast_hits:
    input:
        "results/opsins/{transcriptome}/{transcriptome}.opsins_blast_results",
    output:
        "results/opsins/{transcriptome}/{transcriptome}.candidate_opsins_seq_list",
    shell:
        # column 2 = subject id; drop comment lines, dedupe keeping first occurrence order
        "awk -F'\\t' '!/^#/ {{print $2}}' {input} | awk '!seen[$0]++' > {output}"


rule extract_candidate_opsins:
    input:
        script="workflow/scripts/extract_seqs_by_id.py",
        ids="results/opsins/{transcriptome}/{transcriptome}.candidate_opsins_seq_list",
        pep="resources/{transcriptome}.pep",
    output:
        "results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.fasta",
    shell:
        "python {input.script} {input.ids} {input.pep} {output}"


rule align_candidate_opsins:
    input:
        opsins="results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.fasta",
        alignment=REF_ALIGNMENT,
    output:
        "results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.aln",
    log:
        "logs/opsins/{transcriptome}_mafft.log",
    shell:
        "mafft --add {input.opsins} --reorder {input.alignment} > {output} 2> {log}"


rule aln_to_phylip:
    input:
        script="workflow/scripts/fasta_to_phylip.py",
        aln="results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.aln",
    output:
        "results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.phy",
    shell:
        "python {input.script} {input.aln} {output}"


rule raxml_epa:
    input:
        phy="results/opsins/{transcriptome}/{transcriptome}.candidate_opsins.phy",
        tree=REF_TREE,
    output:
        "results/opsins/{transcriptome}/RAxML_labelledTree.Opsins",
    params:
        outdir=lambda wc: "$PWD/results/opsins/{}".format(wc.transcriptome),
    threads: 8
    log:
        "logs/opsins/{transcriptome}_raxml.log",
    shell:
        # RAxML refuses to run if files from a previous -n Opsins run are
        # still sitting in the working dir, so clear them first.
        "rm -f {params.outdir}/RAxML_*.Opsins && "
        "raxmlHPC-PTHREADS-SSE3 -f v -s {input.phy} -t {input.tree} "
        "-w {params.outdir} -m PROTGAMMAWAG -T {threads} -n Opsins &> {log}"


rule clean_raxml_tree:
    input:
        "results/opsins/{transcriptome}/RAxML_labelledTree.Opsins",
    output:
        "results/opsins/{transcriptome}/RAxML_labelledTree_Corrected.Opsins.tre",
    shell:
        r"sed -E 's/\[I[0-9]+\]//g' {input} > {output}"


rule detach_opsin_clade:
    input:
        script="workflow/scripts/detach_opsin_clade.py",
        tree="results/opsins/{transcriptome}/RAxML_labelledTree_Corrected.Opsins.tre",
    output:
        "results/opsins/{transcriptome}/RAxML_labelledTree_OpsinClade.tre",
    params:
        outgroup=OUTGROUP_LEAF,
        anchor1=CLADE_ANCHOR_1,
        anchor2=CLADE_ANCHOR_2,
    shell:
        "python {input.script} {input.tree} {params.outgroup} "
        "{params.anchor1} {params.anchor2} {output}"


rule list_opsin_queries:
    input:
        script="workflow/scripts/list_opsin_queries.py",
        tree="results/opsins/{transcriptome}/RAxML_labelledTree_OpsinClade.tre",
    output:
        "results/opsins/{transcriptome}/{transcriptome}.opsinQueryList.txt",
    shell:
        "python {input.script} {input.tree} {output}"


rule get_opsins:
    input:
        script="workflow/scripts/extract_seqs_by_id.py",
        ids="results/opsins/{transcriptome}/{transcriptome}.opsinQueryList.txt",
        pep="resources/{transcriptome}.pep",
    output:
        "results/opsins/{transcriptome}.opsins.fasta",
    shell:
        "python {input.script} {input.ids} {input.pep} {output}"
