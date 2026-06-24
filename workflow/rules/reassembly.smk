# This snakefile performs reassembly with each active assembler

# --- Helper function to get correct assembly output file ---
def get_reassembly_fasta(wildcards):
    assembler = wildcards.assembler
    if assembler == "metaflye":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "penguin":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "raven":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    elif assembler == "canu":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "canu_assembly.contigs.fasta")
    elif assembler == "myloasm":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly_primary.fa")
    elif assembler == "metamdbg":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "wtdbg2":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "contigs.fasta")
    elif assembler == "shasta":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "Assembly.fasta")
    elif assembler == "miniasm":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "final_assembly.fasta")
    elif assembler == "hifiasm":
        return os.path.join(REASSEMBLY_DIR, wildcards.sample, assembler, "assembly.fasta")
    # Add other assemblers here if needed in the future

# Place assembly rules here
if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_metaflye:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "metaflye")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "metaflye", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "metaflye", "{sample}.tsv")
        params:
            min_overlap=config["params"]["metaflye_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "metaflye", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (flye --nano-hq {input} --meta --min-overlap {params.min_overlap} -o {output.dir} --deterministic --threads {threads} &> {log}) \
            || \
            (echo "Flye failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_penguin:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "penguin", "contigs.fasta"),
            tmp_dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "penguin", "temp_files")),
            bench=os.path.join(BENCH_DIR, "secondary", "penguin", "{sample}.tsv")
        params:
            min_len=config["params"]["penguin_min_contig_len"],
            min_id=config["params"]["penguin_min_seq_id"],
            min_overlap=config["params"]["penguin_min_overlap"],
            min_clust_id=config["params"]["penguin_min_clust_seq_id"],
            min_clust_cov=config["params"]["penguin_min_clust_cov"]
        threads:
            config["params"]["threads"]		
        log:
            os.path.join(LOG_DIR, "secondary", "penguin", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (penguin guided_nuclassemble {input} {output.fasta} {output.tmp_dir} \
            --min-contig-len {params.min_len} --min-seq-id {params.min_id} --min-aln-len {params.min_overlap} \
            --clust-min-seq-id {params.min_clust_id} --clust-min-cov {params.min_clust_cov} --threads {threads} &> {log}) \
            || \
            (echo "PenguiN failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_raven:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "raven", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "raven", "{sample}.tsv")
        params:
            min_contig_length=config["params"]["raven_min_contig_len"],
            min_seq_id=config["params"]["raven_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "raven", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (raven --threads {threads} -u {params.min_contig_length} --identity {params.min_seq_id} -p 2 {input} > {output.fasta} 2> {log}
            rm raven.cereal) \
            || \
            (echo "Raven failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_canu:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "canu")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "canu", "canu_assembly.contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "canu", "{sample}.tsv")
        params:
            genome_size=config["params"]["canu_genome_size"],
            min_read_length=config["params"]["canu_min_read_len"],
            min_overlap=config["params"]["canu_min_overlap"],
            min_contig_length=config["params"]["canu_min_contig_len"],
            min_coverage=config["params"]["canu_min_cov"],
            # A rule of thumb: 4GB of memory per thread for Canu
            memory=lambda wildcards, threads: threads * 4
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "canu", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (canu -p canu_assembly \
            -d {output.dir} \
            genomeSize={params.genome_size} \
            maxThreads={threads} \
            maxMemory={params.memory} \
            minReadLength={params.min_read_length} \
            minOverlapLength={params.min_overlap} \
            contigFilter="2 {params.min_contig_length} 1.0 0.5 {params.min_coverage}" \
            -nanopore {input} \
            useGrid=false 2> {log}) \
            || \
            (echo "Canu failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_myloasm:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "myloasm")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "myloasm", "assembly_primary.fa"),
            bench=os.path.join(BENCH_DIR, "secondary", "myloasm", "{sample}.tsv")
        params:
            min_coverage=config["params"]["myloasm_min_cov"],
            min_overlap=config["params"]["myloasm_min_overlap"],
            min_seq_id=config["params"]["myloasm_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "myloasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (myloasm {input} \
            -o {output.dir} \
            --quality-value-cutoff {params.min_seq_id} \
            --absolute-coverage-threshold {params.min_coverage} \
            --min-ol {params.min_overlap} \
            -t {threads} 2> {log}) \
            || \
            (echo "Myloasm failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_metamdbg:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "metamdbg")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "metamdbg", "contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "metamdbg", "{sample}.tsv")
        params:
            min_overlap=config["params"]["metamdbg_min_overlap"],
            min_id=config["params"]["metamdbg_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "metamdbg", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (metaMDBG asm --in-ont {input} \
            --out-dir {output.dir} \
            --min-read-overlap {params.min_overlap} \
            --min-read-identity {params.min_id} \
            --density-correction 0.1 \
            --density-assembly 0.05 \
            --threads {threads} 2> {log}

            gzip --decompress -c {output.dir}/contigs.fasta.gz > {output.fasta}) \
            || \
            (echo "metaMDBG failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_wtdbg2:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "wtdbg2")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "wtdbg2", "contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "wtdbg2", "{sample}.tsv")
        params:
            min_read_length=config["params"]["wtdbg2_min_read_len"],
            min_contig_length=config["params"]["wtdbg2_min_contig_len"],
            min_overlap=config["params"]["wtdbg2_min_overlap"],
            min_coverage=config["params"]["wtdbg2_min_cov"],
            min_node_length=config["params"]["wtdbg2_min_node_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "wtdbg2", "{sample}.log")
        shell:
            """
            mkdir -p {output.dir}

            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (wtdbg2 \
            -i {input} \
            -o {output.dir}/dbg \
            -t {threads} \
            -x preset4 \
            -L {params.min_read_length} \
            -l {params.min_overlap} \
            -e {params.min_coverage} \
            --node-len {params.min_node_length} \
            --ctg-min-length {params.min_contig_length} 2> {log}

            wtpoa-cns -t {threads} -i {output.dir}/dbg.ctg.lay.gz -fo {output.fasta}) \
            || \
            (echo "Wtdbg2 failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_shasta:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "shasta")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "shasta", "Assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "shasta", "{sample}.tsv")
        params:
            min_read_length=config["params"]["shasta_min_read_len"],
            min_coverage=config["params"]["shasta_min_cov"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "shasta", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (shasta \
            --input {input} \
            --threads {threads} \
            --config Nanopore-R10-Fast-Nov2022 \
            --Reads.minReadLength {params.min_read_length} \
            --assemblyDirectory {output.dir} \
            --Align.minAlignedMarkerCount 25 \
            --Align.minAlignedFraction 0.80 \
            --MarkerGraph.minCoverage {params.min_coverage} \
            --MarkerGraph.minEdgeCoverage {params.min_coverage} \
            --Kmers.probability 0.25 \
            --Align.maxSkip 10 \
            --Align.maxDrift 10 >> {log}) \
            || \
            (echo "Shasta failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_miniasm:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "miniasm")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "miniasm", "final_assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "miniasm", "{sample}.tsv")
        params:
            min_overlap=config["params"]["miniasm_min_overlap"],
            min_coverage=config["params"]["miniasm_min_cov"],
            min_seq_id=config["params"]["miniasm_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "miniasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            # Create draft assembly graph
            minimap2 -t {threads} -K 1m -x ava-ont {input} {input} \
            | miniasm -s {params.min_overlap} -c {params.min_coverage} -i {params.min_seq_id} -f {input} - > {output.dir}/raw_assembly.gfa 2> {log}

            # --- ROBUSTNESS CHECK ---
            if [ -s {output.dir}/raw_assembly.gfa ]; then
                # Convert graph to fasta
                gfatools gfa2fa {output.dir}/raw_assembly.gfa > {output.fasta}
            else
                echo "Miniasm produced no contigs for sample {wildcards.sample}. Creating empty output." >> {log}
                touch {output.fasta}
            fi
            
            # --- CLEANUP MASSIVE INTERMEDIATE FILES ---
            rm -f {output.dir}/*.paf {output.dir}/*.gfa
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule reassemble_hifiasm:
        input:
            os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            dir=directory(os.path.join(REASSEMBLY_DIR, "{sample}", "hifiasm")),
            fasta=os.path.join(REASSEMBLY_DIR, "{sample}", "hifiasm", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "secondary", "hifiasm", "{sample}.tsv")
        params:
            min_read_length=config["params"]["hifiasm_min_read_len"],
            min_tip_filter=config["params"]["hifiasm_tip_contig_filter"],
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "secondary", "hifiasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (hifiasm -o {output.dir}/asm \
            {input} \
            --ont \
            --rl-cut {params.min_read_length} \
            --chem-c 1 \
            --chem-f 100 \
            -r 3 \
            -a 4 \
            -D 200 \
            --ctg-n {params.min_tip_filter} \
            -t {threads}) 2> {log}

            gfatools gfa2fa {output.dir}/asm.bp.p_ctg.gfa > {output.fasta} \
            || \
            (echo "Hifiasm failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output.fasta})
            '
            """

# Map original contigs to reassembled contigs
if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule map_contigs_to_reassembly:
        message:
            "Map all contigs from the original (all assemblers) to a specific reassembly"
        input:
            reassembled=get_reassembly_fasta,
            combined=os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "combined_contigs_to_reassembly.bam")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "contigs_to_reassembly", "{sample}_{assembler}.log")
        shell:
            """
            # I don't want to grab get_assembly_fasta, I want to grab combined.fasta
            minimap2 -ax asm20 -t {threads} {input.reassembled} {input.combined} \
            | samtools sort -@ {threads} --output-fmt BAM -o {output}
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule extract_and_cluster:
        message:
            "Extracting unmapped contigs and clustering them with MMseqs2"
        input:
            bam=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "combined_contigs_to_reassembly.bam"),
            combined=os.path.join(ASSEMBLY_DIR, "{sample}", "combined", "combined_assemblies.fasta")
        output:
            read_ids=temp(os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "{sample}_unmapped_read_ids.txt")),
            unmapped=temp(os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "{sample}_unmapped_contigs.fasta")),
            rep_seq=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "cluster_rep_seq.fasta"),
            cluster_tsv=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "cluster_cluster.tsv"),
            bench=os.path.join(BENCH_DIR, "final", "{assembler}", "{sample}.tsv")
        params:
            out_prefix=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "cluster"),
            tmp_dir=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "tmp")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "clusters", "{sample}_{assembler}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            # 1. Extract unmapped contigs
            # -f 4 gets unmapped reads
            samtools view -f 4 {input.bam} | cut -f1 | sort -u > {output.read_ids}

            # Extract sequences (seqkit grep needs exact ID matches)
            seqkit grep -f {output.read_ids} {input.combined} -o {output.unmapped}
            
            # 2. Cluster unmapped contigs
            # Ensure the temp dir exists
            mkdir -p {params.tmp_dir}

            # Command syntax: mmseqs easy-cluster <input> <output_prefix> <tmp_dir>
            mmseqs easy-cluster {output.unmapped} {params.out_prefix} {params.tmp_dir} \
            --min-seq-id 0.9 -c 0.8 --cov-mode 1 --remove-tmp-files 1 >> {log} 2>&1
            '
            """

if REASSEMBLY_CONFIG.get("reassemble_contigs", False):
    rule merge_clustered_and_reassembly:
        message:
            "Combine unmapped clustered contigs with reassembled contigs"
        input:
            clusters=os.path.join(CLUSTER_DIR, "{sample}", "{assembler}", "cluster_rep_seq.fasta"),
            reassembly=get_reassembly_fasta
        output:
            final_fasta=os.path.join(FINAL_DIR, "{sample}", "{assembler}", "final_assembly.fasta")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "contigs_to_reassembly", "{sample}_{assembler}.log")
        shell:
            """
            cat {input.clusters} {input.reassembly} > {output.final_fasta}
            """
