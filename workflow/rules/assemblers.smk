# File: workflow/rules/assemblers.smk
# This file contains all rules related to genomic assembly.
# It expects variables like ASSEMBLERS_CONFIG, READ_CLASSIFICATION_DIR,
# ASSEMBLY_DIR, LOG_DIR, and BENCH_DIR to be defined in the main Snakefile.

# --- ASSEMBLY RULES ---

if ASSEMBLERS_CONFIG.get("metaflye", False):
    rule assemble_metaflye:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "metaflye")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "metaflye", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "metaflye", "{sample}.tsv")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "metaflye", "{sample}.log")
        shell:
            """
            # Wrap the flye command in parentheses and use '||' for a fallback.
            # The '&>' redirects both stdout and stderr to the log file.

            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (flye --nano-hq {input} --meta --min-overlap 1000 -o {output.dir} --deterministic --threads {threads} &> {log}) \
            || \
            (echo "Flye failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("penguin", False):
    rule assemble_penguin:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "penguin", "contigs.fasta"),
            tmp_dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "penguin", "temp_files")),
            bench=os.path.join(BENCH_DIR, "primary", "penguin", "{sample}.tsv")
        params:
            min_len=config["params"]["penguin_min_contig_len"],
            min_id=config["params"]["penguin_min_seq_id"],
            min_overlap=config["params"]["penguin_min_overlap"]
        threads:
            config["params"]["threads"]		
        log:
            os.path.join(LOG_DIR, "primary", "penguin", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (penguin guided_nuclassemble {input} {output.fasta} {output.tmp_dir} \
            --min-contig-len {params.min_len} --min-seq-id {params.min_id} --min-aln-len {params.min_overlap} \
            --clust-min-seq-id 0.95 --clust-min-cov 0.85 --threads {threads} &> {log}) \
            || \
            (echo "PenguiN failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("raven", False):
    rule assemble_raven:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "raven", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "raven", "{sample}.tsv")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "raven", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (raven --threads {threads} -u 1000 --identity 0.95 -p 2 {input} > {output.fasta} 2> {log}
            rm raven.cereal) \
            || \
            (echo "Raven failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("canu", False):
    rule assemble_canu:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "canu")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "canu", "canu_assembly.contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "canu", "{sample}.tsv")
        params:
            genome_size=config["params"]["canu_genome_size"],
            # A good rule of thumb: 4GB of memory per thread for Canu
            memory=lambda wildcards, threads: threads * 4
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "canu", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (canu -p canu_assembly \
                -d {output.dir} \
                genomeSize={params.genome_size} \
                maxThreads={threads} \
                maxMemory={params.memory} \
                minReadLength=400 \
                minOverlapLength=250 \
                contigFilter="2 1000 1.0 0.5 5" \
                -nanopore {input} \
                useGrid=false 2> {log}) \
            || \
            (echo "Canu failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("myloasm", False):
    rule assemble_myloasm:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "myloasm")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "myloasm", "assembly_primary.fa"),
            bench=os.path.join(BENCH_DIR, "primary", "myloasm", "{sample}.tsv")
        params:
            min_reads=config["params"]["myloasm_min_reads"],
            min_overlap=config["params"]["myloasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "myloasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (myloasm {input} \
                -o {output.dir} \
                --quality-value-cutoff 95 \
                --absolute-coverage-threshold {params.min_reads} \
                --min-ol {params.min_overlap} \
                -t {threads} 2> {log}) \
            || \
            (echo "Myloasm failed for sample {wildcards.sample}, creating empty output." >> {log} && \
             mkdir -p {output.dir} && \
             touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("metamdbg", False):
    rule assemble_metamdbg:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "metamdbg")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "metamdbg", "contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "metamdbg", "{sample}.tsv")
        params:
            min_overlap=config["params"]["metamdbg_min_overlap"],
            min_id=config["params"]["metamdbg_min_seq_id"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "metamdbg", "{sample}.log")
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

if ASSEMBLERS_CONFIG.get("wtdbg2", False):
    rule assemble_wtdbg2:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "wtdbg2")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "wtdbg2", "contigs.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "wtdbg2", "{sample}.tsv")
        params:
            min_read_length=config["params"]["wtdbg2_min_read_len"],
            min_contig_length=config["params"]["wtdbg2_min_contig_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "wtdbg2", "{sample}.log")
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
            -l 250 \
            -e 5 \
            --node-len 512 \
            --ctg-min-length {params.min_contig_length} 2> {log}

            wtpoa-cns -t {threads} -i {output.dir}/dbg.ctg.lay.gz -fo {output.fasta}) \
            || \
            (echo "Wtdbg2 failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("shasta", False):
    rule assemble_shasta:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "shasta")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "shasta", "Assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "shasta", "{sample}.tsv")
        params:
            min_read_length=config["params"]["shasta_min_read_len"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "shasta", "{sample}.log")
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
            --MarkerGraph.minCoverage 5 \
            --MarkerGraph.minEdgeCoverage 5 \
            --Kmers.probability 0.25 \
            --Align.maxSkip 10 \
            --Align.maxDrift 10 >> {log}) \
            || \
            (echo "Shasta failed for sample {wildcards.sample}, creating empty output." >> {log} && \
            touch {output.fasta})
            '
            """

if ASSEMBLERS_CONFIG.get("miniasm", False):
    rule assemble_miniasm:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "miniasm")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "miniasm", "final_assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "miniasm", "{sample}.tsv")
        params:
            min_overlap=config["params"]["miniasm_min_overlap"]
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "miniasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            # Create draft assembly graph
            minimap2 -t {threads} -K 1m -x ava-ont {input} {input} \
            | miniasm -s {params.min_overlap} -c 5 -i 0.95 -f {input} - > {output.dir}/raw_assembly.gfa 2> {log}

            # Convert graph to fasta
            gfatools gfa2fa {output.dir}/raw_assembly.gfa > {output.dir}/raw_assembly.fasta

            # --- ROBUSTNESS CHECK ---
            # Check if the raw assembly FASTA is empty. The '-s' flag checks if a file has a size greater than zero.
            if [ -s {output.dir}/raw_assembly.fasta ]; then
                # If the file is NOT empty, proceed with polishing
                echo "Miniasm produced contigs. Proceeding with Racon polishing." >> {log}

                # First polishing round
                minimap2 -t {threads} -x map-ont {output.dir}/raw_assembly.fasta {input} \
                | gzip --fast > {output.dir}/polished_overlaps_1.paf.gz 2>> {log}
                racon -t {threads} {input} {output.dir}/polished_overlaps_1.paf.gz {output.dir}/raw_assembly.fasta > {output.dir}/polished_assembly_1.fasta 2>> {log}
                
                # Second polishing round
                minimap2 -t {threads} -x map-ont {output.dir}/polished_assembly_1.fasta {input} \
                | gzip --fast > {output.dir}/polished_overlaps_2.paf.gz 2>> {log}
                racon -t {threads} {input} {output.dir}/polished_overlaps_2.paf.gz {output.dir}/polished_assembly_1.fasta > {output.fasta} 2>> {log}

            else
                # If the file IS empty, skip polishing and create an empty final file
                echo "Miniasm produced no contigs for sample {wildcards.sample}. Skipping polishing." >> {log}
                touch {output.fasta}
            fi

            # --- CLEANUP MASSIVE INTERMEDIATE FILES ---
            rm -f {output.dir}/*.paf.gz {output.dir}/*.gfa {output.dir}/raw_assembly.fasta {output.dir}/polished_assembly_1.fasta
            '
            """

if ASSEMBLERS_CONFIG.get("hifiasm", False):
    rule assemble_hifiasm:
        input:
            os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq")
        output:
            dir=directory(os.path.join(ASSEMBLY_DIR, "{sample}", "hifiasm")),
            fasta=os.path.join(ASSEMBLY_DIR, "{sample}", "hifiasm", "assembly.fasta"),
            bench=os.path.join(BENCH_DIR, "primary", "hifiasm", "{sample}.tsv")
        threads:
            config["params"]["threads"]
        log:
            os.path.join(LOG_DIR, "primary", "hifiasm", "{sample}.log")
        shell:
            """
            /usr/bin/time -f "s\\tmax_rss\\tmean_load\\n%e\\t%M\\t%P" -o {output.bench} \
            bash -c '
            (hifiasm -o {output.dir}/asm \
            {input} \
            --ont \
            --rl-cut 400 \
            --chem-c 1 \
            --chem-f 100 \
            -r 3 \
            -a 4 \
            -D 200 \
            --ctg-n 5 \
            -t {threads}) \

            gfatools gfa2fa {output.dir}/asm.bp.p_ctg.gfa > {output.fasta}
            || \
            (echo "Hifiasm failed for sample {wildcards.sample}, creating empty output. Check log for details." >> {log} && \
             touch {output.fasta})
            '
            """
