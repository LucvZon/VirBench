rule create_sample_overview:
    input:
        raw_reads=expand(os.path.join(QC_DIR, "{sample}.merged.fastq"), sample=SAMPLES),
        qc_reads=expand(os.path.join(QC_DIR, "{sample}.qc.fastq"), sample=SAMPLES),
        target_reads=expand(os.path.join(READ_CLASSIFICATION_DIR, "{sample}.target_reads.fastq"), sample=SAMPLES)
    output:
        tsv = os.path.join(STATS_DIR, "sample_overview.tsv")
    log:
        os.path.join(LOG_DIR, "sample_overview.log")
    params:
        sample_names = SAMPLES
    threads:
        config["params"]["threads"]
    shell:
        """
        # 1. Initialize file with header
        echo -e "sample\\traw_count\\tqc_count\\ttarget_count\\ttarget_percent\\ttarget_avg_length" > {output.tsv}

        # 2. Convert python lists to bash arrays
        RAW_FILES=({input.raw_reads})
        QC_FILES=({input.qc_reads})
        TARGET_FILES=({input.target_reads})
        SAMPLES=({params.sample_names})

        # 3. Loop through indices
        for i in "${{!SAMPLES[@]}}"; do
            
            sample="${{SAMPLES[$i]}}"
            raw_f="${{RAW_FILES[$i]}}"
            qc_f="${{QC_FILES[$i]}}"
            target_f="${{TARGET_FILES[$i]}}"

            # --- Calculations ---

            # Raw Count (wc -l / 4)
            raw_lines=$(wc -l < "$raw_f" || echo 0)
            raw_count=$((raw_lines / 4))

            # QC Count
            qc_lines=$(wc -l < "$qc_f" || echo 0)
            qc_count=$((qc_lines / 4))

            # Target Stats (seqkit)
            # -T: Tabular, tail -n+2 skips header
            target_stats=$(seqkit stats -T -j {threads} "$target_f" 2>> {log} | tail -n+2)
            
            if [ -z "$target_stats" ]; then
                target_count=0
                target_avg_len=0
            else
                target_count=$(echo "$target_stats" | cut -f 4)
                target_avg_len=$(echo "$target_stats" | cut -f 7 | cut -d'.' -f 1)
            fi

            # Calculate Percentage
            # NOTE: Double braces {{ }} are required for awk inside Snakemake shell blocks!
            if [ "$raw_count" -gt 0 ]; then
                target_percent=$(awk -v t="$target_count" -v r="$raw_count" 'BEGIN {{printf "%.2f", (t/r)*100}}')
            else
                target_percent="0.00"
            fi

            # 4. Append row to output
            echo -e "${{sample}}\\t${{raw_count}}\\t${{qc_count}}\\t${{target_count}}\\t${{target_percent}}\\t${{target_avg_len}}" >> {output.tsv}

        done 2>> {log}
        """

rule create_assembly_overview:
    input:
        assembly_summaries = expand(
            os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_assembly_summary.csv"),
            sample=SAMPLES, assembly_type=ASSEMBLY_TYPES
        ),
        benchmark_summaries = expand(
            os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_benchmark_summary.csv"),
            sample=SAMPLES, assembly_type=ASSEMBLY_TYPES
        ),
        viral_gene_counts = expand(
            os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_total_viral_genes.csv"),
            sample=SAMPLES, assembly_type=ASSEMBLY_TYPES
        ),
        inspector_summaries = expand(
            os.path.join(STATS_DIR, "per_sample", "{assembly_type}", "{sample}_inspector_summary.csv"),
            sample=SAMPLES, assembly_type=ASSEMBLY_TYPES
        )
    output:
        tsv = os.path.join(STATS_DIR, "assembly_overview.tsv")
    log:
        os.path.join(LOG_DIR, "create_assembly_overview.log")
    run:
        import csv
        import os

        # --- Helper: Convert Minutes float to HH:MM:SS string ---
        def min_to_hms(val):
            try:
                minutes = float(val)
                seconds = int(minutes * 60)
                m, s = divmod(seconds, 60)
                h, m = divmod(m, 60)
                return f"{h:02d}:{m:02d}:{s:02d}"
            except (ValueError, TypeError):
                return "00:00:00"

        # --- Helper: Safe Float conversion for CSV reading ---
        def get_val(row, key, default="0"):
            return row.get(key, default)

        # Define Output Headers
        headers = [
            "assembly_type", "assembler", "sample", 
            "# Contigs", "Total length (bp)", "N50 (bp)", "L50 (bp)", "Largest contig (bp)",
            "Elapsed time (h:m:s)", "Maximum memory (GB)", 
            "Total viral genes", "Reads mapped (%)",
            "Substitutions", "Collapses", "Collapses (<= 5 bp)", "Collapses (> 5 bp)", "Collapses (> 50 bp)",
            "Expansions", "Expansions (<= 5 bp)", "Expansions (> 5 bp)", "Expansions (> 50 bp)"
        ]

        with open(output.tsv, 'w', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=headers, delimiter='\t')
            writer.writeheader()

            # Iterate through every combination of Type and Sample
            for a_type in ASSEMBLY_TYPES:
                for sample in SAMPLES:
                    
                    # 1. Define paths for this specific sample/type combo
                    base_dir = os.path.join(STATS_DIR, "per_sample", a_type)
                    
                    f_assembly = os.path.join(base_dir, f"{sample}_assembly_summary.csv")
                    f_benchmark = os.path.join(base_dir, f"{sample}_benchmark_summary.csv")
                    f_genes = os.path.join(base_dir, f"{sample}_total_viral_genes.csv")
                    f_inspector = os.path.join(base_dir, f"{sample}_inspector_summary.csv")

                    # 2. Load Viral Genes into a dict: {assembler: count}
                    genes_map = {}
                    if os.path.exists(f_genes):
                        with open(f_genes, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                # The rule 'count_viral_genes' uses column 'assembler'
                                key = row.get('assembler') or row.get('process')
                                genes_map[key] = row['viral_genes']

                    # 3. Load Benchmarks into a dict: {assembler: {time, mem}}
                    bench_map = {}
                    if os.path.exists(f_benchmark):
                        with open(f_benchmark, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                proc = row['process']
                                # Filter out non-assembler steps immediately
                                if proc == "classify_reads_diamond":
                                    continue
                                bench_map[proc] = {
                                    'time': row.get('Run Time (minutes)', '0'),
                                    'mem': row.get('Peak Memory (GB)', '0')
                                }

                    # 4. Load Inspector Stats into a dict: {assembler: counts}
                    insp_map = {}
                    if os.path.exists(f_inspector):
                        with open(f_inspector, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                insp_map[row['assembler']] = row

                    # 5. Read Assembly Stats (Main Loop) and Merge
                    if os.path.exists(f_assembly):
                        with open(f_assembly, 'r') as f:
                            reader = csv.DictReader(f)
                            for row in reader:
                                assembler = row['process']
                                
                                # Skip non-assemblers if they appear here too
                                if assembler == "classify_reads_diamond":
                                    continue

                                # Retrieve auxiliary data
                                b_stats = bench_map.get(assembler, {'time': '0', 'mem': '0'})
                                v_genes = genes_map.get(assembler, '0')
                                i_stats = insp_map.get(assembler, {})

                                # Prepare Output Row
                                out_row = {
                                    "assembly_type": a_type,
                                    "assembler": assembler,
                                    "sample": sample,
                                    "# Contigs": get_val(row, '# Contigs'),
                                    "Total length (bp)": get_val(row, 'Total Length (bp)'),
                                    "N50 (bp)": get_val(row, 'N50 (bp)'),
                                    "L50 (bp)": get_val(row, 'L50 (bp)', "N/A"), 
                                    "Largest contig (bp)": get_val(row, 'Largest Contig (bp)'),
                                    "Elapsed time (h:m:s)": min_to_hms(b_stats['time']),
                                    "Maximum memory (GB)": b_stats['mem'],
                                    "Total viral genes": v_genes,
                                    "Reads mapped (%)": get_val(row, 'Reads Mapped (%)'),
                                    "Substitutions": get_val(i_stats, "Substitutions", "N/A"),
                                    "Collapses": get_val(i_stats, "Collapses", "N/A"),
                                    "Collapses (<= 5 bp)": get_val(i_stats, "Collapses (<= 5 bp)", "N/A"),
                                    "Collapses (> 5 bp)": get_val(i_stats, "Collapses (> 5 bp)", "N/A"),
                                    "Collapses (> 50 bp)": get_val(i_stats, "Collapses (> 50 bp)", "N/A"),
                                    "Expansions": get_val(i_stats, "Expansions", "N/A"),
                                    "Expansions (<= 5 bp)": get_val(i_stats, "Expansions (<= 5 bp)", "N/A"),
                                    "Expansions (> 5 bp)": get_val(i_stats, "Expansions (> 5 bp)", "N/A"),
                                    "Expansions (> 50 bp)": get_val(i_stats, "Expansions (> 50 bp)", "N/A")
                                }
                                
                                writer.writerow(out_row)

rule create_accuracy_overview:
    input:
        # ONLY the valid QUAST jobs via 'get_sample_quast_comparisons' function.
        flags = expand(os.path.join(STATS_DIR, "targeted_comparisons_status", "{sample}.done"), sample=SAMPLES)
    output:
        tsv = os.path.join(STATS_DIR, "accuracy_overview.tsv")
    log:
        os.path.join(LOG_DIR, "create_accuracy_overview.log")
    run:
        import csv
        import os

        # Helper to safely extract value
        def get_val(row, key, default="0"):
            return row.get(key, default)

        headers = [
            "Species", "Assembly type", "Assembler", "Sample",
            # Basic Stats
            "Genome fraction (%)", "Duplication ratio", 
            "Total length", "Total aligned length", "Reference length",
            "NA50", "NGA50", "LGA50", "auNGA", 
            "# Unaligned contigs",
            # Errors / Misassemblies
            "# Misassemblies", "Misassembled contigs length",
            "# c. relocations", "# c. inversions",
            # Mismatches / Indels
            "# mismatches", 
            "# Indels", "Indels length", "# Indels (> 5 bp)", 
            "# indels per 100 kbp"
        ]

        with open(output.tsv, 'w', newline='') as out_f:
            writer = csv.DictWriter(out_f, fieldnames=headers, delimiter='\t')
            writer.writeheader()
            
            # Pattern: STATS_DIR/targeted_quast/{type}/{sample}/{virus}/transposed_report.tsv
            search_pattern = os.path.join(STATS_DIR, "targeted_quast", "*", "*", "*", "transposed_report.tsv")
            found_reports = glob.glob(search_pattern)

            for main_rep_path in found_reports:
                
                # Check 1: Ignore empty dummy files (created by the fallback logic)
                if os.path.getsize(main_rep_path) == 0:
                    continue

                # Parse path to get metadata
                # Path structure: .../targeted_quast/{type}/{sample}/{virus}/transposed_report.tsv
                path_parts = main_rep_path.split(os.sep)
                virus = path_parts[-2]
                sample = path_parts[-3]
                a_type = path_parts[-4]

                # Define path to the corresponding misassemblies file
                base_dir = os.path.dirname(main_rep_path)
                misasm_rep_path = os.path.join(base_dir, "contigs_reports", "transposed_report_misassemblies.tsv")

                # 1. Read Main Stats
                main_stats = {}
                with open(main_rep_path, 'r') as f:
                    reader = csv.DictReader(f, delimiter='\t')
                    for row in reader:
                        assembler = row.get('Assembly')
                        if assembler: main_stats[assembler] = row

                # 2. Read Misassembly Stats (if valid)
                misasm_stats = {}
                if os.path.exists(misasm_rep_path) and os.path.getsize(misasm_rep_path) > 0:
                    with open(misasm_rep_path, 'r') as f:
                        reader = csv.DictReader(f, delimiter='\t')
                        for row in reader:
                            assembler = row.get('Assembly')
                            if assembler: misasm_stats[assembler] = row

                # 3. Write Rows
                for assembler, m_row in main_stats.items():
                    mis_row = misasm_stats.get(assembler, {})

                    out_row = {
                        "Species": virus,
                        "Assembly type": a_type,
                        "Assembler": assembler,
                        "Sample": sample,
                        
                        # --- From transposed_report.tsv ---
                        "Genome fraction (%)": get_val(m_row, 'Genome fraction (%)'),
                        "Duplication ratio": get_val(m_row, 'Duplication ratio'),
                        "# Unaligned contigs": get_val(m_row, '# unaligned contigs'),
                        "NA50": get_val(m_row, 'NA50'),
                        
                        # New Main Stats
                        "Total length": get_val(m_row, 'Total length'),
                        "Total aligned length": get_val(m_row, 'Total aligned length'),
                        "Reference length": get_val(m_row, 'Reference length'),
                        "# indels per 100 kbp": get_val(m_row, '# indels per 100 kbp'),
                        "NGA50": get_val(m_row, 'NGA50'),
                        "LGA50": get_val(m_row, 'LGA50'),
                        "auNGA": get_val(m_row, 'auNGA'),

                        # --- From transposed_report_misassemblies.tsv ---
                        # Logic handles cases where # Misassemblies is in either file
                        "# Misassemblies": m_row.get('# misassemblies') or mis_row.get('# misassemblies', '0'),
                        "# Indels": get_val(mis_row, '# indels'),
                        
                        # Note: Keeping the 4 spaces in the key lookup as requested
                        "# Indels (> 5 bp)": get_val(mis_row, '    # indels (> 5 bp)'),

                        # New Misassembly Stats
                        "# c. relocations": get_val(mis_row, '    # c. relocations'),
                        "# c. inversions": get_val(mis_row, '    # c. inversions'),
                        "Misassembled contigs length": get_val(mis_row, 'Misassembled contigs length'),
                        "# mismatches": get_val(mis_row, '# mismatches'),
                        "Indels length": get_val(mis_row, 'Indels length')
                    }
                    writer.writerow(out_row)
