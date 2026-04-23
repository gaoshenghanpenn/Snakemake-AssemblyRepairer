import os

# ==============================================================================
# 1. Module Config & Global Variables
# ==============================================================================
RUNNING_TYPE = config.get("running_type", "nucflag_run")
OUTDIR = config.get("outdir", "results")
LOGS_DIR = config.get("logs_dir", "logs")
BMKS_DIR = config.get("bmks_dir", "bmks")

THREADS_ALN = int(config.get("threads", 1))
MEM_ALN = str(config.get("mem_aln", "10GB"))

nucflag_filter = config.get("nucflag_filter", "")
merge_distance = config.get("merge_distance", 500000)

# ==============================================================================
# Helper Functions: 动态获取当前 {sample} 的配置
# ==============================================================================
def get_sample_val(w, key, default=None):
    return config.get("samples", {}).get(w.sample, {}).get(key, default)

def get_ref(w):
    if "ref_file_template" in config:
        return config["ref_file_template"].format(sample=w.sample)
    if "ref_file_key" in config:
        return get_sample_val(w, config["ref_file_key"])
    return get_sample_val(w, "target_ref_file") 

def get_coords(w):
    if "coords_file_template" in config:
        return [config["coords_file_template"].format(sample=w.sample)]
    
    if "coords_file_key" in config:
        key = config["coords_file_key"]
        c = get_sample_val(w, key)
        if not c:
            c = config.get(key)
        if c:
            return [c]
    return []

# ==============================================================================
# 2. HiFi Part
# ==============================================================================
def get_hifi_dir(w):
    return get_sample_val(w, "hifi_reads_dir")

def get_hifi_suffix(w):
    return get_sample_val(w, "hifi_reads_suffix", "hifi.fa")

def get_hifi_ids(w):
    d = get_hifi_dir(w)
    s = get_hifi_suffix(w)
    if not d or not os.path.exists(d): 
        return []
    return [f for f in os.listdir(d) if f.endswith(s)]

def get_hifi_aln_to_asm(w):
    ids = get_hifi_ids(w)
    if not ids:
        raise FileNotFoundError(f"No hifi read files found for sample {w.sample}")
    return expand(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "temp_bam", "{id}.hifi.bam"), sample=w.sample, id=ids)

rule align_hifi_reads_to_asm:
    input:
        asm=get_ref,
        reads=lambda w: os.path.join(get_hifi_dir(w), str(w.id))
    output:
        bam=temp(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "temp_bam", "{id}.hifi.bam"))
    threads: THREADS_ALN
    resources:
        mem=MEM_ALN,
        sort_mem="4G"
    params:
        aligner_opts="-ax lr:hqae --eqx"
    conda:
        "../envs/env.yaml"
    log:
        os.path.join(LOGS_DIR, "{sample}", RUNNING_TYPE + "_align_{id}_hifi_reads_to_asm.log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", RUNNING_TYPE + "_align_{id}_hifi_reads_to_asm.tsv")
    shell:
        """
        minimap2 {params.aligner_opts} -t {threads} {input.asm} {input.reads} | \
        samtools view -F 2308 -u - | \
        samtools sort -m {resources.sort_mem} -@ {threads} -o {output.bam} 2> {log}
        """

rule merge_hifi_read_asm_alignments:
    input:
        get_hifi_aln_to_asm
    output:
        alignment=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag.hifi.bam"),
        alignment_idx=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag.hifi.bam.bai")
    threads: THREADS_ALN
    resources:
        mem=MEM_ALN
    conda:
        "../envs/env.yaml"
    log:
        os.path.join(LOGS_DIR, "{sample}", RUNNING_TYPE + "_merge_all_hifi_reads.log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", RUNNING_TYPE + "_merge_all_hifi_reads.tsv")
    shell:
        """
        samtools merge -@ {threads} -o {output.alignment} {input} 2> {log}
        samtools index -@ {threads} {output.alignment} 2>> {log}
        """

# ==============================================================================
# 3. ONT Part
# ==============================================================================
def get_ont_dir(w):
    return get_sample_val(w, "ont_reads_dir")

def get_ont_suffix(w):
    return get_sample_val(w, "ont_reads_suffix", "ont.fa")

def get_ont_ids(w):
    d = get_ont_dir(w)
    s = get_ont_suffix(w)
    if not d or not os.path.exists(d): 
        return []
    return [f for f in os.listdir(d) if f.endswith(s)]

def get_ont_aln_to_asm(w):
    ids = get_ont_ids(w)
    if not ids:
        return []
    return expand(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "temp_bam", "{id}.ont.bam"), sample=w.sample, id=ids)

rule align_ont_reads_to_asm:
    input:
        asm=get_ref,
        reads=lambda w: os.path.join(get_ont_dir(w), str(w.id))
    output:
        bam=temp(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "temp_bam", "{id}.ont.bam"))
    threads: THREADS_ALN
    resources:
        mem=MEM_ALN,
        sort_mem="4G"
    params:
        aligner_opts="-ax lr:hqae --eqx"
    conda:
        "../envs/env.yaml"
    log:
        os.path.join(LOGS_DIR, "{sample}", RUNNING_TYPE + "_align_{id}_ont_reads_to_asm.log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", RUNNING_TYPE + "_align_{id}_ont_reads_to_asm.tsv")
    shell:
        """
        minimap2 {params.aligner_opts} -t {threads} {input.asm} {input.reads} | \
        samtools view -F 2308 -u - | \
        samtools sort -m {resources.sort_mem} -@ {threads} -o {output.bam} 2> {log}
        """

rule merge_ont_read_asm_alignments:
    input:
        get_ont_aln_to_asm
    output:
        alignment=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag.ont.bam"),
        alignment_idx=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag.ont.bam.bai")
    threads: THREADS_ALN
    resources:
        mem=MEM_ALN
    conda:
        "../envs/env.yaml"
    log:
        os.path.join(LOGS_DIR, "{sample}", RUNNING_TYPE + "_merge_all_ont_reads.log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", RUNNING_TYPE + "_merge_all_ont_reads.tsv")
    shell:
        """
        samtools merge -@ {threads} -o {output.alignment} {input} 2> {log}
        samtools index -@ {threads} {output.alignment} 2>> {log}
        """

# ==============================================================================
# 4. NucFlag Execution
# ==============================================================================
rule run_hifi_nucflag:
    input:
        bam_file=rules.merge_hifi_read_asm_alignments.output.alignment,
        regions=get_coords, 
        asm=get_ref,
    output:
        plot_dir=directory(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag_hifi_plots")),
        misassemblies=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "hifi_misassemblies.bed"),
        asm_status=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "hifi_status.bed"),
    params:
        regions_arg=lambda w, input: f"-b {input.regions}" if len(input.regions) > 0 else "",
        nucflag_hifi_preset=lambda w: get_sample_val(w, "nucflag_hifi_preset", "hifi"),
        nucflag_filter=nucflag_filter,
    threads: THREADS_ALN
    conda:
        "../envs/env.yaml"
    resources:
        mem=MEM_ALN
    log:
        os.path.join(LOGS_DIR, "{sample}", "run_hifi_nucflag." + RUNNING_TYPE + ".log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", "run_hifi_nucflag." + RUNNING_TYPE + ".tsv")
    shell:
        """
        nucflag call \
            -i {input.bam_file} \
            -d {output.plot_dir} \
            -f {input.asm} \
            -o {output.misassemblies} \
            -s {output.asm_status} \
            -t {threads} \
            -p {threads} \
            -x {params.nucflag_hifi_preset} \
            {params.regions_arg} --overlap_calls 2> {log}
        """

rule run_ont_nucflag:
    input:
        bam_file=rules.merge_ont_read_asm_alignments.output.alignment,
        regions=get_coords,
        asm=get_ref,
    output:
        plot_dir=directory(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "nucflag_ont_plots")),
        misassemblies=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "ont_misassemblies.bed"),
        asm_status=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "ont_status.bed"),
    params:
        regions_arg=lambda w, input: f"-b {input.regions}" if len(input.regions) > 0 else "",
        nucflag_ont_preset=lambda w: get_sample_val(w, "nucflag_ont_preset", "ont_r9"),
        nucflag_filter=nucflag_filter,
    threads: THREADS_ALN
    conda:
        "../envs/env.yaml"
    resources:
        mem=MEM_ALN
    log:
        os.path.join(LOGS_DIR, "{sample}", "run_ont_nucflag." + RUNNING_TYPE + ".log")
    benchmark:
        os.path.join(BMKS_DIR, "{sample}", "run_ont_nucflag." + RUNNING_TYPE + ".tsv")
    shell:
        """
        nucflag call \
            -i {input.bam_file} \
            -d {output.plot_dir} \
            -f {input.asm} \
            -o {output.misassemblies} \
            -s {output.asm_status} \
            -t {threads} \
            -p {threads} \
            -x {params.nucflag_ont_preset} \
            {params.regions_arg} --overlap_calls 2> {log}
        """

# ==============================================================================
# Merge Errors
# ==============================================================================
def get_nucflag_ont_input(w):
    if get_ont_ids(w):
        return rules.run_ont_nucflag.output.misassemblies
    return []

rule merge_nucflag_errors:
    input:
        hifi_misassemblies=rules.run_hifi_nucflag.output.misassemblies,
        ont_misassemblies=get_nucflag_ont_input,
        regions_file=get_coords,
    output:
        merged_dir=directory(os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "merge_errors")),
        merged_errors=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "merge_errors.bed"),
        merged_error_regions=os.path.join(OUTDIR, "{sample}", RUNNING_TYPE, "merge_errors_regions.bed"),
    params:
        merge_distance=merge_distance,
        nucflag_filter=nucflag_filter,
        asm=get_ref,
        ont_arg=lambda w, input: f"-om {input.ont_misassemblies}" if input.ont_misassemblies else "",
        regions_arg=lambda w, input: f"-b {input.regions_file}" if len(input.regions_file) > 0 else "",
    conda:
        "../envs/env.yaml"
    log:
        os.path.join(LOGS_DIR, "{sample}", "merge_nucflag_errors." + RUNNING_TYPE + ".log")
    shell:
        """
        mkdir -p {output.merged_dir}
        python workflow/script/mergeError.py \
            -hm {input.hifi_misassemblies} \
            {params.ont_arg} \
            {params.regions_arg} \
            -f "{params.nucflag_filter}" \
            -m {params.merge_distance} \
            -r {params.asm} \
            -o {output.merged_dir} \
            -or {output.merged_error_regions} \
            -of {output.merged_errors} 2> {log}
        """

rule all_run_nucflag:
    input:
        rules.merge_nucflag_errors.output
    default_target: True