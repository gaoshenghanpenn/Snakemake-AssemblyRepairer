# Snakemake-AssemblyRepairer
Snakemake-AssemblyRepairer: A snakemake version [AssemblyRepairer](https://github.com/logsdon-lab/AssemblyRepairer.git). Error detection based on [Nucflag v1.0](https://github.com/logsdon-lab/NucFlag.git).

## Dependencies and Installation
Development environment: Linux  
Development tool: VScode
```
conda install bioconda::snakemake=9.5.0
git clone https://github.com/logsdon-lab/AssemblyRepairer.git
```

## Quick start
```
In config:
samples:
  test1:
    # coordinates_file: demo/coordinate_file.txt
    target_ref_file: demo/testasm/NA20355_chr16_haplotype2-0000183.fa
    second_assembly_file: demo/testasm/NA20355_chr16_h2tg000041l#1-37317985.fa
    hifi_reads_dir: demo/reads
    hifi_reads_suffix: hifi.fa
    nucflag_hifi_preset: hifi
    ont_reads_dir: demo/reads
    ont_reads_suffix: ont.fa
    nucflag_ont_preset: ont_r9
  
  test2:
    # coordinates_file: demo/coordinate_file.txt
    target_ref_file: demo/testasm/NA20355_chr16_haplotype2-0000183.fa
    second_assembly_file: demo/testasm/NA20355_chr16_h2tg000041l#1-37317985.fa
    hifi_reads_dir: demo/reads
    hifi_reads_suffix: hifi.fa
    nucflag_hifi_preset: hifi
    ont_reads_dir: demo/reads
    ont_reads_suffix: ont.fa
    nucflag_ont_preset: ont_r9

Run:
snakemake -p --cores 48 --configfile config/config.all.yaml --conda-frontend conda --use-conda
```

## Citing 
**Gao S, Oshima KK**, Chuang SC, Loftus M, Montanari A, Gordon DS, Human Genome Structural Variation Consortium, Human Pangenome Reference Consortium, Hsieh P, Konkel MK, Ventura M, Logsdon GA. A global view of human centromere variation and evolution. bioRxiv. 2025. p. 2025.12.09.693231. [doi:10.64898/2025.12.09.693231](https://doi.org/10.64898/2025.12.09.693231)

