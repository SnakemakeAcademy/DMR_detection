DIR_RAW = config["raw_data_dir"]
SAMPLE = glob_wildcards(DIR_RAW + "/capbs_{sample}.bed.gz").sample
print(SAMPLE)
DIR_FQ = DIR_RAW
TAB4VIEWBS_FILE = "s1_tab4ViewBS/"
TAB4DMR_FILE = "s1_tab4DMR/"
TAB4PHI_FILE = "s1_tab4phi/"
rule all:
    input: 
        expand(TAB4VIEWBS_FILE + "capbs_{sample}.tab", sample=SAMPLE),
        #expand(TAB4DMR_FILE + "capbs_{sample}.{context}.tab", sample=SAMPLE, context = ["CG", "CHG", "CHH"]),
        expand(TAB4DMR_FILE + "capbs4view_{sample}.{context}.tab", sample=SAMPLE, context = ["CG", "CHG", "CHH"]),
        expand(TAB4DMR_FILE + "capbs4view_{sample}.3context.tab.gz", sample=SAMPLE),
        expand(TAB4PHI_FILE + "capbs_{sample}.{context}.tab", sample=SAMPLE, context = ["CG", "CHG", "CHH"]),

rule gener_tab4viewbs:
    input:
        tab=DIR_RAW + "capbs_{sample}.bed.gz"
    output:
        TAB4VIEWBS_FILE + "capbs_{sample}.tab",
    shell:
        '''
python ../../scripts/extra_infor4viewbs_capbs.py -i {input.tab} -d 1 -o {output}
'''

rule gener_tab4DMR:
    input:
        tab=DIR_RAW + "capbs_{sample}.bed.gz"
    params:
        prefix=TAB4DMR_FILE + "capbs_{sample}",
        prefix1=TAB4DMR_FILE + "capbs4view_{sample}",
    output:
        #TAB4DMR_FILE + "capbs_{sample}.CG.tab",
        #TAB4DMR_FILE + "capbs_{sample}.CHG.tab",
        #TAB4DMR_FILE + "capbs_{sample}.CHH.tab",
        TAB4DMR_FILE + "capbs4view_{sample}.CG.tab",
        TAB4DMR_FILE + "capbs4view_{sample}.CHG.tab",
        TAB4DMR_FILE + "capbs4view_{sample}.CHH.tab",
        TAB4DMR_FILE + "capbs4view_{sample}.3context.tab.gz",
    shell:
        '''
#python ../../scripts/gener_tab3DMR.py -i {input.tab} -f ../../data/seqcap_new/Zea_mays.AGPv4.dna.chrom.fa.fai -w 50 -s 50 -p {params.prefix}
python ../../scripts/extr_info4viewbs_capbs.py -f {params.prefix}.CG.tab --context CG --depth 1 -o {params.prefix1}.CG.tab
python ../../scripts/extr_info4viewbs_capbs.py -f {params.prefix}.CHG.tab --context CHG --depth 1 -o {params.prefix1}.CHG.tab
python ../../scripts/extr_info4viewbs_capbs.py -f {params.prefix}.CHH.tab --context CHH --depth 1 -o {params.prefix1}.CHH.tab
cat  {params.prefix1}.CG.tab {params.prefix1}.CHG.tab {params.prefix1}.CHH.tab |sort -k1,1 -k2,2n > {params.prefix1}.3context.tab
bgzip {params.prefix1}.3context.tab
tabix -p vcf {params.prefix1}.3context.tab.gz
'''

rule gener_tab4phi:
    input:
        cg = TAB4DMR_FILE + "capbs_{sample}.CG.tab",
        chg= TAB4DMR_FILE + "capbs_{sample}.CHG.tab",
        chh= TAB4DMR_FILE + "capbs_{sample}.CHH.tab",
    params:
        prefix=TAB4DMR_FILE + "capbs_{sample}"
    output:
        cg=TAB4PHI_FILE + "capbs_{sample}.CG.tab",
        chg=TAB4PHI_FILE + "capbs_{sample}.CHG.tab",
        chh=TAB4PHI_FILE + "capbs_{sample}.CHH.tab",
    shell:
        '''
../../software/bedtools2/bin//intersectBed -a {input.cg} -b ../../data/maize_capture_region/capseqtargets_info_phasRNA_v4.tsv -wa > {output.cg}
../../software/bedtools2/bin//intersectBed -a {input.chg} -b ../../data/maize_capture_region/capseqtargets_info_phasRNA_v4.tsv -wa > {output.chg}
../../software/bedtools2/bin//intersectBed -a {input.chh} -b ../../data/maize_capture_region/capseqtargets_info_phasRNA_v4.tsv -wa > {output.chh}
'''

