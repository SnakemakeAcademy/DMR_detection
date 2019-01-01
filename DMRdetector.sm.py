###software


#shell.prefix("source activate mageck_vispr")
DIR_RAW = config["bismark_report"]
SAMPLE = glob_wildcards(DIR_RAW + "{sample}.bis_rep.cov.CX_report.txt").sample
print(SAMPLE)

DIR_BISREP = "s1_bis_rep/"
DIR_WINMETH = "s2_winmeth/"

############################## DMR comparison
cmp=config["comparison_list"]
DMR_TEST="s3_mer_win/";
#if not DIR_BISREP:
dict_cmp_CG = dict()
dict_cmp_CHG = dict()
dict_cmp_CHH = dict()
with open(cmp) as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        line=line.rstrip()
        ele = line.split("\t")
        key = '_'.join(ele)
        print(ele)
        dict_cmp_CG[key] = DIR_WINMETH + ele[0] + ".CG.tab," + DIR_WINMETH + ele[1] + ".CG.tab"
        dict_cmp_CHG[key] = DIR_WINMETH + ele[0] + ".CHG.tab," + DIR_WINMETH + ele[1] + ".CHG.tab"
        dict_cmp_CHH[key] = DIR_WINMETH + ele[0] + ".CHH.tab," + DIR_WINMETH + ele[1] + ".CHH.tab"


rule all:
    input:
        expand(DIR_BISREP + "{sample}.bis_rep.tab.bgz", sample=SAMPLE),
        expand(DIR_WINMETH + "win_meth_{sample}.done", sample=SAMPLE),
        expand(DMR_TEST + "win_{comparison}.{context}.tab", comparison=list(dict_cmp_CG), context = ["CG", "CHG", "CHH"]),
        expand(DMR_TEST + "win_{comparison}.{context}.fisher", comparison=list(dict_cmp_CG), context = ["CG", "CHG", "CHH"]),
        expand(DMR_TEST + "win_{comparison}.{context}.fil", comparison=list(dict_cmp_CG), context = ["CG", "CHG", "CHH"])
rule tabix:
    input:
        DIR_RAW + "{sample}.bis_rep.cov.CX_report.txt" 
    output: 
        tab = DIR_BISREP + "{sample}.bis_rep.tab.bgz",
        #tabix = DIR_BISREP + "{sample}.bis_rep.tab.bgz.tbi",
    shell:
        '''
cat {input} |bgzip > {output.tab}
#tabix has now created the myfile.bed.bgz.tbi, which is the index.
tabix -p vcf {output.tab} 
'''

rule win_meth:
    input:
        tab = DIR_BISREP + "{sample}.bis_rep.tab.bgz",
        fai = config["fai"] 
    params:
        prefix=DIR_WINMETH + "{sample}"
    output:
        DIR_WINMETH + "win_meth_{sample}.done"
    shell:
        '''
python scripts/gener_tab4DMR.py -i {input.tab} -f {input.fai} -w 200 -s 200 -d 10 -p {params.prefix} 
touch {output}
'''

rule tab3DMR_test:
    """
    Merge windows 
    """
    input:
        expand(DIR_WINMETH + "win_meth_{sample}.done", sample=SAMPLE)
    params: 
        CG = lambda wildcards: dict_cmp_CG[wildcards.comparison],
        CHG = lambda wildcards: dict_cmp_CHG[wildcards.comparison],
        CHH = lambda wildcards: dict_cmp_CHH[wildcards.comparison],
        depth = 10
    output:
        CG=DMR_TEST + "win_{comparison}.CG.tab",
        CHG=DMR_TEST + "win_{comparison}.CHG.tab",
        CHH=DMR_TEST + "win_{comparison}.CHH.tab",
        CG_test=DMR_TEST + "win_{comparison}.CG.fisher",
        CHG_test=DMR_TEST + "win_{comparison}.CHG.fisher",
        CHH_test=DMR_TEST + "win_{comparison}.CHH.fisher",
        CG_fil=DMR_TEST + "win_{comparison}.CG.fil",
        CHG_fil=DMR_TEST + "win_{comparison}.CHG.fil",
        CHH_fil=DMR_TEST + "win_{comparison}.CHH.fil",
        CG_fil=DMR_TEST + "win_{comparison}.CG.fil",
    shell:
        r'''
python scripts//mer_win_lap4DMR.py -i {params.CG} -d {params.depth} -o {output.CG}
python scripts//mer_win_lap4DMR.py -i {params.CHG} -d {params.depth} -o {output.CHG}
python scripts//mer_win_lap4DMR.py -i {params.CHH} -d {params.depth} -o {output.CHH}
Rscript scripts/DMR_fisher_BH_correction.R --input {output.CG} --output {output.CG_test}
Rscript scripts/DMR_fisher_BH_correction.R --input {output.CHG} --output {output.CHG_test}
Rscript scripts/DMR_fisher_BH_correction.R --input {output.CHH} --output {output.CHH_test}
perl -ne 'next if /#/;chomp; @aa=split; my $diff = $aa[4]/($aa[4]+$aa[5]) - $aa[7]/($aa[7]+$aa[8]); print "$_\n" if $aa[10] < 0.05 && abs $diff > 0.4' {output.CG} |sort -k1,1n -k2,2n > {output.CG_fil}
perl -ne 'next if /#/;chomp; @aa=split; my $diff = $aa[4]/($aa[4]+$aa[5]) - $aa[7]/($aa[7]+$aa[8]); print "$_\n" if $aa[10] < 0.05 && abs $diff > 0.2' {output.CHG} |sort -k1,1n -k2,2n > {output.CHG_fil}
perl -ne 'next if /#/;chomp; @aa=split; my $diff = $aa[4]/($aa[4]+$aa[5]) - $aa[7]/($aa[7]+$aa[8]); print "$_\n" if $aa[10] < 0.05 && abs $diff > 0.1' {output.CHH} |sort -k1,1n -k2,2n > {output.CHH_fil}
'''

dict_meth_tab = {}
with open(cmp) as infile:
    for line in infile:
        if line.startswith("#"):
            continue
        line=line.rstrip()
        ele = line.split("\t")
        key = '_'.join(ele)
        print(ele)
        DIR_BISREP + "{sample}.bis_rep.tab.bgz"
        dict_meth_tab[key] = DIR_BISREP + ele[0] + ".bis_rep.tab.bgz," + DIR_BISREP + ele[1] + ".bis_rep.tab.bgz"

rule DMR_res:
    input:
        CHG=DMR_TEST + "win_{comparison}.CHG.fil",
        CHH=DMR_TEST + "win_{comparison}.CHH.fil",
        CG=DMR_TEST + "win_{comparison}.CG.fil",
    params:
        meth_tab = lambda wildcards: dict_meth_tab[wildcards.comparison]
    output: 
        CG=DMR_TEST + "DMR_{comparison}.CG.tab",
        CHG=DMR_TEST + "DMR_{comparison}.CHG.tab",
        CHH=DMR_TEST + "DMR_{comparison}.CHH.tab",
        CG_meth=DMR_TEST + "DMR_{comparison}.CG.meth.tab",
        CHG_meth=DMR_TEST + "DMR_{comparison}.CHG.meth.tab",
        CHH_meth=DMR_TEST + "DMR_{comparison}.CHH.meth.tab",
    shell:
        '''
mergeBed -d 1 -i {inpout.CG} > {output.CG}
mergeBed -d 1 -i {inpout.CHG} > {output.CHG}
mergeBed -d 1 -i {inpout.CHH} > {output.CHH}
python scripts/extract_meth_lev4DMR.py --dmr {output.CG} -t {params.meth_tab} -c CG --methdiff 0.4 > {output.CG_meth}
python scripts/extract_meth_lev4DMR.py --dmr {output.CHG} -t {params.meth_tab} -c CHG --methdiff 0.4 > {output.CHG_meth}
python scripts/extract_meth_lev4DMR.py --dmr {output.CHH} -t {params.meth_tab} -c CHH --methdiff 0.4 > {output.CHH_meth}
'''
