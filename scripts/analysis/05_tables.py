import json,glob,statistics,re
K="/NatureUsers/nsinyavskiy/metaviroms_kgt"; W="/NatureUsers/glitch/mv_review"
exp={"MV-4122_S12":(19.1,5.6),"MV-4125_S18":(12.9,3.7),"MV-4128_S29":(20.3,6.0)}
with open(f"{W}/report/qc_table.tsv","w") as o:
    o.write("library\traw_pairs_M\tclean_pairs_M\tkept_pct\tclean_Gb\tdup_rate_pct\texpected_pairs_M\texpected_Gb\ttoo_short_pct_of_raw_reads\n")
    for f in sorted(glob.glob(f"{K}/qc/fastp_pilot/*.json")+glob.glob(f"{K}/qc/fastp_companions/*.json")):
        j=json.load(open(f)); lib=f.split("/")[-1][:-5]; b=j["summary"]["before_filtering"]; a=j["summary"]["after_filtering"]
        e=exp.get(lib,("",""))
        o.write(f"{lib}\t{b['total_reads']/2e6:.2f}\t{a['total_reads']/2e6:.2f}\t{100*a['total_reads']/b['total_reads']:.1f}\t{a['total_bases']/1e9:.2f}\t{100*j['duplication']['rate']:.1f}\t{e[0]}\t{e[1]}\t{100*j['filtering_result']['too_short_reads']/b['total_reads']:.1f}\n")
SC=["Ribosomal_L2_C","Ribosomal_L3","Ribosomal_L4","Ribosomal_L5","Ribosomal_L6","Ribosomal_L14","Ribosomal_L16","Ribosomal_L18p","Ribosomal_L22","Ribosomal_S3_C","Ribosomal_S8","Ribosomal_S10","Ribosomal_S17","Ribosomal_S19"]
n16={}
for l in open(f"{W}/contam/MV_16S.fna"):
    if l[0]==">":
        s=l[1:].split("|")[0]; L=int(re.search(r"len(\d+)",l).group(1)); n16.setdefault(s,[0,0]); n16[s][0]+=1; n16[s][1]+= L>=1200
with open(f"{W}/report/purity_table.tsv","w") as o:
    o.write("assembly\ttotal_Mb\tMb_in_contigs_ge5kb\tviral_VS2_Mb\tviral_VS2_pct_of_assembly\t16S_loci\t16S_ge1200bp\tribosomal_protein_hits_total\tsingle_copy_RP_median_copies\tRP_per_Mb\n")
    vs2={"MV-4122_S12":2.022702,"MV-4125_S18":2.071048,"MV-4128_S29":1.367457}
    for s in ["MV-4122_S12","MV-4125_S18","MV-4128_S29"]:
        tot=0;ge5=0
        for l in open(f"{K}/assembly/no_NC/{s}/contigs.fasta"):
            if l[0]==">":
                L=int(l.split("_")[3]); tot+=L; ge5+=L if L>=5000 else 0
        hits={}
        for l in open(f"{K}/qc/ribosomal_no_NC/{s}/ribosomal.tblout"):
            if l[0]=="#": continue
            f=l.split()
            if float(f[4])<1e-10: hits.setdefault(f[2],set()).add(f[0])
        allh=sum(len(v) for k,v in hits.items() if k.startswith("Ribosom") or k.startswith("RL"))
        med=statistics.median([len(hits.get(k,())) for k in SC])
        o.write(f"{s}\t{tot/1e6:.1f}\t{ge5/1e6:.1f}\t{vs2[s]:.2f}\t{100*vs2[s]/(tot/1e6):.1f}\t{n16[s][0]}\t{n16[s][1]}\t{allh}\t{med}\t{allh/(tot/1e6):.2f}\n")
