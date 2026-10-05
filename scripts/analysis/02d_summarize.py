import json,glob,os,collections,sys
K="/NatureUsers/nsinyavskiy/metaviroms_kgt"; W="/NatureUsers/glitch/mv_review"; S=f"{W}/step6"
LIBS=["MV-4122_S12","MV-4125_S18","MV-4128_S29","4122-0-22_S11","4122-0-44_S10","4125-0-22_S17","4125-0-44_S16","4128-0-22_S28","4128-0-44_S27","4128_S30","NC1_S64"]
CELL={"4122":["4122-0-22_S11","4122-0-44_S10"],"4125":["4125-0-22_S17","4125-0-44_S16"],"4128":["4128-0-22_S28","4128-0-44_S27","4128_S30"]}
ALLCELL=sum(CELL.values(),[])
tot={}
for f in glob.glob(f"{K}/qc/fastp_*/*.json"):
    lib=os.path.basename(f)[:-5]; tot[lib]=json.load(open(f))["filtering_result"]["passed_filter_reads"]
cov={}
for lib in LIBS:
    p=f"{S}/cov/{lib}.cov_mq0.tsv"
    if not os.path.exists(p): continue
    for l in open(p):
        if l[0]=="#": continue
        f=l.rstrip().split("\t"); c=f[0]
        cov.setdefault(c,{"len":int(f[2])})[lib]=(int(f[3]),float(f[5]),float(f[6]))  # numreads, breadth%, meandepth
libs=[l for l in LIBS if any(l in v for v in cov.values())]
# viral classes from geNomad if available
vir={}
gs=glob.glob(f"{W}/genomad/ref_derep_summary/ref_derep_virus_summary.tsv")
if gs:
    for l in open(gs[0]):
        if l.startswith("seq_name"): continue
        f=l.split("\t"); name=f[0].split("|provirus")[0]
        vir[name]="provirus" if "|provirus" in f[0] else "virus"
def origin(c):
    if c.startswith("MV-"): return "MV_assembly"
    if c.startswith("cell-"): return "Tib_cellular_assembly"
    return "Biragzang_final"
def cls(c):
    d=cov[c]
    nd=lambda l: d.get(l,(0,0,0))[2]/tot[l]*1e6
    mv=max(nd(l) for l in libs if l.startswith("MV-")); ce=max(nd(l) for l in libs if l in ALLCELL)
    v=vir.get(c)
    if v=="virus": base="viral"
    elif v=="provirus": base="cellular_with_provirus"
    else: base="cellular_nonviral"
    if mv==0 and ce==0: e="unmapped"
    elif ce==0 or mv/ce>=100: e="MVspecific(>=100x)"
    elif mv/ce>=10: e="MVenriched(10-100x)"
    elif mv/ce>=0.1: e="shared(0.1-10x)"
    else: e="cellular_enriched(<0.1x)"
    return base+"|"+e
out=open(f"{W}/report/step6_read_fractions.tsv","w")
groups_o=["MV_assembly","Tib_cellular_assembly","Biragzang_final"]
C=collections.Counter(cls(c) for c in cov)
groups_c=sorted(C)
out.write("library\tclean_reads\tmapped_pct\t"+"\t".join("orig:"+g for g in groups_o)+"\t"+"\t".join("class:"+g for g in groups_c)+"\n")
for lib in libs:
    T=tot[lib]; mo=collections.Counter(); mc=collections.Counter()
    for c,d in cov.items():
        n=d.get(lib,(0,0,0))[0]; mo[origin(c)]+=n; mc[cls(c)]+=n
    m=sum(mo.values())
    out.write(f"{lib}\t{T}\t{100*m/T:.3f}\t"+"\t".join(f"{100*mo[g]/T:.4f}" for g in groups_o)+"\t"+"\t".join(f"{100*mc[g]/T:.4f}" for g in groups_c)+"\n")
out.close()
with open(f"{W}/report/step6_contig_classes.tsv","w") as o:
    o.write("contig\tlength\torigin\tclass\t"+"\t".join(f"depth:{l}" for l in libs)+"\n")
    for c,d in cov.items():
        o.write(f"{c}\t{d['len']}\t{origin(c)}\t{cls(c)}\t"+"\t".join(f"{d.get(l,(0,0,0))[2]:.3f}" for l in libs)+"\n")
print("contigs:",len(cov),"classes:",dict(C),"geNomad loaded:",bool(vir))
