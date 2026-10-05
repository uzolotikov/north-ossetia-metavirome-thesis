import sys,re
K="/NatureUsers/nsinyavskiy/metaviroms_kgt"
def readfa(p):
    seqs={};h=None;buf=[]
    for l in open(p):
        if l[0]==">":
            if h: seqs[h]="".join(buf)
            h=l[1:].split()[0];buf=[]
        else: buf.append(l.strip())
    if h: seqs[h]="".join(buf)
    return seqs
comp=str.maketrans("ACGTacgtN","TGCAtgcaN")
out=open(sys.argv[1],"w")
for s in ["MV-4122_S12","MV-4125_S18","MV-4128_S29"]:
    loci=[]
    for m in ["bac","arc"]:
        for l in open(f"{K}/qc/barrnap_no_NC/{s}.{m}.gff"):
            if l.startswith("#") or "Name=16S" not in l: continue
            f=l.split("\t"); loci.append((f[0],int(f[3]),int(f[4]),f[6],m))
    loci.sort(key=lambda x:(x[0],x[1]))
    keep=[]
    for x in loci:
        if keep and keep[-1][0]==x[0] and x[1]<=keep[-1][2]:
            if x[2]-x[1]>keep[-1][2]-keep[-1][1]: keep[-1]=x
        else: keep.append(x)
    seqs=readfa(f"{K}/assembly/no_NC/{s}/contigs.fasta")
    for c,a,b,st,m in keep:
        q=seqs[c][a-1:b]
        if st=="-": q=q.translate(comp)[::-1]
        cl=int(re.search(r"length_(\d+)",c).group(1)); cv=float(re.search(r"cov_([\d.]+)",c).group(1))
        out.write(f">{s}|{c}|{a}-{b}|len{len(q)}|ctg{cl}|cov{cv:.1f}\n{q}\n")
    print(s,"unique 16S loci:",len(keep),">=500bp:",sum(1 for k in keep if k[2]-k[1]>=499))
