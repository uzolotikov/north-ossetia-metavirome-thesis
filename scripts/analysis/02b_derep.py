import sys,collections
fa,bl,out,mapf=sys.argv[1:5]
seqs={};h=None;buf=[]
for l in open(fa):
    if l[0]==">":
        if h: seqs[h]="".join(buf)
        h=l[1:].strip();buf=[]
    else: buf.append(l.strip())
seqs[h]="".join(buf)
iv=collections.defaultdict(list)
for l in open(bl):
    q,s,pid,L,qs,qe,ss,se,ql,sl=l.split("\t")
    if q==s or float(pid)<99: continue
    iv[(q,s)].append((int(qs),int(qe)))
def covered(m):
    m=sorted(m);c=0;s,e=m[0]
    for a,b in m[1:]:
        if a>e: c+=e-s+1;s,e=a,b
        else: e=max(e,b)
    return c+e-s+1
cover={k:covered(v) for k,v in iv.items()}
order=sorted(seqs,key=lambda x:(-len(seqs[x]),x))
kept=[];rep={};keptset=set()
partners=collections.defaultdict(list)
for (q,s),c in cover.items(): partners[q].append((s,c))
for x in order:
    r=None
    for s,c in partners.get(x,[]):
        if s in keptset and len(seqs[s])>=len(seqs[x]) and c/len(seqs[x])>=0.90:
            r=s;break
    if r is None: kept.append(x);keptset.add(x);rep[x]=x
    else: rep[x]=r
with open(out,"w") as o:
    for x in kept: o.write(f">{x}\n{seqs[x]}\n")
with open(mapf,"w") as o:
    o.write("contig\tlength\trepresentative\n")
    for x in order: o.write(f"{x}\t{len(seqs[x])}\t{rep[x]}\n")
print("input",len(seqs),"kept",len(kept),"Mb_kept",round(sum(len(seqs[x]) for x in kept)/1e6,1))
