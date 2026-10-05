"""Кластеризация вирусных последовательностей в vOTU (95 % идентичности, 85 % покрытия более короткой).

Использование: cluster_votu.py <strict_ge5kb.fna> <all_vs_all_ge5kb.tsv> <out_dir>
Таблица BLASTN «все против всех» — 14 колонок:
  qseqid sseqid pident length qlen slen qstart qend sstart send evalue bitscore nident btop
Описание метода — cluster_votu_README.md.
"""
import sys,csv,re,collections,pathlib,json,hashlib
P=pathlib.Path(sys.argv[3]);P.mkdir(parents=True,exist_ok=True)
source=pathlib.Path(sys.argv[1])
seqs={}
for block in source.read_text().split('>')[1:]:
 h,*ls=block.splitlines();seqs[h.split()[0]]=''.join(ls)
assert len(seqs)==119
rs=list(csv.reader(open(sys.argv[2]),delimiter='\t'));pairs=collections.defaultdict(list)
for r in rs:
 assert len(r)==14 and r[0] in seqs and r[1] in seqs
 assert int(r[4])==len(seqs[r[0]]) and int(r[5])==len(seqs[r[1]])
 if r[0]!=r[1]:pairs[(r[0],r[1])].append(r)
assert {r[0] for r in rs}==set(seqs)
# Parse all BTOP traces and remove reused bases on either side, retaining higher-bit-score alignments first.
def metrics(rows):
 usedq=set();useds=set();n=ident=0
 for r in sorted(rows,key=lambda r:(-float(r[11]),-int(r[3]),int(r[6]),int(r[8]))):
  q,s=int(r[6]),int(r[8]);dq=1 if int(r[7])>=q else -1;ds=1 if int(r[9])>=s else -1
  tokens=re.findall(r'\d+|[^\d]{2}',r[13]);assert ''.join(tokens)==r[13]
  cols=[]
  for token in tokens:
   if token.isdigit():
    for _ in range(int(token)):cols.append((q,s,1));q+=dq;s+=ds
   else:
    a,b=token;cols.append((q if a!='-' else None,s if b!='-' else None,0));q+=dq*(a!='-');s+=ds*(b!='-')
  assert len(cols)==int(r[3]) and sum(c[2] for c in cols)==int(r[12])
  assert q==int(r[7])+dq and s==int(r[9])+ds
  for a,b,match in cols:
   if (a is not None and a in usedq) or (b is not None and b in useds):continue
   if a is not None:usedq.add(a)
   if b is not None:useds.add(b)
   ident+=match;n+=1
 return n,ident,len(usedq),len(useds)
edge={};detail=[]
for (q,s),rows in pairs.items():
 choices=[]
 for strand in (1,-1):
  selected=[r for r in rows if (int(r[7])-int(r[6]))*(int(r[9])-int(r[8]))*strand>=0]
  if not selected:continue
  n,i,cq,cs=metrics(selected);af=(cq if len(seqs[q])<=len(seqs[s]) else cs)/min(len(seqs[q]),len(seqs[s]));ani=i/n if n else 0
  choices.append((ani>=.95 and af>=.85,af,ani,n,strand))
 ok,af,ani,n,strand=max(choices)
 detail.append(dict(query=q,subject=s,identity_pct=round(ani*100,5),shorter_coverage_pct=round(af*100,5),retained_alignment_columns=n,strand=strand,passes=ok))
 if ok:
  key=tuple(sorted((q,s)))
  if key not in edge or (af,ani)>edge[key][:2]:edge[key]=(af,ani,q,s)
# Longest-first representatives, no transitive merging.
clusters=[];assigned={}
for q in sorted(seqs,key=lambda x:(-len(seqs[x]),x)):
 if q in assigned:continue
 cid=f'vOTU_{len(clusters)+1:03d}';members=[q];assigned[q]=cid
 for s in sorted(seqs,key=lambda x:(-len(seqs[x]),x)):
  if s not in assigned and tuple(sorted((q,s))) in edge:members.append(s);assigned[s]=cid
 clusters.append((cid,q,members))
assert len(assigned)==119
membership=[];summary=[]
for cid,rep,members in clusters:
 counts=collections.Counter(x.split('__')[0] for x in members)
 missing=[(a,b) for i,a in enumerate(members) for b in members[i+1:] if tuple(sorted((a,b))) not in edge]
 summary.append(dict(votu=cid,representative=rep,representative_length=len(seqs[rep]),members=len(members),samples=len(counts),MV_4122_S12=counts['MV-4122_S12'],MV_4125_S18=counts['MV-4125_S18'],MV_4128_S29=counts['MV-4128_S29'],nonpassing_member_pairs=len(missing)))
 for m in members:
  af,ani=(1,1) if m==rep else edge[tuple(sorted((rep,m)))][:2]
  membership.append(dict(votu=cid,sequence_id=m,sample=m.split('__')[0],length=len(seqs[m]),representative=rep,is_representative=m==rep,identity_to_rep_pct=round(ani*100,5),shorter_coverage_to_rep_pct=round(af*100,5)))
def write(name,rows):
 with (P/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
write('clusters.tsv',summary);write('membership.tsv',membership);write('pairwise_metrics.tsv',detail)
with (P/'representatives_ge5kb.fna').open('w') as f:
 for cid,rep,m in clusters:f.write(f'>{cid} original_id={rep}\n'+ '\n'.join(seqs[rep][i:i+80] for i in range(0,len(seqs[rep]),80))+'\n')
print('clusters',len(clusters),'edges',len(edge),'patterns',collections.Counter(tuple(x for x in ['MV_4122_S12','MV_4125_S18','MV_4128_S29'] if r[x]) for r in summary),'rep10k',sum(len(seqs[r])>=10000 for _,r,_ in clusters))
print('noncliques',[r for r in summary if r['nonpassing_member_pairs']]);print('multimember',[r for r in summary if r['members']>1])
print('inputseq>=10k',sum(len(s)>=10000 for s in seqs.values()))
print('crosscluster_edges',sum(assigned[a]!=assigned[b] for a,b in edge))
