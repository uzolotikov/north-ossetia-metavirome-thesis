"""Отбор вирусных кандидатов по результатам VirSorter2 и CheckV.

Использование: filter_candidates.py <checkv_dir> <virsorter_dir> <out_dir>
  checkv_dir/<образец>/    quality_summary.tsv, contamination.tsv, viruses.fna, proviruses.fna (checkv end_to_end)
  virsorter_dir/<образец>/ final-viral-score.tsv (virsorter run)
Правила отбора — filter_candidates_README.md.
"""
import sys,csv,io,pathlib,re,collections,json,hashlib
CHECKV=pathlib.Path(sys.argv[1]); VS=pathlib.Path(sys.argv[2]); O=pathlib.Path(sys.argv[3]); O.mkdir(parents=True,exist_ok=True)
def table(txt): return list(csv.DictReader(io.StringIO(txt),delimiter='\t'))
def fasta(txt):
 out=[]
 for block in txt.split('>')[1:]:
  h,*seq=block.splitlines();out.append((h,''.join(seq).upper()))
 return out
def write_table(name,rows):
 with (O/name).open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t');w.writeheader();w.writerows(rows)
sets=collections.defaultdict(list); parents=[]; regions=[]; summary=[]
flagged={'MV-4122_S12':{'NODE_2365','NODE_3558'},'MV-4125_S18':{'NODE_1759'},'MV-4128_S29':{'NODE_726'}}
if True:
 for sample in ['MV-4122_S12','MV-4125_S18','MV-4128_S29']:
  read=lambda n:(CHECKV/sample/n).read_text()
  qs=table(read('quality_summary.tsv')); cs={r['contig_id']:r for r in table(read('contamination.tsv'))}
  vs={r['seqname']:r for r in table((VS/sample/'final-viral-score.tsv').read_text())}
  assert {q['contig_id'] for q in qs}==set(vs)
  seqs=collections.defaultdict(list)
  for h,seq in fasta(read('viruses.fna')):
   name=h.split()[0];assert name in vs
   seqs[name].append((name,seq,1,len(seq),False))
  for h,seq in fasta(read('proviruses.fna')):
   name,coords=h.split();parent=name.rsplit('_',1)[0];assert parent in vs
   m=re.fullmatch(r'(\d+)-(\d+)/(\d+)',coords); assert m,h
   a,b,L=map(int,m.groups());assert len(seq)==b-a+1
   seqs[parent].append((name,seq,a,b,True))
  assert set(seqs)==set(vs)
  for q in qs:
   cid=q['contig_id'];v=vs[cid];c=cs[cid];score=float(v['max_score']);hall=float(v['hallmark']);vg=int(q['viral_genes']);hg=int(q['host_genes']);hq=q['checkv_quality'] in ('Complete','High-quality')
   node='_'.join(cid.split('_')[:2]); flag=node in flagged[sample]
   support=(score>=.9 or (score>=.7 and hall>=1))
   if flag: status,reason='exclude_strict','ribosomal_hit_in_candidate;no_CheckV_viral_genes;host_genes_present;no_VS_hallmark'
   elif vg==0 and hg>0:status,reason='exclude_strict','no_CheckV_viral_genes;host_genes_present'
   elif vg==0:status,reason='review','no_CheckV_viral_genes;insufficient_support'
   elif q['warnings']:status,reason='review','CheckV_warning_requires_review'
   elif support or hq:status,reason='strict_supported',('VS_support_and_CheckV_viral_gene' if support else 'CheckV_Complete_or_High_quality_exception')
   else:status,reason='review','VS_support_below_working_threshold'
   rr=seqs[cid]
   if q['provirus']=='Yes':
    expected=[tuple(map(int,co.split('-'))) for typ,co in zip(c['region_types'].split(','),c['region_coords_bp'].split(',')) if typ=='viral']
    assert sorted((r[2],r[3]) for r in rr)==sorted(expected)
    assert sum(len(r[1]) for r in rr)==int(q['proviral_length'])
   else:assert len(rr)==1 and len(rr[0][1])==int(q['contig_length'])
   parentrow=dict(sample=sample,candidate_id=cid,decision=status,reason=reason,vs_score=v['max_score'],vs_hallmark=v['hallmark'],checkv_quality_before_trimming=q['checkv_quality'],completeness_before_trimming=q['completeness'],completeness_method=q['completeness_method'],viral_genes_before_trimming=vg,host_genes_before_trimming=hg,checkv_warnings=q['warnings'],input_length=int(q['contig_length']),retained_length=sum(len(r[1]) for r in rr),output_regions=len(rr),checkv_trimmed=q['provirus'])
   parents.append(parentrow)
   for name,seq,a,b,trimmed in rr:
    uid=f'{sample}__{name}';assert set(seq)<=set('ACGTN'),uid
    region=dict(parentrow,sequence_id=uid,checkv_sequence_id=name,start_in_CheckV_input=a,end_in_CheckV_input=b,sequence_length=len(seq),strict_ge5kb=status=='strict_supported' and len(seq)>=5000,strict_ge10kb=status=='strict_supported' and len(seq)>=10000,short_HQ=status=='strict_supported' and hq and len(seq)<5000)
    regions.append(region);sets['all_after_checkv'].append((uid,seq))
    sets[status].append((uid,seq))
    if region['strict_ge5kb']:sets['strict_ge5kb'].append((uid,seq))
    if region['strict_ge10kb']:sets['strict_ge10kb'].append((uid,seq))
    if region['short_HQ']:sets['short_HQ'].append((uid,seq))
  ps=[r for r in parents if r['sample']==sample];rs=[r for r in regions if r['sample']==sample]
  summary.append(dict(sample=sample,input_candidates=len(ps),trimmed_candidates=sum(r['checkv_trimmed']=='Yes' for r in ps),strict_supported=sum(r['decision']=='strict_supported' for r in ps),review=sum(r['decision']=='review' for r in ps),exclude_strict=sum(r['decision']=='exclude_strict' for r in ps),strict_ge5kb=sum(r['strict_ge5kb'] for r in rs),strict_ge10kb=sum(r['strict_ge10kb'] for r in rs),short_HQ=sum(r['short_HQ'] for r in rs),removed_host_bp=sum(r['input_length']-r['retained_length'] for r in ps)))
for key in ['all_after_checkv','strict_supported','review','exclude_strict','strict_ge5kb','strict_ge10kb','short_HQ']:
 records=sets[key];assert len({r[0] for r in records})==len(records)
 with (O/(key+'.fna')).open('w') as f:
  for name,seq in records:f.write('>'+name+'\n'+'\n'.join(seq[i:i+80] for i in range(0,len(seq),80))+'\n')
write_table('candidate_decisions.tsv',parents);write_table('sequence_manifest.tsv',regions);write_table('summary.tsv',summary)
print(json.dumps(summary,ensure_ascii=False,indent=2));print('TOTAL',len(parents),len(regions));print('strict quality',collections.Counter(r['checkv_quality_before_trimming'] for r in parents if r['decision']=='strict_supported'))
with (O/'input_sha256.txt').open('w') as f:
 for sample in ['MV-4122_S12','MV-4125_S18','MV-4128_S29']:
  for path in [CHECKV/sample/n for n in ('quality_summary.tsv','contamination.tsv','viruses.fna','proviruses.fna')]+[VS/sample/'final-viral-score.tsv']:
   f.write(hashlib.sha256(path.read_bytes()).hexdigest()+'  '+sample+'/'+path.name+'\n')
