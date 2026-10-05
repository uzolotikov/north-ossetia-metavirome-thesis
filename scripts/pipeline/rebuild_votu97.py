#!/usr/bin/env python3
"""Сводная таблица покрытия каталога vOTU по 11 библиотекам.

Входы — таблицы mapping/<группа>/coverage/<библиотека>.{coverage,idxstats}.tsv (coverage_tables.py).
Для каждой vOTU: число отфильтрованных прочтений, RPM, RPKM, средняя глубина, ширина покрытия
и детекция (ширина >=75 % и средняя глубина >=1x). Размеры библиотек — число прочтений после
тримминга (оба мейта пары), по отчётам fastp.

Использование: rebuild_votu97.py <корень_проекта> <папка_вывода>
"""
import csv, pathlib, json, hashlib, math, sys

# Число прочтений после тримминга (оба мейта пары) по отчётам fastp.
LIBRARIES = 'sample\tgroup\tclean_reads\tcoverage_dir\nMV-4122_S12\tMV\t38071806\tmapping/no_NC/coverage\nMV-4125_S18\tMV\t25680390\tmapping/no_NC/coverage\nMV-4128_S29\tMV\t40598794\tmapping/no_NC/coverage\nNC1_S64\tnegative_control\t103764\tmapping/no_NC/coverage\n4122-0-22_S11\tcellular\t1193866\tmapping/companions/coverage\n4122-0-44_S10\tcellular\t30115448\tmapping/companions/coverage\n4125-0-22_S17\tcellular\t4130804\tmapping/companions/coverage\n4125-0-44_S16\tcellular\t19880220\tmapping/companions/coverage\n4128-0-22_S28\tcellular\t1001786\tmapping/companions/coverage\n4128-0-44_S27\tcellular\t3599252\tmapping/companions/coverage\n4128_S30\twhole_water\t6453914\tmapping/companions/coverage\n'

def read(p):
    with p.open(newline='') as f:
        return list(csv.DictReader(f, delimiter='\t'))

def write(p, rows):
    with p.open('w', newline='') as f:
        w=csv.DictWriter(f, fieldnames=list(rows[0]), delimiter='\t')
        w.writeheader(); w.writerows(rows)

def main():
    from types import SimpleNamespace
    a = SimpleNamespace(
        project=pathlib.Path(sys.argv[1]),
        outdir=pathlib.Path(sys.argv[2]),
    )
    if a.outdir.exists():
        raise SystemExit('Output directory already exists; choose a new directory.')
    rows=[]; summaries=[]; provenance={}; catalog=None; samples=set()
    def track(p):
        provenance[str(p)]=hashlib.sha256(p.read_bytes()).hexdigest()
    for m in csv.DictReader(LIBRARIES.splitlines(), delimiter="\t"):
        s=m['sample']; group=m['group']; n=int(m['clean_reads'])
        assert s not in samples and n>0, s
        samples.add(s)
        base=a.project/m['coverage_dir']
        cp=base/(s+'.coverage.tsv'); ip=base/(s+'.idxstats.tsv')
        track(cp); track(ip)
        idx={}
        for line in ip.read_text().splitlines():
            ident,length,mapped,unmapped=line.split('\t')
            if ident=='*': continue
            assert ident not in idx
            idx[ident]=(int(length),int(mapped))
        block=[]; observed={}
        for r in read(cp):
            v=r['votu']; length=int(r['length']); count=int(r['mapped_reads'])
            depth=float(r['mean_depth']); breadth=float(r['breadth_1x_pct'])
            assert r['sample']==s and v not in observed
            assert length>0 and count>=0 and math.isfinite(depth) and depth>=0
            assert math.isfinite(breadth) and 0<=breadth<=100
            assert idx[v]==(length,count), (s,v,'coverage/idxstats mismatch')
            observed[v]=length
            block.append(dict(sample=s,group=group,votu=v,length_bp=length,
                clean_reads=n,mapped_reads_filtered=count,
                RPM=round(count*1e6/n,6),RPKM=round(count*1e9/(length*n),6),
                mean_depth=depth,breadth_1x_pct=breadth,
                detected_75pct_1x=int(breadth>=75 and depth>=1)))
        assert len(observed)==97 and set(observed)==set(idx)
        if catalog is None: catalog=observed
        assert observed==catalog, 'Different reference catalogs'
        block.sort(key=lambda r:r['votu']); rows.extend(block)
        total=sum(r['mapped_reads_filtered'] for r in block)
        summaries.append(dict(sample=s,group=group,clean_reads=n,
            mapped_reads_filtered=total,mapped_pct_filtered=round(total*100/n,6),
            detected_vOTUs=sum(r['detected_75pct_1x'] for r in block)))
    a.outdir.mkdir(parents=True)
    write(a.outdir/'votu97_all_libraries.tsv',rows)
    write(a.outdir/'library_summary.tsv',summaries)
    (a.outdir/'input_sha256.json').write_text(json.dumps(provenance,indent=2)+'\n')
    print(f'Written {len(rows)} rows for {len(samples)} libraries to {a.outdir}')
if __name__=='__main__': main()
