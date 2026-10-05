import sys
def stats(p):
    L=[]
    for l in open(p):
        if l[0]==">": L.append(int(l.split("_length_")[1].split("_")[0]))
    L.sort(reverse=True); T=sum(L); c=0; n50=0
    for x in L:
        c+=x
        if c>=T/2: n50=x; break
    return len(L),T/1e6,max(L),n50,sum(x>=5000 for x in L),sum(x>=10000 for x in L),sum(x for x in L if x>=5000)/1e6
print("sample\tassembler\tcontigs\ttotal_Mb\tlargest\tN50\tn_ge5kb\tn_ge10kb\tMb_ge5kb")
for s in sys.argv[1:]:
    for a,p in [("metaSPAdes",f"/NatureUsers/nsinyavskiy/metaviroms_kgt/assembly/no_NC/{s}/contigs.fasta"),("SPAdes",f"/NatureUsers/glitch/mv_review/assembly_spades/{s}/contigs.fasta")]:
        try: r=stats(p); print(f"{s}\t{a}\t{r[0]}\t{r[1]:.1f}\t{r[2]}\t{r[3]}\t{r[4]}\t{r[5]}\t{r[6]:.1f}")
        except FileNotFoundError: print(f"{s}\t{a}\tNA")
