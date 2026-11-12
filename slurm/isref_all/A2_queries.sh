#!/bin/bash
#SBATCH --job-name=isall_A2
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=10:00:00
#SBATCH --output=%x_%j.out
# STAGE A2 -- query set, rebuilt after the single-linkage version was rejected.
#
# WHAT WAS WRONG.  A1 put sequences in one component when an alignment covered
# >=80% of the SHORTER one, then kept only the longest member.  A 1.3 kb IS
# sitting inside a 5 kb Tn3 covers 100% of the IS, so the two merged and the IS
# stopped being searched in its own right.  Shared fragments were silently
# fusing distinct insertion structures.
#
# WHAT THIS DOES INSTEAD.
#   * Every bag keeps its ORIGINAL representative (the protein-cluster
#     representative's insertion, the same one 232 used).  Nothing replaces it.
#   * A member is promoted to an extra representative unless some already-kept
#     representative of that bag covers it under a TWO-SIDED test: >=90%
#     identity AND >=90% of the member AND >=90% of the representative, using
#     merged HSP blocks on both sides, not one HSP and not the shorter side
#     alone.  Containment inside a larger element therefore does not absorb it.
#   * Nothing is deleted across bags.  Cross-bag DNA equivalence is RECORDED
#     and applied later at the locus level, where collapsing costs no
#     sensitivity, instead of by dropping queries up front.
#   * Both coverages are written out for every decision.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
mkdir -p "$W/tmp"

# A1's selfaln.tsv lacked alignment coordinates, so merged two-sided coverage
# cannot be computed from it.  It is not compatible with this stage and is not
# reused; this recomputes with coordinates.
if [ ! -s "$W/selfaln2.tsv" ]; then
  [ -s "$W/tmp/alldb.nsq" ] || "$BL/makeblastdb" -in "$C/dedup_insertions.fna" \
      -dbtype nucl -out "$W/tmp/alldb" -title all >/dev/null
  "$BL/blastn" -task megablast -query "$C/dedup_insertions.fna" -db "$W/tmp/alldb" \
     -evalue 1e-50 -perc_identity 90 -max_target_seqs 100000 -num_threads 64 \
     -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
     -out "$W/selfaln2.tsv.part" 2>"$W/tmp/selfaln2.log"
  mv "$W/selfaln2.tsv.part" "$W/selfaln2.tsv"
fi
echo "self-alignment rows: $(wc -l < "$W/selfaln2.tsv")"

"$PY" - <<'PY'
import csv, collections, json
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all"
C="/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"
ID, COV = 90.0, 0.90

def mlen(iv):
    iv=sorted(iv); t=0; cs,ce=iv[0]
    for s,e in iv[1:]:
        if s<=ce: ce=max(ce,e)
        else: t+=ce-cs+1; cs,ce=s,e
    return t+ce-cs+1

# --- bag membership and each bag's ORIGINAL representative
bagkeys=collections.defaultdict(set)
for r in csv.DictReader(open(C+"/bags/strict_occurrence_bags.tsv"),delimiter="\t"):
    bagkeys[r["cds_cluster_id"]].add(r["seq_key"])
orig={}
for r in csv.DictReader(open(C+"/bags/insertion_cds_strict.tsv"),delimiter="\t"):
    b=r.get("cds_cluster_id")
    if b and b not in orig:
        orig[b]=r.get("cluster_rep") or r.get("seq_key")
for b in bagkeys:
    if b not in orig or orig[b] not in bagkeys[b]:
        orig[b]=sorted(bagkeys[b])[0]
print("bags %d, original representatives %d"%(len(bagkeys),len(orig)))

seqs={}; n=None; buf=[]
def fl():
    if n: seqs[n.split()[0].replace("SEQ|","")]="".join(buf)
for line in open(C+"/dedup_insertions.fna"):
    if line[0]==">": fl(); n=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
fl()
print("insertion sequences:",len(seqs))

# --- merged two-sided coverage for every aligned pair
qb=collections.defaultdict(list); sb=collections.defaultdict(list); pid={}
for l in open(W+"/selfaln2.tsv"):
    q,s,p,aln,qs,qe,ss,se,ql,sl=l.rstrip("\n").split("\t")
    q=q.replace("SEQ|",""); s=s.replace("SEQ|","")
    if q==s or float(p)<ID: continue
    qs,qe,ss,se=int(qs),int(qe),int(ss),int(se)
    if se<ss: ss,se=se,ss
    qb[(q,s)].append((qs,qe)); sb[(q,s)].append((ss,se))
    pid[(q,s)]=max(pid.get((q,s),0),float(p))
cov={}
for k in qb:
    q,s=k
    cov[k]=(mlen(qb[k])/len(seqs[q]), mlen(sb[k])/len(seqs[s]))
print("aligned pairs with coverage computed:",len(cov))

def covered(member, rep):
    """rep covers member only if BOTH sequences are >=COV covered."""
    c=cov.get((member,rep)) or cov.get((rep,member))
    if not c: return None
    a,b = c if (member,rep) in cov else (c[1],c[0])
    return (a,b) if (a>=COV and b>=COV) else None

# --- per bag: keep the original, promote genuinely distinct variants
reps=collections.defaultdict(list); decisions=[]
for b,keys in sorted(bagkeys.items()):
    o=orig[b]; reps[b]=[o]
    for m in sorted(keys-{o}, key=lambda k:(-len(seqs.get(k,"")),k)):
        if m not in seqs: continue
        hit=None
        for r in reps[b]:
            c=covered(m,r)
            if c: hit=(r,c); break
        if hit:
            decisions.append((b,m,"absorbed",hit[0],round(hit[1][0],3),round(hit[1][1],3)))
        else:
            reps[b].append(m); decisions.append((b,m,"kept_as_variant","",0.0,0.0))
with open(W+"/query_decisions.tsv","w") as o:
    o.write("bag_id\tseq_key\tdecision\tabsorbed_by\tcov_of_member\tcov_of_representative\n")
    for d in decisions: o.write("\t".join(str(x) for x in d)+"\n")

# --- cross-bag DNA equivalence, RECORDED not applied
equiv=[]
allreps={k for v in reps.values() for k in v}
for (q,s),(a,bcov) in cov.items():
    if q in allreps and s in allreps and q<s and a>=COV and bcov>=COV and pid[(q,s)]>=ID:
        equiv.append((q,s,round(a,3),round(bcov,3),pid[(q,s)]))
with open(W+"/query_equivalence.tsv","w") as o:
    o.write("query_a\tquery_b\tcov_a\tcov_b\tpident\n")
    for e in sorted(equiv): o.write("\t".join(str(x) for x in e)+"\n")

# --- write queries; one FASTA record per (bag, representative)
rows=[]
for b,ks in sorted(reps.items()):
    for i,k in enumerate(ks):
        rows.append((("Q_%s"%k), b, k, i==0, len(seqs[k])))
seen=set()
with open(W+"/queries.fna","w") as fa:
    for qid,b,k,isorig,L in rows:
        if k in seen: continue
        seen.add(k); fa.write(">%s\n%s\n"%(qid,seqs[k]))
with open(W+"/queries.tsv","w") as t:
    t.write("query_id\tbag_id\tseq_key\tis_bag_original_rep\tlen_nt\n")
    for qid,b,k,isorig,L in rows:
        t.write("%s\t%s\t%s\t%d\t%d\n"%(qid,b,k,int(isorig),L))
nv=collections.Counter(b for b,ks in reps.items() for _ in ks)
print(json.dumps({
 "bags":len(reps),
 "query_records_bag_x_rep":len(rows),
 "distinct_query_sequences":len(seen),
 "total_query_bp":sum(len(seqs[k]) for k in seen),
 "bags_with_one_representative":sum(1 for b,c in nv.items() if c==1),
 "bags_with_multiple_representatives":sum(1 for b,c in nv.items() if c>1),
 "max_representatives_for_one_bag":max(nv.values()),
 "members_absorbed_two_sided":sum(1 for d in decisions if d[2]=="absorbed"),
 "members_kept_as_distinct_variants":sum(1 for d in decisions if d[2]=="kept_as_variant"),
 "cross_bag_equivalent_query_pairs_recorded":len(equiv)},indent=2))
PY
