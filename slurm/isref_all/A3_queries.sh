#!/bin/bash
#SBATCH --job-name=isall_A3
#SBATCH --account=pc_rubinlab
#SBATCH --partition=cf1
#SBATCH --qos=cf_normal
#SBATCH --nodes=1
#SBATCH --cpus-per-task=64
#SBATCH --mem=120G
#SBATCH --time=12:00:00
#SBATCH --output=%x_%j.out
# STAGE A3 -- query set.  Same SELECTION RULE as A2; cheaper COMPARISON SET.
#
# A2 was cancelled for cost, not for semantics.  It blasted all 63,565 dedup
# insertions against each other (4e9 pairs): 78 GB resident and climbing ~1
# GB/min against a 120 G cap, 750% of 6400% CPU, no output after 62 min.  It
# would have been OOM-killed.  The selection rule only ever consults pairs
# WITHIN one bag -- sum of n^2 over the 6,144 bags is 1.0e8, forty times less
# -- so the cross-bag quadratic was pure waste.
#
# THE RULE IS UNCHANGED FROM A2, which is the point:
#   * Every bag keeps its ORIGINAL representative (cluster_rep from
#     insertion_cds_strict.tsv, the same one 232 used).  Nothing replaces it.
#   * A member is promoted to an extra representative unless an already-kept
#     representative of that bag covers it TWO-SIDED: >=90% identity AND >=90%
#     of the member AND >=90% of the representative, from merged HSP blocks on
#     both sides.  A 1.3 kb IS inside a 5 kb Tn3 covers 100% of the IS and 26%
#     of the Tn3, so it is kept, not absorbed.
#   * Nothing is deleted across bags.  Cross-bag DNA equivalence is RECORDED
#     to query_equivalence.tsv and applied later at the LOCUS level, where
#     collapsing costs no sensitivity, instead of by dropping queries up front.
#
# Pass 2 computes that equivalence on the REPRESENTATIVES only, which is
# exactly the set it is used for, instead of on all members.
#
# RESUMABLE at per-bag grain: a bag is skipped only if its .done marker exists,
# which is written after blastn exits 0, so an interrupted bag is recomputed.
set -euo pipefail
W=/global/scratch/users/kh36969/fna_ins_discovery/isref_all
C=/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1
PY=/global/home/users/kh36969/.conda/envs/claude-env/bin/python3
BL=/global/home/users/kh36969/.conda/envs/fnains_annot/bin
mkdir -p "$W/bagfa" "$W/bagaln" "$W/tmp"
export PATH="$BL:$PATH"

# ---------- prep: one FASTA per bag ----------
if [ ! -s "$W/bagfa/.prepared" ]; then
"$PY" - <<'PY'
import csv, collections, os
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all"
C="/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"
bagkeys=collections.defaultdict(set)
for r in csv.DictReader(open(C+"/bags/strict_occurrence_bags.tsv"),delimiter="\t"):
    bagkeys[r["cds_cluster_id"]].add(r["seq_key"])
seqs={}; n=None; buf=[]
def fl():
    if n: seqs[n.split()[0].replace("SEQ|","")]="".join(buf)
for line in open(C+"/dedup_insertions.fna"):
    if line[0]==">": fl(); n=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
fl()
nb=0; miss=0
for b,ks in bagkeys.items():
    ks=[k for k in sorted(ks) if k in seqs]
    miss+=len(bagkeys[b])-len(ks)
    if not ks: continue
    with open(W+"/bagfa/%s.fna"%b,"w") as f:
        for k in ks: f.write(">%s\n%s\n"%(k,seqs[k]))
    nb+=1
open(W+"/bagfa/.prepared","w").write("%d\n"%nb)
print("bag fastas written: %d   members with no sequence: %d"%(nb,miss))
PY
fi

# ---------- pass 1: within-bag all-vs-all, largest bags first ----------
alnbag() {
  b=$1
  [ -s "$W/bagaln/$b.done" ] && return 0
  fa="$W/bagfa/$b.fna"
  nseq=$(grep -c '^>' "$fa")
  if [ "$nseq" -lt 2 ]; then : > "$W/bagaln/$b.tsv"; echo "$nseq" > "$W/bagaln/$b.done"; return 0; fi
  d=$(mktemp -d "${TMPDIR:-/tmp}/bagdb.XXXXXX")
  if makeblastdb -in "$fa" -dbtype nucl -out "$d/db" -title "$b" >/dev/null 2>&1 &&
     blastn -task megablast -query "$fa" -db "$d/db" -evalue 1e-50 -perc_identity 90 \
            -max_target_seqs 100000 -num_threads 1 \
            -outfmt '6 qseqid sseqid pident length qstart qend sstart send qlen slen' \
            -out "$W/bagaln/$b.tsv.part" 2>>"$W/tmp/A3_blast.err"; then
    mv "$W/bagaln/$b.tsv.part" "$W/bagaln/$b.tsv"
    echo "$nseq" > "$W/bagaln/$b.done"
  else
    rm -f "$W/bagaln/$b.tsv.part"
    echo "FAILED $b" >> "$W/tmp/A3_failed.txt"
  fi
  rm -rf "$d"
}
export -f alnbag; export W

# biggest bags first so the 7,014-member bag is not the straggler
for f in "$W/bagfa"/*.fna; do b=$(basename "$f" .fna); echo -e "$(stat -c%s "$f")\t$b"; done \
  | sort -k1,1nr | cut -f2 > "$W/tmp/bagorder.txt"
echo "bags to align: $(wc -l < "$W/tmp/bagorder.txt")"
xargs -a "$W/tmp/bagorder.txt" -P 60 -n 1 -I{} bash -c 'alnbag "$@"' _ {}

ndone=$(ls "$W/bagaln"/*.done 2>/dev/null | wc -l)
nwant=$(wc -l < "$W/tmp/bagorder.txt")
echo "bags aligned: $ndone / $nwant"
if [ -s "$W/tmp/A3_failed.txt" ]; then echo "FAILURES:"; cat "$W/tmp/A3_failed.txt"; fi
if [ "$ndone" -ne "$nwant" ]; then echo "incomplete -- refusing to select representatives"; exit 1; fi

cat "$W/bagaln"/*.tsv > "$W/selfaln_bag.tsv"
echo "within-bag alignment rows: $(wc -l < "$W/selfaln_bag.tsv")"

# ---------- pass 2: selection, then rep-vs-rep equivalence ----------
"$PY" - <<'PY'
import csv, collections, json, os, subprocess
W="/global/scratch/users/kh36969/fna_ins_discovery/isref_all"
C="/global/scratch/users/kh36969/fna_ins_discovery/cds_pfam_v1"
BL="/global/home/users/kh36969/.conda/envs/fnains_annot/bin"
ID, COV = 90.0, 0.90

def mlen(iv):
    iv=sorted(iv); t=0; cs,ce=iv[0]
    for s,e in iv[1:]:
        if s<=ce: ce=max(ce,e)
        else: t+=ce-cs+1; cs,ce=s,e
    return t+ce-cs+1

bagkeys=collections.defaultdict(set)
for r in csv.DictReader(open(C+"/bags/strict_occurrence_bags.tsv"),delimiter="\t"):
    bagkeys[r["cds_cluster_id"]].add(r["seq_key"])
orig={}
for r in csv.DictReader(open(C+"/bags/insertion_cds_strict.tsv"),delimiter="\t"):
    b=r.get("cds_cluster_id")
    if b and b not in orig:
        orig[b]=r.get("cluster_rep") or r.get("seq_key")
seqs={}; n=None; buf=[]
def fl():
    if n: seqs[n.split()[0].replace("SEQ|","")]="".join(buf)
for line in open(C+"/dedup_insertions.fna"):
    if line[0]==">": fl(); n=line[1:].strip(); buf=[]
    else: buf.append(line.strip())
fl()
# an original representative that is absent from its own bag or has no sequence
# is reported, not silently swapped
repaired=[]
for b in list(bagkeys):
    bagkeys[b]={k for k in bagkeys[b] if k in seqs}
    if not bagkeys[b]: del bagkeys[b]; continue
    if b not in orig or orig[b] not in bagkeys[b]:
        was=orig.get(b,"<none>"); orig[b]=sorted(bagkeys[b])[0]
        repaired.append((b,was,orig[b]))
print("bags %d   original representatives usable %d   repaired %d"
      %(len(bagkeys),len(bagkeys)-len(repaired),len(repaired)))

qb=collections.defaultdict(list); sb=collections.defaultdict(list)
for l in open(W+"/selfaln_bag.tsv"):
    q,s,p,aln,qs,qe,ss,se,ql,sl=l.rstrip("\n").split("\t")
    if q==s or float(p)<ID: continue
    qs,qe,ss,se=int(qs),int(qe),int(ss),int(se)
    if se<ss: ss,se=se,ss
    qb[(q,s)].append((qs,qe)); sb[(q,s)].append((ss,se))
cov={k:(mlen(qb[k])/len(seqs[k[0]]), mlen(sb[k])/len(seqs[k[1]])) for k in qb}
print("within-bag aligned pairs with two-sided coverage:",len(cov))

def covered(member, rep):
    """rep covers member only if BOTH sequences are >=COV covered."""
    c=cov.get((member,rep))
    if c is None:
        c=cov.get((rep,member))
        if c is None: return None
        c=(c[1],c[0])
    return c if (c[0]>=COV and c[1]>=COV) else None

reps=collections.defaultdict(list); decisions=[]
for b,keys in sorted(bagkeys.items()):
    o=orig[b]; reps[b]=[o]
    for m in sorted(keys-{o}, key=lambda k:(-len(seqs[k]),k)):
        hit=None
        for r in reps[b]:
            c=covered(m,r)
            if c: hit=(r,c); break
        if hit: decisions.append((b,m,"absorbed",hit[0],round(hit[1][0],3),round(hit[1][1],3)))
        else:   reps[b].append(m); decisions.append((b,m,"kept_as_variant","",0.0,0.0))
with open(W+"/query_decisions.tsv","w") as o:
    o.write("bag_id\tseq_key\tdecision\tabsorbed_by\tcov_of_member\tcov_of_representative\n")
    for d in decisions: o.write("\t".join(str(x) for x in d)+"\n")
with open(W+"/query_repaired_reps.tsv","w") as o:
    o.write("bag_id\toriginal_rep_as_recorded\tsubstitute\n")
    for d in repaired: o.write("\t".join(d)+"\n")

rows=[]
for b,ks in sorted(reps.items()):
    for i,k in enumerate(ks): rows.append(("Q_%s"%k,b,k,i==0))
seen=set()
with open(W+"/queries.fna","w") as fa:
    for qid,b,k,isorig in rows:
        if k in seen: continue
        seen.add(k); fa.write(">%s\n%s\n"%(qid,seqs[k]))
with open(W+"/queries.tsv","w") as t:
    t.write("query_id\tbag_id\tseq_key\tis_bag_original_rep\tlen_nt\n")
    for qid,b,k,isorig in rows:
        t.write("%s\t%s\t%s\t%d\t%d\n"%(qid,b,k,int(isorig),len(seqs[k])))

# cross-bag equivalence on REPRESENTATIVES ONLY -- recorded, never applied here
os.makedirs(W+"/tmp",exist_ok=True)
subprocess.run([BL+"/makeblastdb","-in",W+"/queries.fna","-dbtype","nucl",
                "-out",W+"/tmp/qdb","-title","q"],check=True,stdout=subprocess.DEVNULL)
subprocess.run([BL+"/blastn","-task","megablast","-query",W+"/queries.fna","-db",W+"/tmp/qdb",
                "-evalue","1e-50","-perc_identity","90","-max_target_seqs","100000",
                "-num_threads","60","-outfmt",
                "6 qseqid sseqid pident length qstart qend sstart send qlen slen",
                "-out",W+"/queries_selfaln.tsv"],check=True)
bag_of={}; 
for qid,b,k,isorig in rows: bag_of.setdefault(qid,b)
eq=collections.defaultdict(list); es=collections.defaultdict(list); ep={}
for l in open(W+"/queries_selfaln.tsv"):
    q,s,p,aln,qs,qe,ss,se,ql,sl=l.rstrip("\n").split("\t")
    if q==s or float(p)<ID: continue
    qs,qe,ss,se=int(qs),int(qe),int(ss),int(se)
    if se<ss: ss,se=se,ss
    eq[(q,s)].append((qs,qe)); es[(q,s)].append((ss,se))
    ep[(q,s)]=max(ep.get((q,s),0),float(p))
equiv=[]
for (q,s) in eq:
    if q>=s: continue
    a=mlen(eq[(q,s)])/len(seqs[q[2:]]); bb=mlen(es[(q,s)])/len(seqs[s[2:]])
    if a>=COV and bb>=COV:
        equiv.append((q,s,bag_of.get(q,""),bag_of.get(s,""),round(a,3),round(bb,3),ep[(q,s)]))
with open(W+"/query_equivalence.tsv","w") as o:
    o.write("query_a\tquery_b\tbag_a\tbag_b\tcov_a\tcov_b\tpident\n")
    for e in sorted(equiv): o.write("\t".join(str(x) for x in e)+"\n")

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
 "original_reps_repaired":len(repaired),
 "cross_bag_equivalent_query_pairs_recorded":len(equiv),
 "cross_bag_equivalent_pairs_in_different_bags":sum(1 for e in equiv if e[2]!=e[3])},indent=2))
PY
