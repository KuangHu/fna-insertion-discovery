"""Local audit: candidate interval groups only, NOT verified events/host loci.

Companion hits table retains all seeds, strand and boundary alternatives.
Deterministic greedy complete-link interval groups are an explicit heuristic.
"""
import argparse
import collections
import csv
import hashlib
import json
from pathlib import Path

csv.field_size_limit(2147483647)


def stable(r):
    return (r['genome'], r['contig'], int(r['elem_start']), int(r['elem_end']),
            json.dumps(r, sort_keys=True))


def iv(r):
    return int(r['elem_start']), int(r['elem_end'])


def same(a, b, cutoff=.9):
    overlap = min(a[1], b[1])-max(a[0], b[0])
    return overlap > 0 and overlap/(a[1]-a[0]) >= cutoff and overlap/(b[1]-b[0]) >= cutoff


def group_hits(rows, cutoff=.9):
    by = collections.defaultdict(list)
    for r in rows:
        s, e = iv(r)
        if s < 0 or e <= s:
            raise ValueError('Invalid interval')
        by[r['genome'], r['contig']].append(r)
    result = []
    for key in sorted(by):
        active = []
        for r in sorted(by[key], key=stable):
            keep = []
            for g in active:
                if max(iv(x)[1] for x in g) <= iv(r)[0]:
                    result.append(g)
                else:
                    keep.append(g)
            active = keep
            matches = [g for g in active if all(same(iv(x), iv(r), cutoff) for x in g)]
            if matches:
                matches[0].append(r)
            else:
                active.append([r])
        result.extend(active)
    return sorted(result, key=lambda g: stable(g[0]))


def representative(g):
    # Display choice, not boundary validation. Final key resolves score ties.
    return min(g, key=lambda r: (-float(r['pct_ident']), -float(r['q_cov']), stable(r)))


def flank_key(r):
    l, q = r.get('left_flank','').upper(), r.get('right_flank','').upper()
    if r.get('flank_complete') != 'yes' or min(len(l),len(q)) < 60:
        return None
    l, q = l[-60:], q[:60]
    if set(l+q)-set('ACGT'):
        return None
    return l, q  # target-contig orientation, no host-homology inference


def summarize(rows):
    groups = group_hits(rows)
    counts = collections.Counter(flank_key(representative(g)) for g in groups)
    records = []
    for i,g in enumerate(groups):
        r = representative(g)
        seeds = sorted(set(x['source_event'] for x in g))
        families = sorted(set(x['cds_cluster'] for x in g))
        starts, ends = [iv(x)[0] for x in g], [iv(x)[1] for x in g]
        k = flank_key(r)
        records.append(dict(site_id='S%06d'%i,genome=r['genome'],contig=r['contig'],
            representative_start=r['elem_start'],representative_end=r['elem_end'],
            representative_seed=r['source_event'],representative_ident=r['pct_ident'],
            representative_qcov=r['q_cov'],coordinate_system='0-based-half-open',
            boundary_evidence='alignment-only',n_hits=len(g),n_seeds=len(seeds),
            n_families=len(families),all_seeds=json.dumps(seeds),all_families=json.dumps(families),
            family_ambiguous='yes' if len(families)>1 else 'no',
            start_min=min(starts),start_max=max(starts),end_min=min(ends),end_max=max(ends),
            boundary_spread_bp=max(starts)-min(starts)+max(ends)-min(ends),
            representative_flank_complete=r.get('flank_complete','.'),
            exact_flank_key=hashlib.sha256('|'.join(k).encode()).hexdigest() if k else '.',
            n_groups_sharing_exact_flank=counts[k] if k else 0,
            representative_left_flank60=r.get('left_flank','')[-60:],
            representative_right_flank60=r.get('right_flank','')[:60],
            strands=','.join(sorted(set(x.get('strand','unknown') for x in g))),
            host_locus_status='not_inferred',pairing_variant_status='not_inferred'))
    return records,groups


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--hits',required=True)
    ap.add_argument('--out',required=True)
    a=ap.parse_args()
    with open(a.hits,newline='') as f:
        reader=csv.DictReader(f,delimiter='\t'); fields=reader.fieldnames; rows=list(reader)
    records,groups=summarize(rows)
    Path(a.out).parent.mkdir(parents=True,exist_ok=True)
    with open(a.out+'.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=list(records[0]) if records else ['site_id'],delimiter='\t')
        w.writeheader();w.writerows(records)
    with open(a.out+'.hits.tsv','w',newline='') as f:
        w=csv.DictWriter(f,fieldnames=['site_id']+fields,delimiter='\t');w.writeheader()
        for record,g in zip(records,groups):
            w.writerows(dict(r,site_id=record['site_id']) for r in g)
    print(json.dumps(dict(input_hits=len(rows),candidate_site_groups=len(groups),
        ambiguous_groups=sum(r['family_ambiguous']=='yes' for r in records),
        exact_flank_classes=len(set(r['exact_flank_key'] for r in records)-{'.'}))))


if __name__=='__main__':
    main()
