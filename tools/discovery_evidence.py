"""Evidence helpers shared by the discovery runner and local comparisons."""
import collections


def shared_intervals(graph, names):
    """Conservative two-path comparison; excludes wraps and reordered anchors."""
    if len(names)!=2 or len(set(names))!=2:
        raise ValueError('Exactly two distinct sequence names required')
    paths={p['name']:p for p in graph['paths'].values()}
    walks={name:sorted((graph['nodes'][str(i)] for i in paths[name]['nodes']),
                       key=lambda n:n['position'][0]) for name in names}
    counts={n:collections.Counter(x['block_id'] for x in w) for n,w in walks.items()}
    core={b for b,c in counts[names[0]].items() if c==1 and counts[names[1]][b]==1}
    a,b=([n for n in walks[name] if n['block_id'] in core] for name in names)
    index={n['block_id']:i for i,n in enumerate(b)}
    rows=[];excluded=0
    for left,right in zip(a,a[1:]):
        i,j=index[left['block_id']],index[right['block_id']]
        if j!=i+1 or left['strand']!=b[i]['strand'] or right['strand']!=b[j]['strand']:
            excluded+=1;continue
        s,e=left['position'][1],right['position'][0]
        u,v=b[i]['position'][1],b[j]['position'][0]
        if e<s or v<u:
            excluded+=1;continue
        delta=(v-u)-(e-s)
        rows.append(dict(id='P%04d'%len(rows),reference_contig=names[0],reference_start=s,
                         reference_end=e,query_contig=names[1],query_start=u,query_end=v,
                         reference_bp=e-s,query_bp=v-u,signed_delta=delta,
                         category='primary' if 500<=abs(delta)<=5000 else 'large' if abs(delta)>5000 else 'below_500',
                         left_block=left['block_id'],right_block=right['block_id']))
    deletions=[d['len'] for block in graph['blocks'].values() for aln in block['alignments'].values() for d in aln['dels']]
    insertions=[len(d['seq']) for block in graph['blocks'].values() for aln in block['alignments'].values() for d in aln['inss']]
    return rows,dict(blocks=len(graph['blocks']),single_copy_shared_blocks=len(core),
        excluded_order_or_wrap=excluded,origin_crossing_interval_not_tested=1,
        max_within_block_deletion=max(deletions,default=0),max_within_block_insertion=max(insertions,default=0))


SEEDED = ('seeded',)


def padded_seed(event, lengths, pad=200):
    """Return (seed, 'seeded') or (None, reason); never drop or crash on a candidate.

    SPEC M1-B: a seed is junction +/- pad, and Arm B gives junction_L/R only
    when the empty reference's own allele was placed exactly
    (frame_status == 'verified', B-6). Any other frame has no verified
    junction, so the candidate is kept as unseedable with the frame as its
    reason (B-8); the raw call interval stays in the Arm B table for anyone
    who wants to look, but it is not a junction and is not used as one."""
    genome, contig = event.get('empty_ref_genome', ''), event.get('empty_ref_contig', '')
    if genome not in lengths or contig not in lengths[genome]:
        return None, 'genome_or_contig_not_in_manifest'
    length = lengths[genome][contig]
    frame = event.get('frame_status')
    if frame != 'verified':
        return None, 'no_frame_status' if frame is None else 'frame_' + frame
    try:
        left, right = int(event['junction_L']), int(event['junction_R'])
    except (KeyError, ValueError):
        return None, 'invalid_junction_frame'
    if left < 0 or right < left or right > length:
        return None, 'invalid_junction_frame'
    start, end = left - pad, right + pad
    if start < 0 or end > length:
        return None, 'seed_padding_outside_contig'
    return dict(locus_id=event['event_id'], seed_genome=genome, contig=contig,
                start=start, end=end), 'seeded'


def candidate_row(event, seed_status, reconstructed):
    """Preserve all Layer 1 candidates, with explicit Layer 2 availability."""
    r=dict(event)
    r['layer2_status']=('unseedable' if seed_status not in SEEDED else
                        'unresolved' if reconstructed is None else
                        'ambiguous' if reconstructed['placement_status']=='ambiguous' else 'resolved')
    r['seed_status']=seed_status
    for key in ('event_class','inserted_len','inserted_md5','junction_L_abs','junction_R_abs',
                'junction_ambiguity_bp','target_bases_lost','placement_status','decomposition_method'):
        r['layer2_'+key]=reconstructed.get(key,'.') if reconstructed else '.'
    return r
