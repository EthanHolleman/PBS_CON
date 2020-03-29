from clstr_parser import read_file, concat_blocks
import os
import graphviz


def find_PBS_canidates(seq_dict, ltr_len=800, reach=1000, start_motif='TGG',
                       seq_len=12, allow_overlap=True):

    # slice sequences in the seq_dict at given indicies (ltr_len is start)
    seq_dict = {key: seq_dict[key][ltr_len:ltr_len+reach] for key in seq_dict}

    # find start motifs in each slice and store as new dictionary
    # called canidate_start_sites key is still record header
    canidate_start_sites = {head: find_start_motif_sites(
        start_motif, seq_dict[head]) for head in seq_dict}

    if not allow_overlap:  # remove overlaps if not allowed
        for head in canidate_start_sites:
            canidates, trimmed, i = canidate_start_sites[head], [], 0
            if canidates:
                for i in range(1, len(canidates)):
                    if canidates[i] > canidates[i-1]+seq_len:
                        trimmed.append(canidates[i-1])
                if canidates[-1] > canidates[-2]+seq_len:
                    trimmed.append(canidates[-1])

    # want to find if there is a generally conserved seqeucne
    canidate_seqs = {}
    for head in seq_dict:
        canidate_seqs[head] = [seq_dict[head][start:start+seq_len]
                               for start in canidate_start_sites[head]]

    occurance_lists = []
    max_canidates = max([len(canidate_seqs[head]) for head in canidate_seqs])

    for i in range(max_canidates):
        occurance = []
        for head in canidate_seqs:
            if len(canidate_seqs[head]) > i:
                occurance.append(canidate_seqs[head][i])
        occurance_lists.append(occurance)

    return [make_score_consensus(occur_list) for occur_list in occurance_lists]


def make_score_consensus(occurance_list):
    BASES, score, con = ['A', 'T', 'G', 'C'], 0, ''
    for i in range(len(occurance_list[0])):
        at_position_i = [l[i] for l in occurance_list]
        base_counts = {at_position_i.count(B): B for B in BASES}
        s = max(base_counts.keys())
        con += base_counts[s]
        score += s
    return score, con, occurance_list


def find_start_motif_sites(start_motif, seq):
    '''
    kmer search based on given start motif. Returns a list of indicies where
    the start motif can be found in the given sequence.
    '''
    canidate_sites = []
    for i in range(len(seq)-len(start_motif)+1):
        kmer = seq[i:i+len(start_motif)]
        if kmer == start_motif:
            canidate_sites.append(i)
    return canidate_sites


def write_consensus_seqs(score_con_list, out_dir, filename='con_seqs.fasta'):
    # score con list is list of tuples consensus score then consensus seq
    with open(os.path.join(out_dir, filename), 'w') as out:
        i = 0
        for score, con in sorted(score_con_list, key=lambda x: x[0], reverse=True):
            out.write(f'>Canidate_{i}_{score}\n{con}\n')
            i+=1



def graph(sequence_list):
    BASES = ['A', 'T', 'G', 'C', '-']
    MIN = 0.15

    g = graphviz.Graph(format='png')
    structure = []
    for i in range(len(sequence_list[0])):  # assume all same length for now
        at_position_i = [l[i] for l in sequence_list]
        counts = {B: at_position_i.count(B) for B in BASES}
        counts = {B: counts[B] / len(at_position_i)
                  for B in counts if counts[B] / len(at_position_i) > MIN}
        structure.append(counts)

    nodes = []
    for j, layer in enumerate(structure):
        nodes.append([])
        for b in layer:
            if nodes[j]:
                nodes[j].append([b, layer[b], 0])
            else:
                nodes[j] = [[b, layer[b], 0]]

    for i in range(len(nodes)):
        max_node, v = 0, 0
        for j in range(len(nodes[i])):
            if nodes[i][j][1] > v:
                v = nodes[i][j][1]
                max_node = j
        nodes[i][max_node][-1] = 1

    # draw the nodes
    for p in range(len(nodes)):
        cur_nodes = nodes[p]
        for l in range(len(cur_nodes)):
            label, size, m = cur_nodes[l]
            print(label, size, m)
            if m == 0:
                g.node(f'{p}{l}', label, **
                       {'width': str(size), 'height': str(size)})
            else:
                g.node(f'{p}{l}', label, style='filled', fillcolor='lightblue',
                       **{'width': str(size), 'height': str(size)})

    for p in range(len(nodes)):
        if p < len(nodes) - 1:
            cur_nodes = nodes[p]
            for l in range(len(cur_nodes)):
                label, size, m = cur_nodes[l]
                for x in range(len(nodes[p+1])):
                    g.edge(f'{p}{l}', f'{p+1}{x}')

    g.render('./top_canidate', view=True)
