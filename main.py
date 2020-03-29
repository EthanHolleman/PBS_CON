from args import get_args
from clstr_parser import concat_blocks, read_file
from finder_grapher import find_PBS_canidates, graph, write_consensus_seqs

def main():
    a = get_args()
    blocks = concat_blocks(read_file(a.c))
    results = find_PBS_canidates(blocks, int(a.s), int(a.e), a.m)
    results = sorted(results, key=lambda x: x[0], reverse=True)  # sort by score
    score_con_list = [(r[0], r[1]) for r in results]
    write_consensus_seqs(score_con_list, a.o)
    graph(results[0][2])


if __name__ == "__main__":
    main()
    
    
    
    