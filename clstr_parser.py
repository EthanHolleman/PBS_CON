def read_file(f):
    '''
    Read a clustal style file and return a list of lists. The outer list
    represent one block where the file is wrapped and the inner lists are
    the actual alignment records.
    '''
    with open(f) as clstr_file:
        blocks, i, j = [[]], 0, 0
        lines = clstr_file.readlines()
        while j < len(lines):
            if lines[j][:4] == 'name':
                record = [index.strip()
                          for index in lines[j].split(' ') if index != '']
                blocks[i].append(record)
                j += 1
            else:
                blocks.append([])
                i += 1
                j += 2
    return blocks


def concat_blocks(blocks):
    '''
    For each unique record in the alignment blocks store the header as a key
    in seq_dict dictionary and concatenate the wrapped sequences for that record
    as the value for that seq_dict entry. Return the seq_dict. Requires that
    headers are unique.
    '''
    seq_dict = {}
    for block in blocks:
        for record in block:
            if record[0] in seq_dict:
                seq_dict[record[0]] += record[1]
            else:
                seq_dict[record[0]] = record[1]
    return seq_dict
