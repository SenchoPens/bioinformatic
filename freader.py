from itertools import islice


def read_fa(file):
    next(file)
    res = []
    cur_line = ''
    for line in file:
        if not line.startswith('>'):
            cur_line += line
        else:
            res.append(cur_line.replace('\n', ''))
    return res

def read_fq(file):
    for read in islice(file, 1, None, 4):
        yield read.replace('\n', '')

def read_fq2(file):
    quality = False
    for line in file:
        starts_quality = line.startswith('+')
        name = line.startswith('@')
        quality = (starts_quality or quality) and not name
        if not (quality or name):
            yield line.replace('\n', '')

def read_nex(file):
    res = []
    for line in file:
        res.append(line[line.find(' '):])
    return res


