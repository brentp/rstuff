import sys
from toolshed import reader
from itertools import izip
from math import log

readers = [reader(f) for f in sys.argv[1:]]


for i, row_set in enumerate(izip(*readers)):
    if i == 0:
        header = ["probe"] + [r['IMAGE_ID'] for r in row_set]
        print "\t".join(header)
    assert len(set(r['SEQ_ID'] for r in row_set)) == 1

    # log normalized
    print row_set[0]['SEQ_ID'] + "\t" + "\t".join("%.3f" % log(float(r['EXPRS']), 2) for r in row_set)
