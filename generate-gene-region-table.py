#!/usr/bin/env python
import sys


def bedfile_iterator(filename):
    with open(filename, 'r') as handle:
        for line in handle:
             # Only keep 4th column of BED file with ID information ...
            region_name = line.rstrip().split('\t')[3]
            id, nr = region_name.split('#')
            yield region_name, id


def main():
    if len(sys.argv) != 2:
        print "Wrong number of input arguments."
        sys.exit(2)
        
    id2regions = dict()
    for region, id in bedfile_iterator(sys.argv[1]):
        id2regions.setdefault(id, set()).add(region)
    for id in sorted(id2regions.keys()):
        print id + "\t" + "\t".join(sorted(id2regions[id])) 


if __name__ == "__main__":
    main()