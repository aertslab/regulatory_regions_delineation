#!/usr/bin/env python3
import sys


def bedfile_iterator(filename):
    with open(filename, "r") as handle:
        for line in handle:
            # Only keep 4th column of BED file with region ID information ...
            region_name = line.rstrip().split("\t")[3]
            gene_name, nr = region_name.split("#")
            yield region_name, gene_name


def main():
    if len(sys.argv) != 2:
        print(
            "Wrong number of input arguments.",
            "Usage: {0:s} bed_filename > gene_name_to_region_names.tsv".format(
                sys.argv[0]
            ),
            sep="\n",
            file=sys.stderr,
        )
        sys.exit(2)

    gene_name2region_names = {}
    for region_name, gene_name in bedfile_iterator(sys.argv[1]):
        gene_name2region_names.setdefault(gene_name, set()).add(region_name)
    for gene_name in sorted(gene_name2region_names.keys()):
        print(gene_name + "\t" + "\t".join(sorted(gene_name2region_names[gene_name])))


if __name__ == "__main__":
    main()
