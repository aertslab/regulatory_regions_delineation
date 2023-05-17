#!/usr/bin/env python
import operator
import sys

INCLUDE_FEATURE_TYPES = (
    "mRNA",
    "miRNA",
    "ncRNA",
    "tRNA",
    "rRNA",
    "snRNA",
    "snoRNA",
    "CDS",
    "exon",
)


def iterate_gff_body(filename):
    with open(filename, "r") as handle:
        for line in handle:
            if line.startswith("###"):
                return
            if line.startswith("##"):
                continue
            (
                seq_id,
                source,
                feature_type,
                start,
                end,
                score,
                strand,
                phase,
                attributes_str,
            ) = line.rstrip().split("\t")
            if feature_type not in INCLUDE_FEATURE_TYPES:
                continue
            attributes = {}
            for element in attributes_str.split(";"):
                key, value = element.split("=")
                if key == "Parent":
                    attributes[key] = value.split(",")
                elif key == "ID" or key == "Name":
                    attributes[key] = value
            # Make intervals half-open and zero-based (in GFF they are closed intervals and 1-based) ...
            yield feature_type, seq_id, int(start) - 1, int(end), strand, attributes


class Entry:
    @staticmethod
    def from_mRNA(feature_id, chromosome, tx, strand, gene_id):
        entry = Entry(feature_id)
        entry.chromosome, entry.tx, entry.strand, entry.gene_id = (
            chromosome,
            tx,
            strand,
            gene_id,
        )
        return entry

    def __init__(self, feature_id):
        (
            self.feature_id,
            self.chromosome,
            self.tx,
            self.strand,
            self.CDSs,
            self.exons,
            self.gene_id,
        ) = (feature_id, "", (), "", [], [], "")

    def isempty(self):
        return len(self.tx) == 0

    def exon_count(self):
        return len(self.exons)

    def exon_starts(self):
        if self.exon_count() == 0:
            return ()
        return sorted(map(operator.itemgetter(0), self.exons))

    def exon_ends(self):
        if self.exon_count() == 0:
            return ()
        return sorted(map(operator.itemgetter(1), self.exons))

    def cds_span(self):
        if len(self.CDSs) == 0:
            return self.tx[1], self.tx[1]
        return min(map(operator.itemgetter(0), self.CDSs)), max(
            map(operator.itemgetter(1), self.CDSs)
        )


def load_data(gff_filename):
    feature_id2entry = {}
    for feature_type, chromosome, start, end, strand, attributes in iterate_gff_body(
        gff_filename
    ):
        if feature_type.endswith("RNA"):
            # Analyze attributes ...
            if "ID" not in attributes:
                print(
                    "'{0:s}' entry ({1:s} {2:d} {3:d} {4:s}) has no ID attribute.".format(
                        feature_type, chromosome, start, end, strand
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)
            feature_id = attributes["ID"]
            if "Parent" in attributes and len(attributes["Parent"]) == 1:
                gene_id = attributes["Parent"][0]
            elif "Name" in attributes:
                gene_id = attributes["Name"]
            else:
                print(
                    "'{0:s}' entry (ID: {1:s}) cannot be associated with a geneID.".format(
                        feature_type, feature_id
                    ),
                    file=sys.stderr,
                )
                sys.exit(1)

            if feature_id in feature_id2entry:
                # if feature_id2entry[feature_id].chromosome != chromosome or feature_id2entry[feature_id].strand != strand:
                #     print("'mRNA' entry not compatible with child entry.", file=sys.stderr)
                #     sys.exit(1)
                feature_id2entry[feature_id].chromosome = chromosome
                feature_id2entry[feature_id].tx = (start, end)
                feature_id2entry[feature_id].strand = strand
                feature_id2entry[feature_id].gene_id = gene_id
            else:
                feature_id2entry[feature_id] = Entry.from_mRNA(
                    feature_id, chromosome, (start, end), strand, gene_id
                )
        elif feature_type == "exon":
            if "Parent" not in attributes:
                print(
                    f"'exon' entry (ID: {feature_id:s}) has no Parent attribute.",
                    file=sys.stderr,
                )
                sys.exit(1)
            exon = (start, end)
            for parent_id in attributes["Parent"]:
                if parent_id not in feature_id2entry.keys():
                    feature_id2entry[parent_id] = Entry(parent_id)
                # elif feature_id2entry[parent_id].chromosome != chromosome or feature_id2entry[parent_id].strand != strand:
                #     print("'exon' entry not compatible with parent 'mRNA' entry {0:s}.".format(parent_id),
                #           file=sys.stderr)
                #     sys.exit(1)
                feature_id2entry[parent_id].exons.append(exon)
                # feature_id2entry[parent_id].chromosome = chromosome
                # feature_id2entry[parent_id].strand = strand
        elif feature_type == "CDS":
            if "Parent" not in attributes:
                print(
                    f"'CDS' entry (ID: {feature_id:s}) has no Parent attribute.",
                    file=sys.stderr,
                )
                sys.exit(1)
            CDS = (start, end)
            for parent_id in attributes["Parent"]:
                if parent_id not in feature_id2entry.keys():
                    feature_id2entry[parent_id] = Entry(parent_id)
                # elif feature_id2entry[parent_id].chromosome != chromosome or feature_id2entry[parent_id].strand != strand:
                #     print("'CDS' entry not compatible with parent 'mRNA' entry {0:s}.".format(parent_id),
                #           file=sys.stderr)
                #     sys.exit(1)
                feature_id2entry[parent_id].CDSs.append(CDS)
                # feature_id2entry[parent_id].chromosome = chromosome
                # feature_id2entry[parent_id].strand = strand
    return feature_id2entry


def main():
    if len(sys.argv) != 2:
        print("Wrong number of input arguments.", file=sys.stderr)
        sys.exit(2)
    gff_filename = sys.argv[1]
    print("Loading GFF into memory ...", file=sys.stderr)
    feature_id2entry = load_data(gff_filename)
    print("Conversion to tbl format ...", file=sys.stderr)
    print(
        "#bin",
        "name",
        "chrom",
        "strand",
        "txStart",
        "txEnd",
        "cdsStart",
        "cdsEnd",
        "exonCount",
        "exonStarts",
        "exonEnds",
        "score",
        "name2",
        "cdsStartStat",
        "cdsEndStat",
        "exonFrames",
        sep="\t",
    )

    for feature_id, entry in feature_id2entry.items():
        if entry.isempty():
            print(f"{feature_id:s} is an empty entry.", file=sys.stderr)
            continue
        if entry.exon_count() == 0:
            entry.exons.append(entry.tx)
        cds = entry.cds_span()
        exon_starts = ",".join(map(str, entry.exon_starts())) + ","
        exon_ends = ",".join(map(str, entry.exon_ends())) + ","
        print(
            "NA",
            feature_id,
            entry.chromosome,
            entry.strand,
            entry.tx[0],
            entry.tx[1],
            cds[0],
            cds[1],
            entry.exon_count(),
            exon_starts,
            exon_ends,
            "NA",
            entry.gene_id,
            "NA",
            "NA",
            "NA",
            sep="\t",
        )


if __name__ == "__main__":
    main()
