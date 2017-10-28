#!/usr/bin/env python

import re
import argparse


parser = argparse.ArgumentParser(description="Annotation script")

parser.add_argument("-in", dest="infile", help="Input gff file")
parser.add_argument("-gff", dest="gff_file", help="Reference gff")

args = parser.parse_args()


def parse_reference_gff(gff_file):
    """Parses a reference gff file and returns a dictionary with the
    annotations

    :gff_file: TODO
    :returns: TODO
    """

    annotations = {}

    with open(gff_file) as fh:

        for line in fh:

            if line.strip() == "" or line.startswith("#"):
                continue

            fields = line.split("\t")

            # Skip non polypeptide lines
            if fields[2] != "polypeptide":
                continue

            line = line.strip()

            # Get gene id
            metadata = fields[-1]
            gene_id = metadata.split(";")[0].split("=")[1]

            # Get go terms
            go_terms = re.findall(r"GO:([0-9]*),", fields[-1])

            # Get range of polypeptide
            gene_range = [int(fields[3]), int(fields[4])]

            # Get chromosome
            chrom = fields[0]

            annotations[gene_id] = {
                    "go": go_terms,
                    "range": gene_range,
                    "chr": chrom
                    }

    return annotations
            

def parse_target_gff(gff_file):

    gff_storage = {}

    with open(gff_file) as fh:

        for line in fh:

            if line.startswith("#") or line.strip() == "":
                continue

            line = line.strip()

            fields = line.split("\t")

            if fields[2] != "transcript":
                continue

            # Get gene id
            gene_id = re.findall(r"gene_id \"MSTRG\.([0-9]*)", fields[-1])[0]

            # Get range
            gene_range = [int(fields[3]), int(fields[4])]

            # Get chrom
            chrom = fields[0]

            if gene_id in gff_storage:

                if gene_range[0] < gff_storage[gene_id]["range"][0]:
                    gff_storage[gene_id]["range"][0] = gene_range[0]

                if gene_range[1] > gff_storage[gene_id]["range"][1]:
                    gff_storage[gene_id]["range"][1] = gene_range[1]

            else:

                gff_storage[gene_id] = {
                        "range": gene_range,
                        "chr": chrom
                        }

    return gff_storage


def test_overlap(rg1, rg2):

    # x = range(rg1[0], rg1[1])
    # y = range(rg2[0], rg2[1])

    # res = set(x).intersection(y)

    if rg2[0] < rg1[0] < rg2[1] or rg2[0] < rg1[1] < rg2[1]:
        return True
    else:
        return False


def get_annotations(gff_storage, gff_reference):

    annotations = {}

    for gene_id, vals in gff_storage.items():

        print("Scanning gene {}".format(gene_id), end="\r")

        for ref_id, metadata in gff_reference.items():

            if vals["chr"] != metadata["chr"]:
                continue

            if vals["range"][1] < metadata["range"][0]:
                continue

            if test_overlap(vals["range"], metadata["range"]):
                annotations[gene_id] = metadata["go"]
                break
        
        # If no overlapp was found, populate with  none
        if gene_id not in annotations:
            annotations[gene_id] = None

    print(len(annotations))
                

def main():
    """Main execution function

    """
    
    infile = args.infile
    gff_file = args.gff_file

    ref = parse_reference_gff(gff_file)
    target = parse_target_gff(infile)

    get_annotations(target, ref)


main()
