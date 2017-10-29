#!/usr/bin/env python

import re
import plotly
import argparse


parser = argparse.ArgumentParser(description="Annotation script")

parser.add_argument("-in", dest="infile", help="Input gff file")
parser.add_argument("-gff", dest="gff_file", help="Reference gff")
parser.add_argument("-go", dest="go_file", help="GO terms file")
parser.add_argument("-inter", dest="inter_file", help="Interpro2GO mapping")

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
            
            # Get interpro if no go terms are found
            if not go_terms:
                interpro = re.findall(r"InterPro:(IPR[0-9]*),", fields[-1])
            else:
                interpro = None

            # Get range of polypeptide
            gene_range = [int(fields[3]), int(fields[4])]

            # Get chromosome
            chrom = fields[0]

            annotations[gene_id] = {
                    "go": go_terms,
                    "interpro": interpro,
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


def parse_interpro(interpro_file):

    inter_map = {}

    with open(interpro_file) as fh:

        for line in fh:
            if line.startswith("InterPro:"):
                inter_id = re.findall(r"InterPro:(IPR[0-9]*) ",
                                      line.strip())[0]
                go_term = re.findall(r"GO:([0-9]*)", line.strip())[1]

                inter_map[inter_id] = go_term

    return inter_map


def test_overlap(rg1, rg2):

    if rg2[0] < rg1[0] < rg2[1] or rg2[0] < rg1[1] < rg2[1] or \
            rg1[0] < rg2[0] < rg1[1] or rg1[0] < rg2[1] < rg1[1]:
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

            if test_overlap(vals["range"], metadata["range"]):
                if metadata["go"]:
                    annotations[gene_id] = metadata["go"]
                else:
                    annotations[gene_id] = metadata["interpro"]
                break
 
        # If no overlapp was found, populate with  none
        if gene_id not in annotations:
            annotations[gene_id] = None

    return annotations


def parse_go_terms(go_file):

    go_storage = {}

    go_id = None
    go_cat = None
    go_namespace = None
    go_is_a = {}

    with open(go_file) as fh:

        for line in fh:

            if line.startswith("[Term]"):

                if go_id and go_cat:
                    go_storage[go_id] = {
                            "cat": go_cat,
                            "namespace": go_namespace,
                            "is_a": go_is_a,
                            }

                go_id = None
                go_cat = None
                go_namespace = None
                go_is_a = {}

            try:
                if line.startswith("id:"):
                    go_id = line.strip().split()[1].split(":")[1]

                if line.startswith("name:"):
                    go_cat = line.split(":")[1].strip()

                if line.startswith("namespace:"):
                    go_namespace = line.split(":")[1].strip()

                if line.startswith("is_a:"):
                    is_a_temp = line.strip().split(": ")[1]
                    is_a = is_a_temp.split("!")
                    go_is_a[is_a[0].split(":")[1].strip()] = is_a[1].strip()
            except IndexError:
                pass

    return go_storage


def get_go_hierarchy(annotation_dic, go_terms):

    pass



def main():
    """Main execution function

    """
    
    infile = args.infile
    gff_file = args.gff_file
    go_file = args.go_file
    inter_file = args.inter_file

    ref = parse_reference_gff(gff_file)
    target = parse_target_gff(infile)
    inter_map = parse_interpro(inter_file)

    print(list(inter_map.items())[:10])

    annotations = get_annotations(target, ref)
    print(list(annotations.items())[:10])
    check = [x for x in annotations.values() if not x]
    print(len(check))

    go_terms = parse_go_terms(go_file)
    print(list(go_terms.items())[:10])


main()
