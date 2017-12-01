#!/usr/bin/env python

import re
import matplotlib.pyplot as plt
import argparse
import pickle


parser = argparse.ArgumentParser(description="Annotation script")

parser.add_argument("-in", dest="infile", help="Input gff file")
parser.add_argument("-gff", dest="gff_file", help="Reference gff")
parser.add_argument("-go", dest="go_file", help="GO terms file")
parser.add_argument("-inter", dest="inter_file",
                    help="Interpro2GO mapping")
parser.add_argument("-annot", dest="annotations", help="Pickle with "
                    "annotations dict")

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

            # Get product
            product = re.findall(r"PRODUCT:(.*),CC_f", fields[-1])[0]

            # Get range of polypeptide
            gene_range = [int(fields[3]), int(fields[4])]

            # Get chromosome
            chrom = fields[0]

            annotations[gene_id] = {
                    "go": go_terms,
                    "interpro": interpro,
                    "product": product,
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
        annotations[gene_id] = {}

        for ref_id, metadata in gff_reference.items():

            if vals["chr"] != metadata["chr"]:
                continue

            if test_overlap(vals["range"], metadata["range"]):

                annotations[gene_id]["product"] = metadata["product"]
                annotations[gene_id]["chr"] = metadata["chr"]
                annotations[gene_id]["range"] = vals["range"]
                if metadata["go"]:
                    annotations[gene_id]["annot"] = metadata["go"]
                else:
                    annotations[gene_id]["annot"] = metadata["interpro"]
                break
 
        # If no overlap was found, populate with  none
        if not annotations[gene_id]:
            annotations[gene_id] = {"annot": None,
                                    "chr": metadata["chr"],
                                    "product": metadata["product"]}

    return annotations


def parse_go_terms(go_file):

    go_storage = {}

    go_id = None
    go_cat = None
    go_namespace = None
    go_alt_id = []
    go_cons = None
    go_is_a = {}

    with open(go_file) as fh:

        for line in fh:

            if line.startswith("[Term]"):

                if go_id and go_cat:
                    go_storage[go_id] = {
                            "cat": go_cat,
                            "namespace": go_namespace,
                            "is_a": go_is_a,
                            "alt_id": go_alt_id,
                            "cons": go_cons
                            }

                go_id = None
                go_cat = None
                go_namespace = None
                go_cons = None
                go_alt_id = []
                go_is_a = {}

            try:
                if line.startswith("id:"):
                    go_id = line.strip().split()[1].split(":")[1]

                if line.startswith("alt_id:"):
                    go_alt_id.append(line.strip().split(":")[-1])

                if line.startswith("name:"):
                    go_cat = line.split(":")[1].strip()

                if line.startswith("namespace:"):
                    go_namespace = line.split(":")[1].strip()

                if line.startswith("consider:"):
                    go_cons = line.strip().split(":")[-1]

                if line.startswith("replaced_by:"):
                    go_cons = line.strip().split(":")[-1]

                if line.startswith("is_a:"):
                    is_a_temp = line.strip().split(": ")[1]
                    is_a = is_a_temp.split("!")
                    go_is_a[is_a[0].split(":")[1].strip()] = is_a[1].strip()
            except IndexError:
                pass

    return go_storage


def get_go_hiearchy(go, go_terms):

    scan = True
    go_tree = [go_terms[go]["cat"]]

    try:
        go_id = list(go_terms[go]["is_a"].keys())[0]
    except IndexError:
        return None

    while scan:

        if not go_terms[go_id]["is_a"]:
            scan = False
            break

        go_tree.append(go_terms[go_id]["cat"])
        go_id = list(go_terms[go_id]["is_a"].keys())[0]

    go_tree.append(go_terms[go]["namespace"])

    return go_tree


def get_alt_go(go, go_terms):

    for k, v in go_terms.items():
        if go in v["alt_id"]:
            return k


def get_go_annotation(annotation_dic, go_terms, inter_map):

    go_terms_dic = {}
    bad = 0

    for gene_id, dic in annotation_dic.items():

        annot = dic["annot"]

        print("Scanning gene {}".format(gene_id), end="\r")
        go_terms_dic[gene_id] = []

        if not annot:
            continue

        for i in annot:

            try:
                if i.startswith("IPR"):
                    i = inter_map[i]
            except KeyError:
                bad += 1
                continue

            # Check if the GO id is in the mapping. If not, check the
            # alternative ids
            if i not in go_terms:
                i = get_alt_go(i, go_terms)

            # Some GO ids may be absolte. Check for this, and update the id
            if go_terms[i]["cons"]:
                i = go_terms[i]["cons"]

            go_tree = get_go_hiearchy(i, go_terms)

            if go_tree:
                go_terms_dic[gene_id].append({
                        "product": dic["product"],
                        "chr": dic["chr"],
                        "range": dic["range"],
                        "go_id": i,
                        "go_tree": go_tree
                        })

    return go_terms_dic


def write_annotations(final_annotation):

    with open("annotations.pck", "wb") as fh:
        pickle.dump(final_annotation, fh)


def load_annotations(path):

    with open(path, "rb") as fh:
        return pickle.load(fh)


def missing_annotations(final_annotation):

    missing = len([x for x in final_annotation.values() if not x])


def histogram_go_terms(final_annotation, level=-1):

    data = [x["go_tree"][level] for x in final_annotation.values() if x]


def main():
    """Main execution function

    """
    
    infile = args.infile
    gff_file = args.gff_file
    go_file = args.go_file
    inter_file = args.inter_file
    annot = args.annotations

    if not annot:
        ref = parse_reference_gff(gff_file)
        target = parse_target_gff(infile)
        inter_map = parse_interpro(inter_file)
        go_terms = parse_go_terms(go_file)

        annotations = get_annotations(target, ref)

        final_annotations = get_go_annotation(annotations, go_terms, inter_map)
        print(list(final_annotations.items())[:10])
        write_annotations(final_annotations)
    else:
        final_annotations = load_annotations(annot)

    # missing_annotations(final_annotations)
    # histogram_go_terms(final_annotations, -2)


main()
