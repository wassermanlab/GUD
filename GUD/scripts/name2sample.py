#!/usr/bin/env python

import argparse
import coreapi
import json
import os

# Import from OnTarget module
from . import OnTargetUtils

usage_msg = """
usage: %s --name STR [-h] [options]
""" % os.path.basename(__file__)

help_msg = """%s
fuzzy search for samples.

  --name STR          sample name for string matching

optional arguments:
  -h, --help          show this help message and exit
  -j, --json          output in JSON format
""" % usage_msg

#-------------#
# Functions   #
#-------------#

def parse_args():
    """
    This function parses arguments provided via the command line and returns an {argparse} object.
    """

    parser = argparse.ArgumentParser(add_help=False)

    # Mandatory args
    parser.add_argument("--name")

    # Optional args
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-h", "--help", action="store_true")
    optional_group.add_argument("-j", "--json", action="store_true")
    
    args = parser.parse_args()

    check_args(args)

    return(args)

def check_args(args):
    """
    This function checks an {argparse} object.
    """

    # Print help
    if args.help:
        print(help_msg)
        exit(0)

    # Check mandatory arguments
    if not args.name:
        error = ["%s\n%s" % (usage_msg, os.path.basename(__file__)), "error", "argument \"--name\" is required\n"]
        print(": ".join(error))
        exit(0)

def main():

    # Parse arguments
    args = parse_args()

    # Fuzzy search for samples
    s = search(args.name, args.json)

    print(s)

def search(name, as_json=False):
    """
    e.g. python -m GUD.scripts.name2samples --name "endothelial cells"
    """

    # Initialize
    answers = []
    client = coreapi.Client()
    codec = coreapi.codecs.CoreJSONCodec()
    params = "name=%s" % name # I assume it's going to be name

    # i.e. expected output for endothelial cells
    tmp_answers = [
        ("endothelial cell (microvasculature)", 90.),
        ("endothelial cell (thorax)", 95.),
        ("endothelial cell (vein)", 97.),
        ("endothelial cell (artery)", 95.),
        ("endothelial cell (aorta)", 96.),
        ("endothelial cell (glomerulus)", 93.),
        ("endothelial cell (renal glomerulus)", 90.),
        ("endothelial cell (pulmonary artery)", 90.),
        ("endothelial cell (lung microvasculature)", 87.),
        ("endothelial cell (umbilical vein)", 91.),
        ("endothelial cell (liver sinusoid)", 91.),
        ("endothelial cell (lymph node)", 93.),
        ("endothelial cell (brain microvasculature)", 86.),
        ("endothelial cell (blood vessel dermis)", 89.),
        ("endothelial cell (microvascular lymphatic vessel dermis)", 79.),
        ("endothelial progenitor cell (derived from CD14-positive monocyte)", 75.)
    ]

    try:

        # ------------------- #
        # Here GUD is queried #
        # ------------------- #

        # Get first page
        response = client.get(os.path.join(OnTargetUtils.gud, parameters))
        page = json.loads(codec.encode(response))

        # While there are more pages...
        while page["url_of_next_page"]:

            # For each feature in page...
            for feat in page["features"]:
                # I don't know how they look like, but should return name + score from fuzzy search
                answers.append(feat.name, feat.score)

            # Go to next page
            response = client.get(page["url_of_next_page"])
            page = json.loads(codec.encode(response))

        # Do last page...
        for feat in page["features"]:
            answers.append((feat.name, feat.score))

    except:

        answers = tmp_answers

    # Sort
    answers.sort(key=lambda x: x[-1], reverse=True)

    if answers:

        if as_json:
            return(json.dumps(answers, indent=4))
        else:
            return("\n".join("%s\t%s" % (a, s) for a, s in answers))

    else:
        raise ValueError("No samples found!!!")

#-------------#
# Main        #
#-------------#

if __name__ == "__main__":

    main()