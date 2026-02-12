#!/usr/bin/env python3

import argparse
from typing import Union

import pysam
from pysam import AlignmentHeader


def retrieve_info_from_header(header: Union[dict, AlignmentHeader]) -> str:
    line_with_info = header['RG'][0]

    results_elements = []
    if 'PL' in line_with_info:
        results_elements.append("PL:"+line_with_info["PL"])
    if 'PU' in line_with_info:
        results_elements.append("PU:"+line_with_info["PU"])
    if 'DS' in line_with_info:
        results_elements.append("DS:"+line_with_info["DS"])
    if 'LB' in line_with_info:
        results_elements.append("LB:"+line_with_info["LB"])

    return "\t".join(results_elements)


def main(input_path: str):
    sam_file = pysam.AlignmentFile(input_path, "rb")
    print(retrieve_info_from_header(sam_file.header))


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('input_path', help='path to bam file to process')
    args = vars(parser.parse_args())
    main(args['input_path'])

