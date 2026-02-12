#!/usr/bin/env python3

import argparse
from typing import Set, List, Dict

from Bio.Seq import Seq
from Bio.SeqIO.FastaIO import FastaWriter
from Bio.SeqIO.FastaIO import SeqRecord
from pysam import AlignmentFile
from Bio import SeqIO
from pyfaidx import Fasta


def get_sequences_needed(bam_path: str, current_set: List[str]) -> None:
    if len(current_set) == 0:
        with AlignmentFile(bam_path) as input:
            for reference in input.header.references:
                current_set.append(reference)
    else:
        with AlignmentFile(bam_path) as input:
            count: int = 0
            for reference in input.header.references:
                if count < len(current_set):
                    if reference != current_set[count]:
                        exit(1)
                else:
                    current_set.append(reference)
                count  += 1


def get_sequences_available(fasta_path: str) -> Set[str]:
    fasta_input = SeqIO.parse(fasta_path, "fasta")
    output: Set[str] = set()
    for sequence in fasta_input:
        output.add(sequence.id)
    return output


def copy_content(output_fasta: FastaWriter, base_fasta_path: str, sequence_name : str = None):
    if sequence_name is None:
        fasta_input = SeqIO.parse(base_fasta_path, "fasta")
        for sequence in fasta_input:
            output_fasta.write_record(sequence)
    else:
        sequences = Fasta(base_fasta_path)
        sequence = sequences[sequence_name]
        record: SeqRecord = SeqRecord(Seq(str(sequence[:])), id=sequence.name, name=sequence.name, description=sequence.long_name)
        output_fasta.write_record(record)

def parse_contig_mappings(contig_mappings_path: str) -> Dict[str, str]:
    result: Dict[str, str] = {}
    with open(contig_mappings_path, "r") as contig_mappings:
        for contig_mapping in contig_mappings:
            elements = contig_mapping.split()
            result[elements[0].strip()] = elements[1].strip()
    return result


def create_expanded_reference(bam_paths: List[str], base_fasta_path: str, mappings_path: str, output_file: str ):
    sequences_needed: List[str] = list()
    for bam_path in bam_paths:
        get_sequences_needed(bam_path, sequences_needed)

    sequences_available: Set[str] = get_sequences_available(base_fasta_path)

    mappings = parse_contig_mappings(mappings_path)
    for sequence_available in sequences_available:
        mappings[sequence_available] = base_fasta_path

    with open(output_file, "w") as output:
        writer = FastaWriter(output)

        for sequence_needed in sequences_needed:
            if sequence_needed not in mappings.keys():
                exit(1)
            else:
                copy_content(writer, mappings[sequence_needed], sequence_needed)



if __name__ == '__main__':
    ap = argparse.ArgumentParser(description="Adds to a base fasta file, new sequences from a database, in order to cover"
                                             "all the sequences needed by the input bams")

    ap.add_argument("-b", "--bams", nargs='+', required=True)
    ap.add_argument("-f", "--fasta", required=True)
    ap.add_argument("-m", "--mappings", required=True)
    ap.add_argument("-o", "--output", required=True)

    args = vars(ap.parse_args())

    create_expanded_reference(args['bams'], args['fasta'], args['mappings'], args['output'])
