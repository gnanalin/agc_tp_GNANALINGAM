#!/bin/env python3
# -*- coding: utf-8 -*-
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#    A copy of the GNU General Public License is available at
#    http://www.gnu.org/licenses/gpl-3.0.html

"""OTU clustering"""

import argparse
import sys
import os
import gzip
import statistics
import textwrap
from pathlib import Path
from collections import Counter
from typing import Iterator, Dict, List
# https://github.com/briney/nwalign3
# ftp://ftp.ncbi.nih.gov/blast/matrices/
import nwalign3 as nw

__author__ = "Stéphanie Gnanalingam"
__copyright__ = "Universite Paris Cité"
__credits__ = ["Stéphanie Gnanalingam"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Stéphanie Gnanalingam"
__email__ = "stephanie.gnanalingam@etu.u-paris.fr"
__status__ = "Developpement"



def isfile(path: str) -> Path:  # pragma: no cover
    """Check if path is an existing file.

    :param path: (str) Path to the file

    :raises ArgumentTypeError: If file does not exist

    :return: (Path) Path object of the input file
    """
    myfile = Path(path)
    if not myfile.is_file():
        if myfile.is_dir():
            msg = f"{myfile.name} is a directory."
        else:
            msg = f"{myfile.name} does not exist."
        raise argparse.ArgumentTypeError(msg)
    return myfile


def get_arguments(): # pragma: no cover
    """Retrieves the arguments of the program.

    :return: An object that contains the arguments
    """
    # Parsing arguments
    parser = argparse.ArgumentParser(description=__doc__, usage=
                                     "{0} -h"
                                     .format(sys.argv[0]))
    parser.add_argument('-i', '-amplicon_file', dest='amplicon_file', type=isfile, required=True, 
                        help="Amplicon is a compressed fasta file (.fasta.gz)")
    parser.add_argument('-s', '-minseqlen', dest='minseqlen', type=int, default = 400,
                        help="Minimum sequence length for dereplication (default 400)")
    parser.add_argument('-m', '-mincount', dest='mincount', type=int, default = 10,
                        help="Minimum count for dereplication  (default 10)")
    parser.add_argument('-o', '-output_file', dest='output_file', type=Path,
                        default=Path("OTU.fasta"), help="Output file")
    return parser.parse_args()


def read_fasta(amplicon_file: Path, minseqlen: int) -> Iterator[str]:
    """Read a compressed fasta and extract all fasta sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :return: A generator object that provides the Fasta sequences (str).
    """
    with gzip.open(amplicon_file, "rt") as  file_read:
        sequence = ""
        for line in file_read:
            line = line.strip()
            if line.startswith(">"):
                if len(sequence) >= minseqlen:
                    yield sequence
                sequence = ""
            else:
                sequence += line
        yield sequence



def dereplication_fulllength(amplicon_file: Path, minseqlen: int, mincount: int) -> Iterator[List]:
    """Dereplicate the set of sequence

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length
    :param mincount: (int) Minimum amplicon count
    :return: A generator object that provides a (list)[sequences, count] of sequence with a count >= mincount and a length >= minseqlen.
    """
    sequences = list(read_fasta(amplicon_file, minseqlen))
    dict_occ_sequences = {}
    for sequence in sequences:
        dict_occ_sequences[sequence] = dict_occ_sequences.get(sequence, 0) + 1
    list_occ_sequences = sorted(dict_occ_sequences.items(), key=lambda x: x[1], reverse=True)
    for sequence, occurrence in list_occ_sequences:
        if occurrence >= mincount:
            yield [sequence, occurrence]


def get_identity(alignment_list: List[str]) -> float:
    """Compute the identity rate between two sequences

    :param alignment_list:  (list) A list of aligned sequences in the format ["SE-QUENCE1", "SE-QUENCE2"]
    :return: (float) The rate of identity between the two sequences.
    """
    sequence1 = alignment_list[0]
    sequence2 = alignment_list[1]
    nb_identity, seq_align_len = 0, len(sequence1)
    for i in range(seq_align_len):
        if sequence1[i] == sequence2[i]:
            nb_identity += 1
    return (nb_identity/seq_align_len)*100

def abundance_greedy_clustering(amplicon_file: Path, minseqlen: int, mincount: int, chunk_size: int, kmer_size: int) -> List:
    """Compute an abundance greedy clustering regarding sequence count and identity.
    Identify OTU sequences.

    :param amplicon_file: (Path) Path to the amplicon file in FASTA.gz format.
    :param minseqlen: (int) Minimum amplicon sequence length.
    :param mincount: (int) Minimum amplicon count.
    :param chunk_size: (int) A fournir mais non utilise cette annee
    :param kmer_size: (int) A fournir mais non utilise cette annee
    :return: (list) A list of all the [OTU (str), count (int)] .
    """
    seq_with_occurrence = list(dereplication_fulllength(amplicon_file, minseqlen, mincount))
    seq_OTU = [[seq_with_occurrence[0][0], seq_with_occurrence[0][1]]]
    for i in range(1, len(seq_with_occurrence)):
        count_identity_seq = 0
        for sequences in seq_OTU:
            seq1, seq2, count = seq_with_occurrence[i][0], sequences[0], seq_with_occurrence[i][1]
            seq1_align, seq2_align = nw.global_align(seq1, seq2, gap_open=-1, gap_extend=-1, matrix=str(Path(__file__).parent / "MATCH"))
            print(seq1_align, seq2_align)
            if get_identity([seq1_align, seq2_align]) >= 97:
                count_identity_seq += 1
        if count_identity_seq == 0:
            seq_OTU.append([seq1, count])
    return seq_OTU

def write_OTU(OTU_list: List, output_file: Path) -> None:
    """Write the OTU sequence in fasta format.

    :param OTU_list: (list) A list of OTU sequences
    :param output_file: (Path) Path to the output file
    """
    with open(output_file, "w") as file_write:
        counter = 1
        for element in OTU_list:
            file_write.write(f">OTU_{counter} occurrence:{element[1]}\n")
            for i in range(0, len(element[0]), 80):
                file_write.write(f"{element[0][i:i+80]}\n")
            counter += 1


#==============================================================
# Main program
#==============================================================
def main(): # pragma: no cover
    """
    Main program function
    """
    # Get arguments
    args = get_arguments()
    # Votre programme ici



if __name__ == '__main__':
    main()
