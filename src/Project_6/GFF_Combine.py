from importlib import import_module
import argparse
import collections
from datetime import date
import sys
try:
    from ORForise.utils import sortORFs  # Calling from ORForise via pip
    from .Constants import *



########################################


def main():
    print("Thank you for using ORForise\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='ORForise ' + ORForise_Version + ': GFF-Adder Run Parameters.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-dna', dest='genome_DNA', required=True, help='Genome DNA file (.fa) which both annotations '
                                                                    'are based on')
    required.add_argument('-ref', dest='reference_annotation', required=True,
                        help='Which reference annotation file to use as reference?')
    required.add_argument('-at', dest='additional_tool', required=True,
                        help='Which format to use for additional annotation?')
    required.add_argument('-add', dest='additional_annotation', required=True,
                        help='Which annotation file to add to reference annotation?')
    required.add_argument('-o', dest='output_file', required=True,
                        help='Output filename')

    optional = parser.add_argument_group('Optional Arguments')
    optional.add_argument('-rt', dest='reference_tool', required=False,
                        help='Which tool format to use as reference? - If not provided, will default to the '
                             'standard GFF format and will only look for "CDS" features')
    optional.add_argument('-gene_ident', action='store', dest='gene_ident', default='CDS',
                        help='Identifier used for identifying genomic features in reference annotation "CDS,rRNA,tRNA"')
    optional.add_argument('-mc', dest='mark_consensus', default=False, type=bool, required=False,
                        help='Default - False: Mark reference annotations which where present in the additional tool annotation')
    optional.add_argument('-olap', dest='overlap', default=50, type=int, required=False,
                        help='Maximum overlap between reference and additional genic regions (CDS,rRNA etc) - Default: 50 nt')

    options = parser.parse_args()

    gff_adder(options)


if __name__ == "__main__":
    main()
    print("Complete")
