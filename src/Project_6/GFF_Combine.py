from importlib import import_module
import argparse
import collections
from datetime import date
import sys

try:
    from ORForise.utils import sortORFs  # Calling from ORForise via pip
    from ORForise.GFF_Adder import gff_adder  # Calling from ORForise via pip
    from .Constants import *
except (ModuleNotFoundError, ImportError, NameError, TypeError) as error:
    sys.path.insert(0, '../../../ORForise/src/ORForise/') # Calling from ORForise locally (StORF_Reporter and ORForise in same dir)
    from utils import sortORFs
    from GFF_Adder import gff_adder
    from Constants import *


########################################


def main():
    print("Thank you for using GF_Finisher\nPlease report any issues to: https://github.com/NickJD/ORForise/issues\n#####")

    parser = argparse.ArgumentParser(description='GFF Combiner: Run Parameters.')
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
    optional.add_argument('-pred', dest='deepgoplus', required=False,
                        help='Use DeepGoPlus to predict unnanotated protein functions')
    
    options = parser.parse_args()

    gff_adder(options)


if __name__ == "__main__":
    main()
    print("Complete")
