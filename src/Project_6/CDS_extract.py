import argparse
import collections
from datetime import date
import gzip
import sys
import os
import pathlib

###################
gencode = {
      'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
      'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
      'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
      'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
      'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
      'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
      'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
      'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
      'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
      'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
      'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
      'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
      'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
      'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
      'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
      'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W'}

def translate_frame(sequence):
    translate = ''.join([gencode.get(sequence[3 * i:3 * i + 3], 'X') for i in range(len(sequence) // 3)])
    return translate


def reverseCorrectLoci(seq_length,first,second,third): # here for the negative loci correction

    if second == None:
        corrected_start = max(seq_length - int(third),1)
        corrected_stop = max(seq_length - int(first-1),1)
        return corrected_start, corrected_stop
    else: # Needs to be checked
        corrected_start = max(seq_length - int(third),1)
        corrected_mid = max(seq_length - int(second-3),1)
        corrected_stop = max(seq_length - int(first),1)
        return corrected_start, corrected_mid, corrected_stop
#################################
def revCompIterative(watson): #Gets Reverse Complement

    complements = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                   'R': 'Y', 'Y': 'R', 'S': 'S', 'W': 'W', 'K': 'M',
                   'M': 'K', 'V': 'B', 'B': 'V', 'H': 'D', 'D': 'H'}
    watson = watson.upper()
    watsonrev = watson[::-1]
    crick = ""
    for nt in watsonrev: # make dict to catch bad nts - if more than 1 output them to std error
        try:
            crick += complements[nt]
        except KeyError:
            crick += nt
            #ns_nt[nt] +=1
    return crick


def write_fasta(dna_regions, fasta_outfile):
    for dna_region, dna_region_ur in dna_regions.items():
      #  fasta_outfile.write('##sequence-region\t' + dna_region + ' 1 ' + str(len(dna_region_ur[0])) + '\n')
        if dna_region_ur[3]:
            for storf, seq in dna_region_ur[3].items():
                aa_seq = translate_frame(seq)
                fasta_outfile.write('>' + storf + '\n' + aa_seq + '\n')
    fasta_outfile.close()

def write_gff(dna_regions,options,gff_outfile, gff):
    gff_outfile.write("##gff-version\t3\n#\tStORF-Extractor \n#\tRun Date:" + str(date.today()) + '\n')
    gff_outfile.write('##StORF-Reporter ' + '\n')
    for seq_reg in dna_regions:
        gff_outfile.write('##sequence-region\t' + seq_reg + ' 1 ' + str(dna_regions[seq_reg][1]) + '\n')
    gff_outfile.write("##Original File: " + gff.split('/')[-1] + '\n\n')
    for dna_region, orf_data in dna_regions.items():
        #ur_ident = dna_region + options.ident
        if orf_data[3]:
            for orf, seq in orf_data[3].items():

                entry = (dna_region + '\tStORF-Reporter\tCDS\t' + orf.split('_')[0] + '\t' + orf.split('_')[1]+ '\t.\t' +
                         orf.split('_')[2].split(';')[0]  + '\t.\t' + orf.split(';',1)[1] + '\n')
                gff_outfile.write(entry)
    gff_outfile.close()

def fasta_load(fasta_in):
    dna_regions = collections.OrderedDict()
    first = True
    if '>' in fasta_in.readline().rstrip():
        fasta_in.seek(0)
        #### Default for when presented with standard fasta file
        for line in fasta_in:
            line = line.strip()
            if line.startswith('>') and first == False:  # Check if first seq in file
                dna_region_length = len(seq)
                dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            elif line.startswith('>'):
                seq = ''
                dna_region_id = line.split()[0].replace('>', '')
            else:
                seq += str(line)
                first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
    elif '##' in fasta_in.readline().rstrip(): # Clunky and may fall over
        fasta_in.seek(0)
        #### Called when presented with Prokka GFF file so must get fasta from inside it
        ### Get to genome seq
        at_FASTA = False
        for line in fasta_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
            if line.startswith('##FASTA'):  # Not to crash on empty lines in GFF
                at_FASTA = True
            elif at_FASTA == True:
                line = line.strip()
                if line.startswith('>') and first == False:  # Check if first seq in file
                    dna_region_length = len(seq)
                    dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                elif line.startswith('>'):
                    seq = ''
                    dna_region_id = line.split()[0].replace('>', '')
                else:
                    seq += str(line)
                    first = False
        dna_region_length = len(seq)
        dna_regions.update({dna_region_id: (seq, dna_region_length, list(), None)})

    return dna_regions


def gff_load(options,gff_in,dna_regions):
    has_storfs = False
    for line in gff_in:  # Get gene loci from GFF - ID=Gene will also classify Pseudogenes as genes
        line_data = line.split()
        if line.startswith('\n') or line.startswith('#'):  # Not to crash on empty lines in GFF
            continue
        else:
            try:
                if line_data[0] in dna_regions:
                    if 'CDS' in line_data[2]:
                        annotations = line_data[1]
                        has_storfs = True
                        pos = line_data[3] + '_' + line_data[4] + '_' + line_data[6] + '_' + annotations
                        if pos not in dna_regions[line_data[0]][2]:
                            dna_regions[line_data[0]][2].append(pos) # This will add to list

            except IndexError:
                continue
    if has_storfs == False:
        dna_regions = None
    return dna_regions




def fasta_extractor(options):

    fasta_in = open(options.fasta,'r',encoding='unicode_escape')
    dna_regions = fasta_load(fasta_in)
    gff_in = open(options.gff,'r',encoding='unicode_escape')
    dna_regions = gff_load(options,gff_in,dna_regions)

    if dna_regions == None and options.verbose == True:
        print("No seqs to extract from " + options.gff)
        return


    for (key,(seq,seq_length,posns,tmp))  in dna_regions.items(): #Extract URs from 1 dna_region at a time
        Extracted_ORFs = collections.OrderedDict()
        seq_rev = revCompIterative(seq)
        if posns: # If UR has a pos
            for pos in posns: # Iterate over GFF loci and measure flanking regions for potential URs
                start = int(pos.split('_')[0])
                stop = int(pos.split('_')[1])
                frame = pos.split('_')[2].split(';')[0]
                ###### This hack is to get over GFF errors where genome-long annotations
                if stop-start >= 100000:
                    if options.verbose == True:
                        print("UR " + pos + " is more than 100,000 kbs - Please Check Annotation")
                    continue
                if frame == '+':
                    ORF_seq = seq[start-1:stop]
                elif frame == '-':
                    rev_corrected_start, rev_corrected_stop = reverseCorrectLoci(seq_length, start, None, stop)
                    ORF_seq = seq_rev[rev_corrected_start:rev_corrected_stop]

                Extracted_ORFs[pos] = ORF_seq
        dna_regions.update({key: (seq, seq_length, posns, Extracted_ORFs)})



    write_fasta(dna_regions, options.fasta_outfile)

def main():

    parser = argparse.ArgumentParser(description='None for now.')
    parser._action_groups.pop()

    required = parser.add_argument_group('Required Arguments')
    required.add_argument('-fasta', action='store', dest='fasta', required=False,

                        help='fasta?\n')
    required.add_argument('-gff', action='store', dest='gff', default='', required=False,
                        help='gff')




    options = parser.parse_args()
    filename_tmp = options.gff.replace('.gff','_Extracted_Fasta.fasta')
    options.fasta_outfile = open(filename_tmp,'w')
    fasta_extractor(options)



if __name__ == "__main__":
    main()
    print("Complete")
