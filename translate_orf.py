#! /usr/bin/env python3

import sys
import re
import find_orf
import translate

#Import package

def main():
    import argparse

# Create a command-line parser object
    parser = argparse.ArgumentParser()

    start = ['AUG']
    stop = ['UGA', 'UAA', 'UAG']
    genetic_code = {'GUC': 'V', 'ACC': 'T', 'GUA': 'V', 'GUG': 'V', 'ACU': 'T', 'AAC': 'N', 'CCU': 'P', 'UGG': 'W', 'AGC': 'S', 'AUC': 'I', 'CAU': 'H', 'AAU': 'N', 'AGU': 'S', 'GUU': 'V', 'CAC': 'H', 'ACG': 'T', 'CCG': 'P', 'CCA': 'P', 'ACA': 'T', 'CCC': 'P', 'UGU': 'C', 'GGU': 'G', 'UCU': 'S', 'GCG': 'A', 'UGC': 'C', 'CAG': 'Q', 'GAU': 'D', 'UAU': 'Y', 'CGG': 'R', 'UCG': 'S', 'AGG': 'R', 'GGG': 'G', 'UCC': 'S', 'UCA': 'S', 'UAA': '*', 'GGA': 'G', 'UAC': 'Y', 'GAC': 'D', 'UAG': '*', 'AUA': 'I', 'GCA': 'A', 'CUU': 'L', 'GGC': 'G', 'AUG': 'M', 'CUG': 'L', 'GAG': 'E', 'CUC': 'L', 'AGA': 'R', 'CUA': 'L', 'GCC': 'A', 'AAA': 'K', 'AAG': 'K', 'CAA': 'Q', 'UUU': 'F', 'CGU': 'R', 'CGC': 'R', 'CGA': 'R', 'GCU': 'A', 'GAA': 'E', 'AUU': 'I', 'UUG': 'L', 'UUA': 'L', 'UGA': '*', 'UUC': 'F'}
    
    # Tell the parser what command-line arguments this script can receive
    parser.add_argument('sequence',
            metavar = 'SEQUENCE',
            type = str,
            help = ('The sequence to search for an open-reading frame. '
                    'If the path flag (\'-p\'/\'--path\') is specified, '
                    'then this should be a path to a file containing the '
                    'sequence to be searched.'))
    parser.add_argument('-p', '--path',
            action = 'store_true',
            help = ('The sequence argument is a path to a file'
                    'containing the sequence to be searched.'))
    parser.add_argument('-s', '--start-codon',
            type = str,
            action = 'append', # append each argument to a list
            default = None,
            help = ('One or more start codon(s).'
                   'Default: {0}.'.format(" ".join(start))))
    parser.add_argument('-x', '--stop-codon',
            type = str,
            action = 'append', # one or more arguments
            default = None,
            help = ('One or more stop codon(s).'
                   'Default: {0}.'.format(" ".join(stop))))

    # Parse the command-line arguments into a 'dict'-like container
    args = parser.parse_args()

    # Check to see if the path option was set to True by the caller. If so, parse
    # the sequence from the path
    if args.path:
        sequence = parse_sequence_from_path(args.sequence)
    else:
        sequence = args.sequence

    #Check to see if start/stop codons were provided by caller. If not, use the defaults
    if not args.start_codon:
        args.start_codon = start
    if not args.stop_codon:
        args.stop_codon = stop

    orf = find_orf.find_first_orf(sequence = sequence,
            start_codons = args.start_codon,
            stop_codons = args.stop_codon)


    translation = translate.translate_sequence(rna_sequence = orf, genetic_code = genetic_code)
    sys.stdout.write('{}\n'.format(translation))


if __name__ == '__main__':
    main()
