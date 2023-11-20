#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Mar 29 12:14:25 2023

@author: tomslab
"""

import os
from argparse import ArgumentParser
import csv

def escritura(lista, nombre):

    archivo = open(nombre, 'w')
    archivo.close()    
    
    for linea in lista:
        
        with open(nombre, mode='a') as result_file:
            line_writer = csv.writer(result_file, delimiter='\t', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
            line_writer.writerow(linea)      
               
#enddef

def escritura_comma_separated(lista, nombre):

    archivo = open(nombre, 'w')
    archivo.close()    
    
    for linea in lista:
        
        with open(nombre, mode='a') as result_file:
            line_writer = csv.writer(result_file, delimiter=',', quotechar='"', quoting=csv.QUOTE_MINIMAL)
        
            line_writer.writerow(linea)      
               
#enddef

def liftoff(reference_gff, output_folder, types_file, target_genome, reference_genome):
    
    print('\nRunning liftoff\n')
    
    liftoff_output = output_folder + '/liftoff_output.gff3'
    
    liftoff_command = 'liftoff -g ' + reference_gff + ' -o ' + liftoff_output + ' -f ' + types_file + ' ' + target_genome + ' ' + reference_genome + ' -copies'
    os.system(liftoff_command)
    
# enddef

def mybannotator(mybannotator_folder, mybannotator_out_folder, target_genome_proteins):
    
    print('\nRunning MYBannotator\n')
    
    baits = mybannotator_folder + '/MYB_baits.fasta'
    info = mybannotator_folder + '/MYB_baits.txt'
    hmm = mybannotator_folder + '/MYB_baits.hmm'
    program = mybannotator_folder + '/MYB_annotator.py'
    
    command = 'python3 ' + program + ' --baits ' + baits + ' --info ' + info + ' --out ' + mybannotator_out_folder + ' --subject ' + target_genome_proteins + ' --hmm ' + hmm
    print(command)
    os.system(command)

# enddef

def extracting_conversions(output_folder, target_genome_anno):
    
    print('\nGenerating conversions between reference and target genes\n')
    
    '''Reading the liftoff over the haplotig haplotype'''

    haplotig_cromo = []
    haplotig_starts = []
    haplotig_ends = []
    haplotig_ids = []
    
    liftoff_gff_output = output_folder + '/liftoff_output.gff3'
    file = open(liftoff_gff_output)
    for line in file:
        line = line.strip().split('\t')
        try:
            if line[2] == 'gene' or line[2] == 'pseudogene':
                
                if line[0] not in haplotig_cromo:
                    haplotig_cromo.append(line[0])
                    haplotig_starts.append([int(line[3])])
                    haplotig_ends.append([int(line[4])])
                    haplotig_ids.append([line[8].split(';')[0].split('=')[1]])
                    
                else:
                    index = haplotig_cromo.index(line[0])
                    haplotig_starts[index].append(int(line[3]))
                    haplotig_ends[index].append(int(line[4]))
                    haplotig_ids[index].append(line[8].split(';')[0].split('=')[1])
                    
        except:
            continue
            
    '''Reading the target genome gff file'''

    semillon_cromo = []
    semillon_starts = []
    semillon_ends = []
    semillon_ids = []

    file = open(target_genome_anno)
    for line in file:
        line = line.strip().split('\t')
        try:
            if line[2] == 'gene' or line[2] == 'pseudogene':
                if line[0] not in semillon_cromo:
                    semillon_cromo.append(line[0])
                    semillon_starts.append([int(line[3])])
                    semillon_ends.append([int(line[4])])
                    semillon_ids.append([line[8].split(';')[0].split('=')[1]])
                    
                else:
                    index = semillon_cromo.index(line[0])
                    semillon_starts[index].append(int(line[3]))
                    semillon_ends[index].append(int(line[4]))
                    semillon_ids[index].append(line[8].split(';')[0].split('=')[1])
        except:
            continue

    '''Looking for overlapping genes between the 2 annotations'''

    result = []
    not_found_genes = []
    added_genes = []
    repeated_genes = []

    for chromo_counter in range(len(haplotig_cromo)):
        try:
            index_semillon = semillon_cromo.index(haplotig_cromo[chromo_counter])
            
            for pn in range(len(haplotig_ids[chromo_counter])):
                
                found = False
                pn_id = haplotig_ids[chromo_counter][pn]
                pn_start = haplotig_starts[chromo_counter][pn]
                pn_end = haplotig_ends[chromo_counter][pn]
                
                for sm in range(len(semillon_ids[index_semillon])):
                    
                    sm_id = semillon_ids[index_semillon][sm]
                    sm_start = semillon_starts[index_semillon][sm]
                    sm_end = semillon_ends[index_semillon][sm]
                    
                    if ((pn_start < sm_end) and (pn_end >= sm_end)) or ((pn_start <= sm_start) and (pn_end > sm_start)) or ((pn_start <= sm_start) and (pn_end >= sm_end)) or ((pn_start >= sm_start) and (pn_end <= sm_end)):
                        
                        found = True
                        result.append([pn_id, sm_id])
                        if pn_id in added_genes:
                            repeated_genes.append(pn_id)
                        else:
                            added_genes.append(pn_id)   
                        
                if found == False:
                    not_found_genes.append(pn_id)
            
        except:
            for entry in haplotig_ids[chromo_counter]:
                not_found_genes.append(entry)
                
    total_genes_lifted = 0
    for entry in haplotig_ids:
        for x in entry:
            total_genes_lifted +=1
            
    print('\nOut of ', total_genes_lifted, ' lifted genes, ', len(added_genes) - len(repeated_genes), ' overlap with a unique gene in the target genome, ', len(repeated_genes), ' overlap with more than one gene and ', len(not_found_genes), ' do not overlap.')

    result = [['Reference ID', 'Target ID']] + result
    output = output_folder + '/conversions_between_ref_and_targets.csv'
    escritura(result, output)

# enddef

def rbh_blast(output_folder, target_genome_cds, reference_genome_cds):
    
    output = output_folder + '/RBH_blast/output.tsv'
    os.system('mkdir -p ' + output_folder + '/RBH_blast')
    os.system('crb-blast -q ' + reference_genome_cds + ' -t ' + target_genome_cds + ' --threads 8 -o '+ output)
    
#enddef

def main(reference_gff, output_folder, types_file, target_genome, reference_genome, mybannotator_folder, mybannotator_out_folder, target_genome_proteins, target_genome_anno, myb_list_file, target_genome_cds, reference_genome_cds):
    
    liftoff(reference_gff, output_folder, types_file, target_genome, reference_genome)
    
    mybannotator(mybannotator_folder, mybannotator_out_folder, target_genome_proteins)
    
    extracting_conversions(output_folder, target_genome_anno)
    
    rbh_blast(output_folder, target_genome_cds, reference_genome_cds)
    
# enddef


'''MAIN PROGRAM'''


if __name__ == '__main__':
    
    parser = ArgumentParser()
    
    parser.add_argument(
        '--target_genome',
        dest='target_genome',
        action='store',
        required=True,
        help='Path of the target genome.'
        )
    
    parser.add_argument(
        '--target_genome_anno',
        dest='target_genome_anno',
        action='store',
        required=True,
        help='Path of the gff annotations file of the target genome.'
        )
    
    parser.add_argument(
        '--target_genome_proteins',
        dest='target_genome_proteins',
        action='store',
        required=True,
        help='Path of the fasta file that contains the sequence of all the proteins described in the target annotations file.'
        )
    
    parser.add_argument(
        '--target_genome_cds',
        dest='target_genome_cds',
        action='store',
        required=True,
        help='Path of the fasta file that contains the CDS sequence of all the genes described in the target annotations file.'
        )
    
    parser.add_argument(
        '--reference_genome',
        dest='reference_genome',
        action='store',
        required=True,
        help='Path of the reference genome.'
        )
    
    parser.add_argument(
        '--reference_gff',
        dest='reference_gff',
        action='store',
        required=True,
        help='Path of the types file required by liftoff.'
        )

    parser.add_argument(
        '--reference_genome_proteins',
        dest='reference_genome_proteins',
        action='store',
        required=True,
        help='Path of the fasta file that contains the sequence of all the proteins described in the reference annotations file.'
        )

    parser.add_argument(
        '--reference_genome_cds',
        dest='reference_genome_cds',
        action='store',
        required=True,
        help='Path of the fasta file that contains the CDS sequence of all the genes described in the reference annotations file.'
        )
    
    parser.add_argument(
        '--output_folder',
        dest='output_folder',
        action='store',
        required=True,
        help='Path of the folder where outputs of this pipeline will be storaged.'
        )
    
    parser.add_argument(
        '--types_file',
        dest='types_file',
        action='store',
        required=True,
        help='Path of the types file required by liftoff.'
        )
    
    parser.add_argument(
        '--mybannotator_folder',
        dest='mybannotator_folder',
        action='store',
        required=True,
        help='Path of the folder where mybannotator tool is installed.'
        )
    
    parser.add_argument(
        '--myb_list_file',
        dest='myb_list_file',
        action='store',
        required=True,
        help='Path of the tsv file that contains in the first columnn the reference MYB gene IDs and in the second column their gene names.'
        )
    
    args = parser.parse_args ()
    
    target_genome = args.target_genome
    target_genome_anno = args.target_genome_anno
    target_genome_proteins = args.target_genome_proteins
    target_genome_cds = args.target_genome_cds
    reference_genome = args.reference_genome
    reference_gff = args.reference_gff
    reference_genome_proteins = args.reference_genome_proteins
    reference_genome_cds = args.reference_genome_cds
    output_folder = args.output_folder
    types_file = args.types_file
    mybannotator_folder = args.mybannotator_folder
    myb_list_file = args.myb_list_file
    
    mybannotator_out_folder = output_folder + '/mybanno_out'
    
    main(reference_gff, output_folder, types_file, target_genome, reference_genome, mybannotator_folder, mybannotator_out_folder, target_genome_proteins, target_genome_anno, myb_list_file, target_genome_cds, reference_genome_cds)