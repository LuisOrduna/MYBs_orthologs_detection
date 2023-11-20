#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Oct 31 16:24:32 2023

@author: tomslab
"""

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

def integration(myb_list_file, output_folder):

    mybs_reference_IDs = []
    mybs_reference_symbols = []
    file = open(myb_list_file)
    for line in file:
        line = line.strip().split('\t')
        mybs_reference_IDs.append(line[0])
        mybs_reference_symbols.append(line[1])
    
    mybannotator_out = []
    analyzed = []
    mybannotator_out_file = output_folder + '/mybanno_out/RESULTS/04a_MYB_domain_check.txt'
    file = open(mybannotator_out_file)
    for line in file:
        line = line.strip().split('\t')
        if line[2] == 'R2R3':
            identifier = line[0].split('.')
            removed = identifier.pop(-1)
            identifier = '.'.join(identifier)
            if identifier not in analyzed:
                analyzed.append(identifier)
                mybannotator_out.append([identifier, 'R2R3'])
                
    rbh_blast = []
    path = output_folder + '/RBH_blast/output.tsv'
    file = open(path)
    for line in file:
        line = line.strip().split('\t')
        ref_id = line[0].split('.')
        removed = ref_id.pop(-1)
        ref_id = '.'.join(ref_id)
        target_id = line[1].split('.')
        removed = target_id.pop(-1)
        target_id = '.'.join(target_id)
        if [ref_id, target_id] not in rbh_blast:
            rbh_blast.append([ref_id, target_id])
        
    liftoff = []
    path = output_folder + '/conversions_between_ref_and_targets.csv'
    file = open(path)
    for line in file:
        line = line.strip().split('\t')
        liftoff.append(line)
    removed = liftoff.pop(-1)    
    
    orthofinder = []
    path = output_folder + '/orthofinderSequences/results/Results_orthofinder_out/Orthologues/reference_genome.tsv'
    file = open(path)
    for line in file:
        line = line.strip().split('\t')
        orthofinder.append([line[2], line[3]])
    removed = orthofinder.pop(0)
    
    common_mybs = []
    
    exclusive_mybs = []
    
    for entry in mybannotator_out:
        liftoff_conversion = '-'
        rbh_blast_conversion = '-'
        orthofinder_conversion = '-'
        mybanno_classification = entry[1]
        liftoff_symbol = ''
        rbh_blast_symbol = ''
        orthofinder_symbol = ''
        target_id = entry[0]
        
        '''Adding liftoff reference IDs and gene symbols'''
        
        for conversion in liftoff:
            if target_id in conversion[1]:
                
                if liftoff_conversion != '-':
                    liftoff_conversion = liftoff_conversion + '; ' + conversion[0]
                else:
                    liftoff_conversion = conversion[0]
        
        temporal = liftoff_conversion.split('; ')
        for gene in temporal:
            if gene in mybs_reference_IDs:
                if liftoff_symbol == '':
                    liftoff_symbol = mybs_reference_symbols[mybs_reference_IDs.index(gene)]
                else:
                    liftoff_symbol = liftoff_symbol + '; ' + mybs_reference_symbols[mybs_reference_IDs.index(gene)]
            else:
                if liftoff_symbol == '':
                    liftoff_symbol = '-'
                else:
                    liftoff_symbol = liftoff_symbol + '; -'
        
        '''Adding RBH reference IDs and gene symbols'''
            
        for conversion in rbh_blast:
            
            if target_id in conversion[1]:
                
                if rbh_blast_conversion != '-':
                    rbh_blast_conversion = rbh_blast_conversion + '; ' + conversion[0]
                else:
                    rbh_blast_conversion = conversion[0]
        
        temporal = rbh_blast_conversion.split('; ')
        for gene in temporal:
            if gene in mybs_reference_IDs:
                if rbh_blast_symbol == '':
                    rbh_blast_symbol = mybs_reference_symbols[mybs_reference_IDs.index(gene)]
                else:
                    rbh_blast_symbol = rbh_blast_symbol + '; ' + mybs_reference_symbols[mybs_reference_IDs.index(gene)]
            else:
                if rbh_blast_symbol == '':
                    rbh_blast_symbol = '-'
                else:
                    rbh_blast_symbol = rbh_blast_symbol + '; -'
                    
        '''Adding Orthofinder reference IDs and gene symbols'''
        
        for conversion in orthofinder:
            
            if target_id in conversion[1]:
                
                genes = conversion[0].split(', ')
                unique_genes = []
                for gene in genes:
                    gene = gene.split('.')
                    removed = gene.pop(-1)
                    gene = '.'.join(gene)
                    if gene not in unique_genes:
                        unique_genes.append(gene)
                        
                for gene in unique_genes:
                    
                    if orthofinder_conversion == '-':
                        orthofinder_conversion = gene
                    else:
                        if gene not in orthofinder_conversion:
                            orthofinder_conversion = orthofinder_conversion + '; ' + gene
                            
        temporal = orthofinder_conversion.split('; ')
        for gene in temporal:
            if gene in mybs_reference_IDs:
                if orthofinder_symbol == '':
                    orthofinder_symbol = mybs_reference_symbols[mybs_reference_IDs.index(gene)]
                else:
                    orthofinder_symbol = orthofinder_symbol + '; ' + mybs_reference_symbols[mybs_reference_IDs.index(gene)]
            else:
                if orthofinder_symbol == '':
                    orthofinder_symbol = '-'
                else:
                    orthofinder_symbol = orthofinder_symbol + '; -'
            
        
        '''Integration of all methods'''
        
        temp_orthofinder = orthofinder_conversion.split('; ')
        temp_rbh_blast = rbh_blast_conversion.split('; ')
        temp_liftoff = liftoff_conversion.split('; ')
        
        if (temp_liftoff[0] != '-') and (len(temp_liftoff) == 1) and (temp_liftoff == temp_rbh_blast) and (temp_liftoff == temp_orthofinder):
            
            unique_detection = 'YES'
        
        else:
            
            unique_detection = 'NO'
        
        '''Splitting between common and exclusive MYBs'''
        
        '''Shared MYBs are defined as CS IDs that have at least one conversion to 
        a R2R3-MYB by at least one method. If it has a conversion but to a gene
        that is not a R2R3-MYB in the ref genome, the CS ID is considered as 
        exclusive from the target genome, since its ortholog is not a R2R3.-MYB
        in the ref genome'''
        
        temp_ortho_symbols = list(set(orthofinder_symbol.split('; ')))
        temp_rbh_symbols = list(set(rbh_blast_symbol.split('; ')))
        temp_liftoff_symbols = list(set(liftoff_symbol.split('; ')))
        
        if (temp_liftoff_symbols == ['-']) and (temp_rbh_symbols == ['-']) and (temp_ortho_symbols == ['-']):
            exclusive_mybs.append([target_id, mybanno_classification, unique_detection, liftoff_conversion, liftoff_symbol, rbh_blast_conversion, rbh_blast_symbol, orthofinder_conversion, orthofinder_symbol])
        else:
            common_mybs.append([target_id, mybanno_classification, unique_detection, liftoff_conversion, liftoff_symbol, rbh_blast_conversion, rbh_blast_symbol, orthofinder_conversion, orthofinder_symbol])
                
    # output = output_folder + '/shared_R2R3_MYBs.tsv'
    # escritura(common_mybs, output)
    # output = output_folder + '/exclusive_R2R3_MYBs.tsv'
    # escritura(exclusive_mybs, output)
    
    result_def = [['Target genome ID',
                    'MYB classification according to MYBannotator',
                    'R2R3-MYB ortholog classification',
                    'Reference gene ID according to liftoff',
                    'Reference gene symbol according to liftoff',
                    'Reference gene ID according to CRBH blast',
                    'Reference gene symbol according to CRBH blast',
                    'Reference gene ID according to Orthofinder',
                    'Reference gene symbol according to Orthofinder']]
    
    for entry in common_mybs:
        
        if entry[2] == 'YES':
            entry[2] = 'Unanimous R2R3-MYB ortholog'
        if entry[2] == 'NO':
            entry[2] = 'Not unanimous R2R3-MYB ortholog'
        result_def.append(entry)
        
    for entry in exclusive_mybs:
        entry[2] = 'Target genome exclusive R2R3-MYB'
        result_def.append(entry)
    
    output = output_folder + '/pipeline_final_classification.tsv'
    escritura(result_def, output)
    
# enddef

'''MAIN PROGRAM'''

if __name__ == '__main__':
    
    parser = ArgumentParser()
    
    parser.add_argument(
        '--output_folder',
        dest='output_folder',
        action='store',
        required=True,
        help='Path of the folder where outputs of this pipeline will be storaged.'
        )
    
    parser.add_argument(
        '--myb_list_file',
        dest='myb_list_file',
        action='store',
        required=True,
        help='Path of the tsv file that contains in the first columnn the reference MYB gene IDs and in the second column their gene names.'
        )
    
    args = parser.parse_args ()
    
    output_folder = args.output_folder
    myb_list_file = args.myb_list_file
    
    integration(myb_list_file, output_folder)
                    
    
        
    
        
                
