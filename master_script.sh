#!/bin/bash

shopt -s expand_aliases
source ~/.bashrc

# Defining paths for the execution

target_genome='/media/tomslab/hard_drive/Luis/anotaciones/vitis/CS/haplo1/VITVvi_vCabSauv08_v1.1.pseudomolecules.hap1.fasta'
target_genome_anno='/media/tomslab/hard_drive/Luis/anotaciones/vitis/CS/haplo1/VITVvi_vCabSauv08_v1.1_haplo1.gff3'
target_genome_proteins='/media/tomslab/hard_drive/Luis/anotaciones/vitis/CS/haplo1/VITVvi_vCabSauv08_v1.1.pseudomolecules.haplo1.protein.fasta'
target_genome_cds='/media/tomslab/hard_drive/Luis/anotaciones/vitis/CS/haplo1/VITVvi_vCabSauv08_v1.1.pseudomolecules.haplo1.CDS.fasta'

reference_genome='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN/VCost_genome.fasta'
reference_gff='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN/VCost.v3_27.gff3'
reference_genome_proteins='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN/vitviv2.pep.fasta'
reference_genome_cds='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN/vitviv2.cds.fasta'

output_folder='/media/tomslab/hard_drive/Luis/paper_MYBs/CS_haplo1_results/'
types_file='/media/tomslab/hard_drive/Luis/anotaciones/vitis/liftoff_PN40X_4.3_on_Semillon/types.txt'
mybannotator_folder='/home/tomslab/MYB_annotator'
myb_list_file='/media/tomslab/hard_drive/Luis/paper_MYBs/pipeline/CS_haplo1/mybs_file.tsv'

# Running the execution

mkdir -p $output_folder

python3 master_script.py --target_genome $target_genome --target_genome_anno $target_genome_anno --target_genome_proteins $target_genome_proteins --target_genome_cds $target_genome_cds --reference_genome $reference_genome --reference_gff $reference_gff --reference_genome_proteins $reference_genome_proteins --reference_genome_cds $reference_genome_cds --output_folder $output_folder --types_file $types_file --mybannotator_folder $mybannotator_folder --myb_list_file $myb_list_file

mkdir $output_folder/orthofinderSequences
cp $reference_genome_proteins $output_folder/orthofinderSequences/reference_genome.fa
cp $target_genome_proteins $output_folder/orthofinderSequences/target_genome.fa

orthofinder --fewer-files -f $output_folder/orthofinderSequences -o $output_folder/orthofinderSequences/results -n orthofinder_out

python3 pipeline_integration.py --output_folder $output_folder --myb_list_file $myb_list_file
