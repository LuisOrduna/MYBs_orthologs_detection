# MYBs_orthologs_detection
Pipeline for R2R3-MYB orthologs detection

The pipeline takes as input the following files:

target_genome='/path/to/target/cultivar/reference/fastafile'
target_genome_anno='/path/to/target/cultivar/annotations/gff3'
target_genome_proteins='/path/to/target/cultivar/proteins/fastafile'
target_genome_cds='/path/to/target/cultivar/CDSs/fastafile'

reference_genome='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN40X/PN40024_40X_REF_chloro_mito.chr_renamed.fasta'
reference_gff='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN40X/original_PN40024_pseudomolecules.v4.3.BETA.gff3'
reference_genome_proteins='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN40X/original_PN40024.4.3_proteins_modified.fasta'
reference_genome_cds='/media/tomslab/hard_drive/Luis/anotaciones/vitis/PN40X/original_PN40024.4.3_CDSs.fasta'

output_folder='/media/tomslab/hard_drive/Luis/paper_MYBs/Semillon_haplo1_results/'
types_file='/media/tomslab/hard_drive/Luis/anotaciones/vitis/liftoff_PN40X_4.3_on_Semillon/types.txt'
mybannotator_folder='/home/tomslab/MYB_annotator'
myb_list_file='/media/tomslab/hard_drive/Luis/paper_MYBs/pipeline/Semillon_haplo1/mybs_file.tsv'
