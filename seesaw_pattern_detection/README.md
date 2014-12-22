Scripts used to perform the seesaw pattern detection algorithm from the RNA-Seq and FUS iCLIP data

# File seesaw_step1.R  

Input files:
* data/Final_150K_junction_analysis.csv which contains the information about the split read junctions identified by the split read analysis
* data/transcript_info.tab information about each transcript

Output file:
* data/ensembl_annotations_complete.csv which contains the augmented annotation data, combining transcript information and split read data



File seesaw_step2.R