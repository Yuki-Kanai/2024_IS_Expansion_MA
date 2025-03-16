## Directory Structure and Descriptions

This directory contains several folders related to the analysis and classification of Insertion Sequences (IS) events, their positions, and visualization of results. Each folder serves a specific purpose in the analysis pipeline.

### `classify_IS_events` (Main Datasets)
Contains comparisons of two consecutive rounds of sequencing to identify which Insertion Sequences (IS) correspond to which, by comparing the predecessor and descendant genomes.
It returns the ISs of the predecessor genome in the coordinates of the descendant genome.

### `classify_IS_events_Anc`
Focuses on comparing the descendant genomes with the ancestor of the MA, allowing for the identification of IS events specifically relative to the ancestral state (the genome just before MA).

### `classify_IS_events_MDS`
Compares the descendant genomes with MDS42 + IS1, allowing for the comparison of IS events in all lines by mapping the insertions to the MDS42 genome coordinates.

### `is_position_files`
Describes the position of IS elements within the coordinates of the corresponding genome.
Results of blast hits. `Crossing`: Whether the IS crosses ori (should be False). `Complete`: Whether the IS is full length. `SyRI_genome_id`: the generation ID (1:Anc. 2: MA8, 3: MA20). `SyRI_prefix`: ID of the line.

### `merged` 
Merged output of plotsr-isrm and plotsr.
"local" indicates plotsr using data of nested syntenies.