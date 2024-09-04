## File Descriptions and Column Details

### `IS_positions.csv`
The position of the IS in its corresponding genome.
- **`start, end`**: The start and end positions of the IS on the chromosome.
- **`IS_strand`**: The strand orientation of the IS.
- **`max_alignment_length`**: The largest match by a BLAST search.
- **`length`**: The length of the IS.
- **`cluster_id`**: The ID of the specific IS.
- **`Line, Gen`**: Identifiers for the MA line and generation number.

### `IS_positions_in_ref_genome.csv` and `IS_positions_in_ref_genome_identical.csv`
Positions of the ends of the ISs in the reference genome coordinates.
- **`pos_id`**: Identifier for the position of the IS as in `IS_positions.csv`: `cluster_id`.
- **`IS_strand`**: Strand orientation of the IS.
- **`flank`**: Which flank of the IS is being referred to?
- **`pos`**: The coordinates in the reference genome.
- **`genome`**: Indicates whether the IS is of the reference or the descendant genome.
- **`cluster_id, insert_id`**: Identifiers for clusters of inserts and reordered cluster ID.
- **`flankMatchDir`**: Indicates whether the two flanking sequences match directions, used to identify inversions.
- **`position_status`**: Whether the insertion is new or not.
- **`sv`**: Corresponding structural variant of the flank.
- **`Line, Gen`**: Identifiers for the MA line and generation number.
- For `IS_positions_in_ref_genome_identical.csv`, includes `pos_cluster, cluster_id_cluster, insert_id_cluster` to identify identical new IS insertions across multiple lines. Only identical IS insertions are included and the dataframe should be empty in the final version of the analysis.

### `depth_interval.csv`
Result of depth analysis using `pysam` based on genome sequences with IS hard masks.
- **`group`**: ID of the interval.
- **`start, end`**: Start and end positions of the interval in hard-masked genome coordinates.
- **`depth`**: Average sequencing depth across the interval.
- **`length`**: Length of the interval.
- **`start_fix, end_fix`**: Coordinates after recovering the ISs from the hard-masked genome sequence.
- **`is_start, is_end`**: ID of the IS corresponding to the ends of the interval; -1 if there is no corresponding IS.
- **`line, gen_id`**: Identifiers for the line and generation.

### `excision_status.3000.csv` and `excision_flank_blast.3000.csv`
Output of `20231213_is_excision_pipeline.ipynb` with details on excision events and BLAST results.
- **`Line, Gen`**: Identifiers for the MA line and generation number.
- **`is_cluster_id`**: ID of the IS.
- **`dir`**: Direction of the sequence read.
- **`read_name, sseqid`**: Read and subject sequence identifiers.
- **`qstart, qend, sstart, send`**: Start and end positions of the query and subject in the alignment.
- **`length`**: Length of the alignment.
- **`evalue, bitscore, pident`**: BLAST alignment metrics.
- For `excision_status.3000.csv`, includes `min_dist, head_qend, tail_qstart, min_i, min_j, raw_dist, min_start_raw, min_end_raw, flank_dist` to detail the excision analysis, focusing on the minimal distance between flanks and adjustments for IS lengths. When the `min_dist` is smaller than the IS length, manual checks were made to confirm the excision event.

### `genome_stats.csv` and `genome_stats_merged.csv`
Data including genome size and the number of IS-related events.
- **`IS_Detect_ID`**: Unique identifier for the sample.
- **`ParentLine, SubLine, gen`**: Identifiers for the lineage and generation of the genome.
- **`file_name, File`**: Names of the files associated with the genome data.
- **`Anc, RecA`**: Indicates if the genome is ancestral or has undergone recombination.
- **`Prefix, sample_name_raw`**: Prefix and raw sample names used in the analysis.
- **`Contig_Date, Complete`**: Date of contig assembly and completeness status.
- **`gen_id`**: The generation ID (1: Ancestor, 2: MA8, 3: MA20).
- **`genome_size`**: Total size of the genome in base pairs.
- **`is_cnt`**: Count of IS elements detected.
- **`is_length`**: Total length of all IS elements combined.
- **`ins_pos_cnt_raw`**: Number of insertion sites, including duplicates.
- **`ins_pos_cnt`**: Number of insertion sites, excluding duplicates.
- **`new_cnt, lost_end_cnt, simple_insertion_cnt`**: Counts of new IS insertions, lost ends, and simple insertions, respectively.
- Additional columns in `genome_stats_merged.csv` relate to classified IS events in terms of IS insertion sites based on the coordinates of the MA ancestor and the MDS42 + IS1 reference genome.

### `simple_insertion_distances.csv`
Data frame used to analyze local hopping of IS elements.
- **`insert_id`**: Unique identifier for the IS insertion site.
- **`pos`**: Position of the insertion in the genome.
- **`position_status`**: Status of the insertion (e.g., new or not).
- **`sv`**: simple insertion or not.
- **`closest_original_id`**: ID of the closest original IS element in terms of reference genome coordinates.
- **`closest_original_distance`**: Distance to the closest original IS element.

### `sv_df.csv`
Output of SV (Structural Variation) detection, including details from SYRI analysis.
- **`chr1, start1, end1, REF, ALT, chr2, start2, end2, ID, parent, SVTYPE, MISC`**: SYRI output.
- **`*_is_rm`**: Indicates chromosome coordinates in terms of IS-hardmasked genomes.
- **`start1_is, end1_is, start2_is, end2_is`**: Corresponding IS `pos_id` if exists.
- **`SVLEN_ISDEL`**: Length of the SV in the hardmasked genome.
- **`*_fixed`**: Coordinates accounting for the boundaries of ISs according to the SV type.
- **`xx_status`**: Whether the corresponding IS is new or not.
- **`new_cnt`**: Count of new structural variations detected.
- **`SVLEN, SVLEN_fixed`**: Length of the SV, and its adjusted length considering IS boundaries.
