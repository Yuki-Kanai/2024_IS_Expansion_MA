import logging
import os
import pandas as pd
import ISDetect, SyRI_IS, classify_IS_events

def genome_is_stats(cluster_is_position_file_dir):
	cluster_df = pd.read_csv(cluster_is_position_file_dir, sep=',')
	is_count = len(cluster_df)
	total_is_length = cluster_df.length.sum()
	return is_count, total_is_length

def analyze_is_events(tmpdir, prefix, genome_id):
	r_fasta = os.path.join(tmpdir, prefix, prefix + '_genome.'+str(genome_id-1)+'.fasta')
	q_fasta = os.path.join(tmpdir, prefix, prefix + '_genome.'+str(genome_id)+'.fasta')
	is_classified_df, merged_cluster_df = classify_IS_events.create_merged_cluster_df(r_fasta, q_fasta, os.path.join(tmpdir, prefix), prefix, genome_id)
	is_insertion_cnt = sum(is_classified_df.status == 'simple_insertion')
	is_increasing_event_cnt = merged_cluster_df.query('position_status=="new"')['insert_id'].nunique()
	return is_insertion_cnt, is_increasing_event_cnt

def get_merged_clusters(tmpdir, prefixs, genome_ids):
	assert len(prefixs) == len(genome_ids), 'len(prefixs) != len(genome_ids)'
	merged_cluster_df = pd.DataFrame()
	merged_classification_df = pd.DataFrame()
	for i in range(len(prefixs)):
		prefix, genome_id = prefixs[i], genome_ids[i]
		r_fasta = os.path.join(tmpdir, prefix, prefix + '_genome.'+str(genome_id-1)+'.fasta')
		q_fasta = os.path.join(tmpdir, prefix, prefix + '_genome.'+str(genome_id)+'.fasta')
		is_classified_df_, merged_cluster_df_ = classify_IS_events.create_merged_cluster_df(r_fasta, q_fasta, os.path.join(tmpdir, prefix), prefix, genome_id)
		is_classified_df_['Prefix'] = prefix
		is_classified_df_['genome_id'] = genome_id
		merged_cluster_df_['Prefix'] = prefix
		merged_cluster_df_['genome_id'] = genome_id
		merged_cluster_df = pd.concat([merged_cluster_df, merged_cluster_df_])
		merged_classification_df = pd.concat([merged_classification_df, is_classified_df_])
	return merged_cluster_df, merged_classification_df

def get_genome_sizes(file_df, fastadir, used_fasta_files, lines):
	genome_size_df = pd.DataFrame()
	for i, sample in file_df.iterrows():
		if sample.gen == 'Anc':
			genome_size_df = pd.concat([genome_size_df,
						pd.DataFrame({'Prefix': lines, 'Gen':0, 'File': sample.file_name})])
		elif sample.gen == 'FACS':
			genome_size_df = pd.concat([genome_size_df,
						pd.DataFrame({'Prefix': sample.Prefix, 'Gen': 1, 'File': sample.file_name}, index=[0])])
		else:
			genome_size_df = pd.concat([genome_size_df,
						pd.DataFrame({'Prefix': sample.Prefix, 'Gen': 2, 'File': sample.file_name}, index=[0])])
	genome_size_df = genome_size_df.reset_index(drop=True)
	genome_size_df['File'] = genome_size_df.File.apply(lambda x: os.path.join(fastadir, x))
	genome_size_df.loc[genome_size_df.Gen > 0, 'File'] = used_fasta_files
	genome_size_df['genome_size'] = genome_size_df.File.apply(lambda x: SyRI_IS.get_fasta_length(x, rmN=True))
	return genome_size_df

def get_insertion_distances(lines, exportdir, gids):
	insertion_distance_file = pd.DataFrame()
	classify_IS_event_export = os.path.join(exportdir, 'classify_IS_events')
	for line in lines:
		for gid in gids:
			df_ = pd.read_csv(os.path.join(classify_IS_event_export, line, f'{line}_genome.{gid}.simple_insertion_distances.csv'))
			df_['Line'] = line
			df_['Gen'] = gid
			insertion_distance_file = pd.concat([insertion_distance_file, df_])
	return insertion_distance_file

def get_is_positions(lines, tmpdir, gids):
	is_positions_df = pd.DataFrame()
	for line in lines:
		for gid in gids:
			df_dir = os.path.join(tmpdir, line, f'{line}_clustered_IS_positions_in_subjectgenome.{gid}.csv')
			if not os.path.isfile(df_dir):
				continue
			df_ = pd.read_csv(df_dir)
			df_['Line'] = line
			df_['Gen'] = gid
			is_positions_df = pd.concat([is_positions_df, df_])
	return is_positions_df

def get_is_positions_in_ref_genome(lines, exportdir, gids):
	is_in_r_df = pd.DataFrame()
	for line in lines:
		for gid in gids:
			df_dir = os.path.join(exportdir, line, f'{line}_genome.{gid}.is_flank_seq_pos_in_ref_genome.csv')
			if not os.path.isfile(df_dir):
				continue
			df_ = pd.read_csv(df_dir)
			df_['Line'] = line
			df_['Gen'] = gid
			is_in_r_df = pd.concat([is_in_r_df, df_])
	is_in_r_df.reset_index(drop=True, inplace=True)
	return is_in_r_df

def get_identical_inserts(is_in_r_df, parents):
	identical_inserts_df = pd.DataFrame()
	for parent_ in parents:
		is_in_r_df_ = is_in_r_df[is_in_r_df.Line.str.startswith(parent_) & (is_in_r_df.Gen == 2)].reset_index(drop=True)
		is_in_r_df_clus_ = ISDetect.GetClusters(is_in_r_df_[['pos']], key = 'pos', distance_threshold=20, minSupportNum=1, ignoreStrand=True)
		is_in_r_df_clus_ = pd.concat([is_in_r_df_, is_in_r_df_clus_[1].add_suffix('_cluster')], axis=1)
		# count unique Line per insert_id_cluster
		non_unique_inserts = is_in_r_df_clus_.query('position_status in ["new", "loss"]').\
			groupby('insert_id_cluster').Line.nunique().reset_index().\
			query('Line > 1').insert_id_cluster.unique().astype(int)
		identical_inserts_df = pd.concat([identical_inserts_df,
					is_in_r_df_clus_.query('insert_id_cluster in @non_unique_inserts and position_status in ["new", "loss"]').\
						sort_values('insert_id_cluster')])
	return identical_inserts_df, is_in_r_df_clus_