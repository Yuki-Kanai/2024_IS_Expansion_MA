import pandas as pd
import numpy as np

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import pysam
from statistics import mean

import os
import logging

import pyarrow as pa
import pyarrow.parquet as pq
import pyarrow.feather as feather


import ISDetect, SyRI_IS

def detect_overlap(s1, e1, s2, e2):
	s1_ = min(s1, e1)
	e1_ = max(s1, e1)
	s2_ = min(s2, e2)
	e2_ = max(s2, e2)
	return min(e1_, e2_) - max(s1_, s2_) if (s1_ <= e2_ and s2_ <= e1_) else 0

def range_distance(s1, e1, s2, e2):
	s1_ = min(s1, e1)
	e1_ = max(s1, e1)
	s2_ = min(s2, e2)
	e2_ = max(s2, e2)
	return min(abs(s1_ - e2_), abs(s2_ - e1_)) if not (s1_ <= e2_ and s2_ <= e1_) else 0

# df input is the merged flank df
def annotate_simple_insertions_of_q_is(df, flank_length=100):
	# Split dataframes based on flank
	new_is_lengths = df.apply(lambda x: abs(x['s_start_x'] - x['s_end_x'])+1, axis=1)
	df['new_is_length'] = new_is_lengths
	df['is_matched'] =  df.apply(lambda x: x['new_is_length']*0.9 <= x['is_match_length'], axis=1)
	df_f = df[df['flank'] == 'f'].copy()
	df_r = df[df['flank'] == 'r'].copy()

	# Drop unnecessary columns
	cols_to_drop = ['flank', 'flank_id', 'subject_id', 'new_is_length']
	df_f.drop(columns=cols_to_drop, inplace=True)
	df_r.drop(columns=cols_to_drop, inplace=True)

	# Rename columns
	df_f.columns = [col + '_f' if col not in ['is_id_desc'] else col for col in df_f.columns]
	df_r.columns = [col + '_r' if col not in ['is_id_desc'] else col for col in df_r.columns]

	# Merge dataframes on is_id_desc
	merged_df = df_f.merge(df_r, on='is_id_desc', how='outer')

	merged_df['_is_length'] = merged_df.apply(lambda x: abs(x['s_start_x_f'] - x['s_end_x_f'])+1, axis=1)
	merged_df['_same_strand'] = merged_df.apply(lambda x: True if x['IS_strand_f'] == x['IS_strand_r'] else False, axis=1)
	merged_df['_range_overlap'] = merged_df.apply(lambda x: detect_overlap(x['s_start_y_f'], x['s_end_y_f'], x['s_start_y_r'], x['s_end_y_r']), axis=1)
	merged_df['_range_distance'] = merged_df.apply(lambda x: range_distance(x['s_start_y_f'], x['s_end_y_f'], x['s_start_y_r'], x['s_end_y_r']), axis=1)
	merged_df['_same_ins_pos'] = merged_df.apply(lambda x: abs(x['is_ins_pos_f'] - x['is_ins_pos_r']) <= 20, axis=1)
	merged_df['_detected_is_start'] = merged_df.apply(lambda x: min(x['is_ins_pos_f'], x['is_ins_pos_r']), axis=1)
	merged_df['_detected_is_end'] = merged_df.apply(lambda x: max(x['is_ins_pos_f'], x['is_ins_pos_r']), axis=1)

	merged_df['status'] = 'unknown' # initialize status
	# split the dataframe to avoid calculating nans
	merged_df_wo_nan = merged_df[merged_df.apply(lambda x: not (np.isnan(x['s_start_x_f']) or np.isnan(x['s_start_x_r'])), axis=1)].copy()
	merged_df_w_nan = merged_df[merged_df.apply(lambda x: np.isnan(x['s_start_x_f']) or np.isnan(x['s_start_x_r']), axis=1)]
	# all boolean columns should be boolean
	merged_df_wo_nan['_same_strand'] = merged_df_wo_nan['_same_strand'].astype(bool)
	merged_df_wo_nan['_same_ins_pos'] = merged_df_wo_nan['_same_ins_pos'].astype(bool)
	merged_df_wo_nan['is_matched_f'] = merged_df_wo_nan['is_matched_f'].astype(bool)
	merged_df_wo_nan['is_matched_r'] = merged_df_wo_nan['is_matched_r'].astype(bool)

	# those insertions with flanks significantly matching the IS seqeunces are considered as pristine 
	pristine_mask = (merged_df_wo_nan['_same_strand'] &
				 (merged_df_wo_nan['_range_distance'] < 20) &
				 merged_df_wo_nan['is_matched_f'] &
				 merged_df_wo_nan['is_matched_r'])
	merged_df_wo_nan.loc[pristine_mask, 'status'] = 'pristine'

	# those insertion with flanks partially matching the IS seqeunces are considered as sv within ISs
	is_mutation_criteria = '_same_strand and _range_distance < 20 and (is_match_length_f > _is_length*0.05 or is_match_length_r> _is_length*0.05) and status != "pristine"'
	merged_df_wo_nan.loc[merged_df_wo_nan.eval(is_mutation_criteria), 'status'] = 'mutation_in_is'

	# those insertions with the flanking pairs matched to adjacent sequence in the right direction is considered as a simple insertion
	simple_ins_mask = (merged_df_wo_nan['_same_strand'] & 
					(merged_df_wo_nan['_range_distance'] < 20) & 
					(~merged_df_wo_nan['is_matched_f']) & # sequence similar to IS not found in the inserted position
					(~merged_df_wo_nan['is_matched_r']))
	merged_df_wo_nan.loc[simple_ins_mask, 'status'] = 'simple_insertion'

	merged_df = pd.concat([merged_df_wo_nan, merged_df_w_nan], ignore_index=True).sort_values(by=['is_id_desc']).reset_index(drop=True)

	return merged_df

# Have to change col names to merge
def merge_ref_is_cluster_with_query_is_pos_df_(is_cluster_in_ref_df, ref_pos_of_query_is_df):
	# check all necesary col names exist
	cols_in_ref = ['cluster_id', 'IS_strand', 'start', 'end']
	cols_in_query = ['is_id_desc', 'flank', 'is_ins_pos', 'IS_strand']
	assert all([col in is_cluster_in_ref_df.columns for col in cols_in_ref]), 'Missing columns in is_cluster_in_ref_df: ' + str([col for col in cols_in_ref if col not in is_cluster_in_ref_df.columns])
	assert all([col in ref_pos_of_query_is_df.columns for col in cols_in_query]), 'Missing columns in ref_pos_of_query_is_df: ' + str([col for col in cols_in_query if col not in ref_pos_of_query_is_df.columns])
	# merge the two dfs
	df1 = is_cluster_in_ref_df.melt(id_vars = ['cluster_id', 'IS_strand'], value_vars=['start', 'end'], var_name='flank', value_name='pos') 
	df1['genome'] = 'Ref'
	df1 = df1.rename(columns={'cluster_id': 'pos_id'})
	df1['IS_strand'] = df1['IS_strand'].apply(lambda x: x == 'forward')
	df1['flank'] = df1['flank'].apply(lambda x: 'f' if x == 'start' else 'r')
	df2 = pd.DataFrame({'pos_id': ref_pos_of_query_is_df.is_id_desc, 'flank': ref_pos_of_query_is_df.flank,
		'pos': ref_pos_of_query_is_df.is_ins_pos, 'genome': 'Query', 'IS_strand': ref_pos_of_query_is_df.IS_strand})
	df = pd.concat([df1, df2], ignore_index=True)
	return df

def analyze_whether_is_flank_seq_match_pos_is_new(merged_cluster_df):
	ref_insert_cluster_id = merged_cluster_df[merged_cluster_df['genome'] == 'Ref'].cluster_id.unique()
	# those query IS flanks without matching ref IS flanks are considered as new
	merged_cluster_df['position_status'] = merged_cluster_df.apply(lambda x:
		'original' if x['cluster_id'] in ref_insert_cluster_id else 'new', axis=1)
	query_insert_cluster_id = merged_cluster_df[merged_cluster_df['genome'] == 'Query'].cluster_id.unique()
	# those ref IS flanks without matching query IS flanks are considered as lost
	merged_cluster_df['position_status'] = merged_cluster_df.apply(lambda x:
		'lost_fragment' if (x['position_status'] == 'original' and x['cluster_id'] not in query_insert_cluster_id) else x['position_status'], axis=1)

	ref_genome_is_pos_status = merged_cluster_df[['pos','insert_id', 'position_status']].\
		groupby(['insert_id']).agg({'pos': 'median', 'position_status': 'first'}).reset_index()
	ref_genome_is_pos_status['pos'] = ref_genome_is_pos_status['pos'].astype(int)
	return merged_cluster_df, ref_genome_is_pos_status

def merge_ins_id_info_to_ref_genome_is_pos_status(is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status):
	f_strand_info = is_flank_seq_pos_in_ref_genome_df.query('genome == "Query" and flank =="f"')[['pos_id', 'insert_id', 'position_status']]
	f_strand_info.columns = ['is_id_desc', 'insert_id_f', 'target_locus_status_f']
	r_strand_info = is_flank_seq_pos_in_ref_genome_df.query('genome == "Query" and flank =="r"')[['pos_id', 'insert_id', 'position_status']]
	r_strand_info.columns = ['is_id_desc', 'insert_id_r', 'target_locus_status_r']
	is_classified_df = pd.merge(is_classified_df, f_strand_info, on='is_id_desc', how='left')
	is_classified_df = pd.merge(is_classified_df, r_strand_info, on='is_id_desc', how='left')

	ref_genome_is_pos_status['sv'] = 'unknown'
	pristine_iss = is_classified_df.query('status in ["pristine", "mutation_in_is"]')
	pristine_iss = set(pristine_iss.insert_id_f.tolist() + pristine_iss.insert_id_r.tolist())
	ref_genome_is_pos_status['sv'] = ref_genome_is_pos_status.apply(lambda x: 'pristine' if x['insert_id'] in pristine_iss else x['sv'], axis=1)
	insertion_iss = is_classified_df.query('status == "simple_insertion"')
	insertion_iss = set(insertion_iss.insert_id_f.tolist() + insertion_iss.insert_id_r.tolist())
	ref_genome_is_pos_status['sv'] = ref_genome_is_pos_status.apply(lambda x: 'simple_insertion' if x['insert_id'] in insertion_iss else x['sv'], axis=1)
	is_flank_seq_pos_in_ref_genome_df = pd.merge(is_flank_seq_pos_in_ref_genome_df, ref_genome_is_pos_status[['insert_id', 'sv']], on='insert_id', how='left')
	return is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status

def get_is_flank_seq_pos_in_ref_genome(r_fasta, q_fasta, r_is_cluster, q_is_cluster, tmpdir, prefix, gid_q): 
	# get the flanking sequence of query ISs
	q_flank_fasta_dir = q_fasta.replace('.fasta', '.ISflank.fasta')
	flank_df, q_flank_fasta = SyRI_IS.write_fasta_of_flank_is_pair(q_is_cluster, q_fasta, q_flank_fasta_dir, flank_length=100)

	# map the flanking sequence of query ISs to reference genome
	gid = str(gid_q)
	os.makedirs(os.path.join(tmpdir, 'q_is_in_r_genome', gid), exist_ok=True)
	q_is_in_r_genome = ISDetect.get_insertion_position(flank_df, q_flank_fasta, r_fasta, prefix, tmpdir, distance_threshold=20, savesubdir = os.path.join('q_is_in_r_genome', gid))
	q_is_in_r_genome_raw = q_is_in_r_genome[1][0]

	# merge the query ISs with reference ISs to detect which is the same
	is_merged_by_ref_pos_df = merge_ref_is_cluster_with_query_is_pos_df_(r_is_cluster, q_is_in_r_genome_raw)
	merged_cluster_df = ISDetect.GetClusters(is_merged_by_ref_pos_df, key='pos', distance_threshold=20, minSupportNum=1, ignoreStrand=True)
	merged_cluster_df = merged_cluster_df[1]
	merged_cluster_df, ref_genome_is_pos_status = analyze_whether_is_flank_seq_match_pos_is_new(merged_cluster_df)

	return merged_cluster_df, q_is_in_r_genome_raw, ref_genome_is_pos_status

# NOT USED for creating hypothetical genome with IS insretions
def get_new_is_insertion_specs(ref_mapped_cluster_df, query_mapped_cluster_df, ref_seq, query_seq):
	new_df = ref_mapped_cluster_df[ref_mapped_cluster_df['position_status'] == 'new'].copy()
	new_inserts = new_df.insert_id.unique()
	
	insert_specs = pd.DataFrame(columns=['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO', 'FORMAT', 'SAMPLE'])

	# add dummy specifications to pass bcftools consensus
	insert_specs['QUAL'] = [60] * len(new_inserts)
	insert_specs['FILTER'] = ['PASS'] * len(new_inserts)
	insert_specs['FORMAT'] = ['GT:GQ:DR:DV'] * len(new_inserts) 
	insert_specs['SAMPLE'] = ['0/1:60:0:0'] * len(new_inserts)

	# add IS positions 
	ref_record = [record for record in SeqIO.parse(ref_seq, 'fasta')][0]
	query_record = [record for record in SeqIO.parse(query_seq, 'fasta')][0]
	insert_specs['CHROM'] = ref_record.id
	is_lengths = [0] * len(new_inserts)
	for i, insert_id in enumerate(new_inserts):
		# assume the insert is the sequence of the longest IS match
		df_ = new_df[new_df['insert_id'] == insert_id]
		is_pos_in_ref = int(mean(df_.pos))
		new_is_id_v = df_.pos_id
		q_df = query_mapped_cluster_df.query('cluster_id in @new_is_id_v')
		max_is_length = q_df['length'].max()
		is_lengths[i] = max_is_length
		row = q_df[q_df['length'] == max_is_length].iloc[0]
		insert_specs.loc[i, 'ID'] = f'IS_{insert_id}'
		insert_specs.loc[i, 'POS'] = is_pos_in_ref
		insert_specs.loc[i, 'REF'] = ref_record.seq[is_pos_in_ref-1]
		is_seq = query_record.seq[(row['start']-1):row['end']]
		is_seq = is_seq if row['IS_strand'] == 'forward' else is_seq.reverse_complement()
		assert len(df_.IS_strand.unique()) == 1, 'Multiple IS strands detected for insert ' + str(insert_id) + '. Did you cluster the ISs in query genome?'
		is_strand_in_ref = df_.IS_strand.iloc[0]
		insert_specs.loc[i, 'ALT'] = ref_record.seq[is_pos_in_ref-1] + \
			str(is_seq if is_strand_in_ref else is_seq.reverse_complement())
		insert_specs.loc[i, 'INFO'] = f'END={is_pos_in_ref};SVLEN={max_is_length};SVTYPE=INS'
	is_record_df = pd.DataFrame({'insert_id': new_inserts, 'is_length': is_lengths, 'insert_positions': insert_specs['POS']})
	return insert_specs, is_record_df

#NOT USED for creating hypothetical genome with IS insretions
def create_is_pos_map(is_record_df, is_inserted_genome_size, is_range_buffer = 20):
	is_record_df = is_record_df.sort_values(by=['insert_positions'])
	is_record_df = is_record_df.reset_index(drop=True)
	ref_pos = np.zeros((is_inserted_genome_size,), dtype=int)
	cur_ref_pos = 0
	cur_query_pos = 0
	is_range = np.full(False, is_inserted_genome_size)
	for i, row in is_record_df.iterrows():
		is_pos = row['insert_positions']
		syntenic_pos = range(cur_ref_pos+1, is_pos+1)
		synteny_size = len(syntenic_pos)
		ref_pos[cur_query_pos:(cur_query_pos+synteny_size)] = syntenic_pos
		cur_query_pos += synteny_size
		cur_ref_pos += synteny_size
		is_size = row['is_length']
		ref_pos[(cur_query_pos+1):(cur_query_pos+is_size+2)] = cur_ref_pos
		is_range[min(0, cur_query_pos-is_range_buffer+1):(cur_query_pos+is_size+is_range_buffer+2)] = True
		cur_query_pos += is_size
	is_pos_map = pd.DataFrame({'ref_pos': ref_pos, 'is_range': is_range})	
	return is_pos_map

def create_is_pos_vec(is_df, genome_size, is_range_buffer = 20, default=-1):
	# is_df: insert_id, start, end of ISs. Non-IS sequence is assigned -1
	assert all([col in is_df.columns for col in ['insert_id', 'start', 'end']]), 'Missing columns in is_df: ' + str([col for col in ['insert_id', 'start', 'end'] if col not in is_df.columns])
	is_pos = np.full(genome_size+1, default, dtype=int)
	is_df = is_df.sort_values(by=['start']).reset_index(drop=True)
	for i, row in is_df.iterrows():
		start = max(1, row['start']-is_range_buffer)
		end = min(genome_size, row['end']+is_range_buffer)
		is_pos[start:end+1] = row['insert_id']
	return is_pos

def gen_coord_transformation_vecs(is_df, genome_size, default=-1):
	# is_df: insert_id, start, end of ISs. Non-IS sequence is assigned -1
	assert all([col in is_df.columns for col in ['insert_id', 'start', 'end']]), 'Missing columns in is_df: ' + str([col for col in ['insert_id', 'start', 'end'] if col not in is_df.columns])
	is_df = is_df.sort_values(by=['start']).reset_index(drop=True)
	# check no overlapping region
	assert all(is_df['start'][1:].values > is_df['end'][:-1].values), 'Overlapping ISs detected'

	# predefine the vectors
	fw = np.full(genome_size+1, default, dtype=int)
	is_df['length'] = is_df['end'] - is_df['start'] + 1
	assert all(is_df['length'] > 0), 'IS length <= 0 detected'
	is_rm_genome_size = genome_size - sum(is_df['length'])
	rv = np.full(is_rm_genome_size+1, default, dtype=int)

	fw_pos, rv_pos = 0, 0
	for i, row in is_df.iterrows():
		# add value for non-IS region
		start = row['start']
		end = row['end']
		fw_pos += 1
		rv_pos += 1
		non_is_length = start - fw_pos
		fw[fw_pos:start] = range(rv_pos, rv_pos+non_is_length)
		rv[rv_pos:(rv_pos+non_is_length)] = range(fw_pos, start)
		# add value for IS region
		fw_pos = start
		rv_pos += non_is_length-1
		fw[fw_pos:(end+1)] = rv_pos
		fw_pos = end
	fw[fw_pos+1:] = range(rv_pos+1, is_rm_genome_size+1)
	rv[rv_pos+1:] = range(fw_pos+1, genome_size+1)
	return fw, rv

def save_coord_transformation_vecs(is_df, fasta, default=-1, rewrite=False):
	#max_is_coord = max(is_df['end']) why?
	assert max(is_df['end']) <= len([record for record in SeqIO.parse(fasta, 'fasta')][0].seq), 'IS end position exceeds the length of the genome'
	fw_dir = fasta.replace('.fasta', '.fw_coord_transformation.arrow')
	rv_dir = fasta.replace('.fasta', '.rv_coord_transformation.arrow')
	if os.path.exists(fw_dir) and os.path.exists(rv_dir) and not rewrite:
		logging.info('Coordinate transformation vectors already exist. Skipping... ' + fw_dir + ' and ' + rv_dir)
		return fw_dir, rv_dir
	else:
		fw, rv = gen_coord_transformation_vecs(is_df, SyRI_IS.get_fasta_length(fasta), default=default)
		feather.write_feather(pd.DataFrame({"coord_in_IS_less_genome": fw}), fw_dir, compression='zstd')
		feather.write_feather(pd.DataFrame({"coord_in_original_genome": rv}), rv_dir, compression='zstd')
		logging.info('Coordinate transformation vectors saved to ' + fw_dir + ' and ' + rv_dir)
		return fw_dir, rv_dir

def convert_cluster_to_simple_df(cluster_df):
	return pd.DataFrame({'insert_id': cluster_df['cluster_id'], 'start': cluster_df['start'], 'end': cluster_df['end']})

def convert_merged_cluster_to_simple_df(merged_cluster_df, r_is_cluster):
	assert all([col in r_is_cluster.columns for col in ['insert_id', 'start', 'end']]), 'Missing columns in r_is_cluster: ' + str([col for col in ['insert_id', 'start', 'end'] if col not in r_is_cluster.columns])
	new_is_df_ = merged_cluster_df.query('position_status == "new"')[['insert_id', 'pos']].\
		groupby('insert_id').agg(start=('pos', 'min'), end=('pos', 'max')).reset_index()
	new_is_df_['status'] = 'new'
	r_is_cluster['status'] = 'original'
	r_is_cluster.rename(columns={'insert_id': 'pos_id'}, inplace=True)
	r_is_cluster = pd.merge(merged_cluster_df.query('genome == "Ref" and flank =="f"')[['pos_id', 'insert_id']],
		r_is_cluster, on='pos_id', how='left').drop(columns=['pos_id'])
	merged_is_cluster_ = pd.concat([r_is_cluster, new_is_df_]).sort_values('insert_id').reset_index(drop=True)
	return merged_is_cluster_

def gen_is_pos_conversion_vecs(r_is_cluster, q_is_cluster, is_flank_seq_pos, r_genome_size, q_genome_size, buffer_sizes = [0, 20, 100]):
	r_df = pd.DataFrame()
	r_is_cluster_ = convert_cluster_to_simple_df(r_is_cluster)
	for b in buffer_sizes:
		is_mark_vec = create_is_pos_vec(r_is_cluster_, r_genome_size, is_range_buffer=b)
		r_df['buffer_' + str(b)] = is_mark_vec

	merged_is_cluster_ = convert_merged_cluster_to_simple_df(is_flank_seq_pos, r_is_cluster_)
	m_df = pd.DataFrame()
	for b in buffer_sizes:
		is_mark_vec = create_is_pos_vec(merged_is_cluster_, r_genome_size, is_range_buffer=b)
		m_df['buffer_' + str(b)] = is_mark_vec

	q_df = pd.DataFrame()
	q_is_cluster_ = convert_cluster_to_simple_df(q_is_cluster)
	for b in buffer_sizes:
		is_mark_vec = create_is_pos_vec(q_is_cluster_, q_genome_size, is_range_buffer=b)
		q_df['buffer_' + str(b)] = is_mark_vec

	return r_df, m_df, q_df

def save_is_pos_conversion_vecs(r_is_cluster, q_is_cluster, is_flank_seq_pos, r_genome_size, q_genome_size, save_dir, prefix, rid, qid, buffer_sizes = [0, 20, 100], skip_wrote=True):
	os.makedirs(save_dir, exist_ok=True)
	is_marked_coords = gen_is_pos_conversion_vecs(r_is_cluster, q_is_cluster, is_flank_seq_pos, r_genome_size, q_genome_size, buffer_sizes = buffer_sizes)
	labels = [f'ref_{rid}', f'merged_{rid}_{qid}', f'query_{qid}']
	all_saved = True
	file_dirs = [os.path.join(save_dir, prefix + '_' + label + '_is_pos_conversion.arrow') for label in labels]
	for i, df in enumerate(is_marked_coords):
		file_dir = file_dirs[i]
		if not skip_wrote or not os.path.exists(file_dir):
			feather.write_feather(df, file_dir, compression='zstd')
			all_saved = False
	if all_saved:
		logging.info('IS position conversion vectors already exist. Skipping... ')
	else:
		logging.info('IS position conversion vectors saved to ' + save_dir)
	file_dirs = dict(zip(labels, file_dirs))
	return file_dirs

def compare_length_(start1, end1, start2, end2, min_sv_length):
	cond1 = False 
	if start1 != "-" and end1 != "-" :
		cond1 = int(end1) - int(start1) >= min_sv_length

	cond2 = False 
	if start2 != "-" and end2 != "-" :
		cond2 = abs(int(end2) - int(start2)) >= min_sv_length
	return cond1 or cond2

def remove_short_sv_syri_out(syri_out, min_sv_length):
	min_sv_length = min_sv_length  - 1
	syri_out = syri_out[syri_out['SVTYPE'] != "SNP"]
	condition2 = syri_out.apply(lambda row: compare_length_(row['start1'], row['end1'], row['start2'], row['end2'], min_sv_length), axis=1)
	syri_out = syri_out[condition2]
	return syri_out

def classify_sv_boundary_is_in_SyRI_out(syri_out_df, is_pos_vec_r, is_pos_vec_q):
	is_pos_vec_r = pd.read_feather(is_pos_vec_r).iloc[:,1].values # use most strict buffered IS position
	is_pos_vec_q = pd.read_feather(is_pos_vec_q).iloc[:,1].values # use buffered is position
	syri_out_df = syri_out_df.copy()
	syri_out_df['start1_is'] = syri_out_df.start1.apply(lambda x: is_pos_vec_r[int(x)] if str(x).isdigit() else -1)
	syri_out_df['end1_is'] = syri_out_df.end1.apply(lambda x: is_pos_vec_r[int(x)] if str(x).isdigit() else -1)
	syri_out_df['start2_is'] = syri_out_df.start2.apply(lambda x: is_pos_vec_q[int(x)] if str(x).isdigit() else -1)
	syri_out_df['end2_is'] = syri_out_df.end2.apply(lambda x: is_pos_vec_q[int(x)] if str(x).isdigit() else -1)
	return syri_out_df

def read_syri_out(syri_out_dir):
	syri_out = pd.read_csv(syri_out_dir, sep='\t', header=None, dtype=str)
	syri_out.columns = ['chr1', 'start1', 'end1', 'REF', 'ALT', 'chr2', 'start2', 'end2', 'ID', 'parent', 'SVTYPE', 'MISC']
	return syri_out

# is in the reference genome has two insertion positions assigned. this function recovers the other one given list of is ids interms of merged cluster ids 
def recover_matching_ref_is_pos(is_ids, ref_is_pos_status):
	is_in_ref = ref_is_pos_status.query('insert_id in @is_ids and genome == "Ref"').pos_id
	return list(set((list(is_ids) + list(ref_is_pos_status.query('pos_id in @is_in_ref and genome == "Ref"').insert_id.values))))

def get_closest_index_and_distance_(df, value, pos='pos', id='insert_id'):
	diff = np.abs(df[pos] - value)
	idx = diff.idxmin()
	return df[id][idx], diff.iloc[idx]

def analyze_simple_insertion_distances(ref_genome_is_pos_status):
	original_pos = ref_genome_is_pos_status.query('position_status in ["original", "lost_fragment"]').reset_index(drop=True)
	new_insertions = ref_genome_is_pos_status.query('position_status == "new"').reset_index(drop=True)
	new_insertions[['closest_original_id', 'closest_original_distance']] =\
		new_insertions.apply(lambda x: get_closest_index_and_distance_(original_pos, x['pos']), axis=1, result_type='expand')
	return new_insertions

def annotate_simple_transpositions(is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status, syri_out):
	for row in syri_out.query('SVTYPE in ["TRANSAL", "INVTRAL", "DUPAL"]').itertuples():
		if not (row.start2_is != -1 and row.end2_is != -1):
			continue
			# miss when the destination is is lost
		# get the is position in the ref genome
		from_is = [row.start1_is, row.end1_is]
		from_is_bounds = recover_matching_ref_is_pos(from_is, is_flank_seq_pos_in_ref_genome_df)
		f_is_in_ref = is_classified_df.iloc[row.start2_is,][['insert_id_f', 'insert_id_r']].values.flatten()
		r_is_in_ref = is_classified_df.iloc[row.end2_is,][['insert_id_f', 'insert_id_r']].values.flatten()
		f_notmatching = list(set(f_is_in_ref) - set(from_is_bounds))
		r_notmatching = list(set(r_is_in_ref) - set(from_is_bounds))
		# both ends of the aligning region have one end inside the IS and the other end
		# outside the IS is identical between the pair of ends
		# and this identical outside IS does not match ISs that were present in the reference genome
		if len(f_notmatching) == 1 and len(r_notmatching) == 1 and f_notmatching[0] == r_notmatching[0] and f_notmatching[0] != -1:
			# the case where two iss are two close and not aligned is also removed
			ref_genome_is_pos_status['sv'] = ref_genome_is_pos_status.apply(lambda x: 'simple_transposition' if x['insert_id'] == f_notmatching[0] else x['sv'], axis=1)
			is_flank_seq_pos_in_ref_genome_df['sv'] = is_flank_seq_pos_in_ref_genome_df.apply(lambda x: 'simple_transposition' if x['insert_id'] == f_notmatching[0] else x['sv'], axis=1)
			is_classified_df['status'] = is_classified_df.apply(lambda x: 'simple_transposition' if x['is_id_desc'] in [row.start2_is, row.end2_is] else x['status'], axis=1)
	return 0

def aggregate_ref_is_flank_pair_of_merged_df(merged_df, reference_cluster_df):
	retain_cols = ['pos_id', 'IS_strand', 'pos', 'insert_id', 'cluster_id']
	new_df = pd.merge(merged_df.query('genome=="Ref" and flank=="f"')[retain_cols],
		merged_df.query('genome=="Ref" and flank=="r"')[retain_cols],
		on='pos_id', how='outer')
	new_df = pd.merge(new_df, reference_cluster_df[['cluster_id', 'start', 'end', 'IS_strand', 'length']],
		left_on='pos_id', right_on='cluster_id', how='right').drop(columns=['cluster_id'])
	return new_df

def roundDepth_(depth, rangelength, depthThresholdRatio=0.05, depthThresholdLength=100, depthThresholdRatio2=0.2):
	threshold_depth = min(depthThresholdRatio2, max(float(depthThresholdLength)/rangelength, depthThresholdRatio))
	depth_ = -1 
	if depth < threshold_depth:
		depth_ = 0
	elif 1 - threshold_depth < depth < 1 + threshold_depth:
		depth_ = 1
	if  2 - threshold_depth * 2 < depth:
		depth_ = round(depth)
	return depth_

def calc_cn_changed_intervals(merged_flank_insert_df, r_is_clusters_df, r_masked_fasta, q_masked_fasta, prefix, tmpdir):
	insert_id_pos_df = merged_flank_insert_df[['insert_id', 'pos']]
	insert_id_pos_df.columns = ['insert_id', 'medianpos']
	insert_id_pos_df = insert_id_pos_df.groupby('insert_id').median().reset_index()
	insert_id_pos_df['medianpos'] = insert_id_pos_df['medianpos'].astype(int)
	samdir = q_masked_fasta.replace('.fasta', '.sam')
	intervalDepth, depthDF = ISDetect.CalcDepthInterval(insert_id_pos_df, samdir, prefix, os.path.join(tmpdir, prefix), q_masked_fasta, r_masked_fasta,
		refgenomeLength = SyRI_IS.getFastaLength(r_masked_fasta), skipalignment=True)
	intervalDepth['RoundedDepth'] = intervalDepth.apply(lambda x: roundDepth_(x.depth, x.length), axis=1)

	# determine which interval corresponds to is_sequences preexisting in reference genome
	intervalDepth['r_is_id'] = -1
	for i, row in r_is_clusters_df.iterrows():
		if row['insert_id_y'] - row['insert_id_x'] != 1:
			print(row)
			logging.warn('insert_id_y - insert_id_x != 1')
			continue
		intervalDepth.loc[intervalDepth['interval_id'] == row['insert_id_y'], 'r_is_id'] = row['pos_id']
	return intervalDepth

def idenfity_cn_changed_range_status(intervalDepth, merged_cluster_df):
	cn_changed_intervals = intervalDepth.query('r_is_id == -1 and RoundedDepth != 1')
	cn_changed_intervals['new_cnt'] = 0
	cn_changed_intervals['q_is_id'] = -1
	for i, row in cn_changed_intervals.iterrows():
		id = row['interval_id']
		max_id = max(merged_cluster_df['insert_id'])
		new_cnt = 0
		q_is_id = merged_cluster_df.query('insert_id == @id and genome=="Query"').pos_id.values[0]
		cn_changed_intervals.loc[i, 'q_is_id1'] = q_is_id
		q_is_id = merged_cluster_df.query('insert_id == @id and genome=="Query"').pos_id.values[0]
		cn_changed_intervals.loc[i, 'q_is_id2'] = q_is_id
		# get the number of new is insertions in terms of ref genome positions
		if id != 1:
			new_cnt += len(merged_cluster_df.query('insert_id in [@id - 1, @id] and position_status=="new"'))
		else:
			new_cnt += len(merged_cluster_df.query('insert_id == [@id, @max_id] and position_status=="new"'))
		cn_changed_intervals.loc[i, 'new_cnt'] = new_cnt
	return cn_changed_intervals

def create_N_rm_genome(input, output):
	if not os.path.isfile(output):
		masked_fasta = SeqIO.parse(input, 'fasta')
		masked_fasta = [SeqRecord(Seq(str(rec.seq).replace('N', '')), id=rec.id, description='') for rec in masked_fasta]
		SeqIO.write(masked_fasta, output, 'fasta')
		logging.info(f'created {output}')
	return output

def delete_is_positions(is_df,input, output, rewrite=False, contig_name='contig_1'):
	if not os.path.isfile(output) or rewrite:
		is_df = is_df.copy().sort_values(by=['start']).reset_index(drop=True)
		is_df['length'] = is_df['end'] - is_df['start'] + 1
		masked_fasta = SeqIO.parse(input, 'fasta')
		seq = [rec.seq for rec in masked_fasta][0]
		new_seq = ''
		coor_ref = 0
		for i, row in is_df.iterrows():
			new_seq += seq[coor_ref:(row['start']-1)]
			coor_ref = row['end']
		new_seq += seq[coor_ref:]
		assert len(new_seq) == len(seq) - sum(is_df['length']), 'Length of the new sequence is not equal to the original sequence minus the length of ISs'
		SeqIO.write(SeqRecord(Seq(new_seq), id=contig_name, description=''), output, 'fasta')
		logging.info(f'created {output}')
	return output


def reverse_lift_up_syri_out_(syri_out_df, ref_trans_rv, query_trans_rv, suffix = '', add_rows = False):
	if add_rows:
		syri_out_df['start1' + suffix] = syri_out_df['start1']
		syri_out_df['end1' + suffix] = syri_out_df['end1']
		syri_out_df['start2' + suffix] = syri_out_df['start2']
		syri_out_df['end2' + suffix] = syri_out_df['end2']
	syri_out_df['start1'] = syri_out_df['start1'].apply(lambda x: ref_trans_rv.iloc[int(x),0] if x.isdigit() else x)
	syri_out_df['end1'] = syri_out_df['end1'].apply(lambda x: ref_trans_rv.iloc[int(x),0] if x.isdigit() else x)
	syri_out_df['start2'] = syri_out_df['start2'].apply(lambda x: query_trans_rv.iloc[int(x),0] if x.isdigit() else x)
	syri_out_df['end2'] = syri_out_df['end2'].apply(lambda x: query_trans_rv.iloc[int(x),0] if x.isdigit() else x)
	return syri_out_df

def reverse_lift_up_syri_out(syri_out_dir, ref_trans_rv, query_trans_rv):
	syri_out_df = pd.read_csv(syri_out_dir, sep='\t', header=None, dtype='str')
	syri_out_df.columns = ['chr1', 'start1', 'end1', 'REF', 'ALT', 'chr2', 'start2', 'end2', 'ID', 'parent', 'SVTYPE', 'MISC']
	syri_out_df = reverse_lift_up_syri_out_(syri_out_df, ref_trans_rv, query_trans_rv)
	new_syri_dir = syri_out_dir.replace('.syri.out', '.rv_coord.syri.out')
	syri_out_df.to_csv(new_syri_dir, sep='\t', header=False, index=False)
	return new_syri_dir

def fix_syri_sv_boundaries_of_syri_out(syri_out_is_annotated, r_is_cluster, q_is_cluster, is_flank_seq_pos_in_ref_genome_df, is_classified_df, ref_genome_is_pos_status):
	# Get only the SVs that are parental
	subset_df = syri_out_is_annotated.query('not SVTYPE in ["SYN"]').\
		query('not SVTYPE.str.endswith("AL") or SVTYPE == "NOTAL"').\
			query('chr1 != "-"')
	# start end of reference should be int
	subset_df['start1'] = subset_df.start1.astype(int)
	subset_df['end1'] = subset_df.end1.astype(int)

	sv_df = subset_df.assign(
		start1_fixed=subset_df.start1,
		end1_fixed=subset_df.end1,
		start2_fixed=subset_df.start2,
		end2_fixed=subset_df.end2
	)

	# we consider NOTALs as DELs that were not resolved by SyRI.
	# Typically it should be those deletions next to other rearrangements
	# such as inversions or translocations (not checked).
	sv_df['SVTYPE'] = sv_df.SVTYPE.replace('NOTAL', 'DEL')
	sv_df.reset_index(drop=True, inplace=True)

	# note that for ISs found in the reference genome,
	# only the status of the left flank is considered
	# the right flank is considered as the status of the next IS flank ordered by the position in the reference genome
	def get_status(x, ref_is_pos_stat, stat_str='position_status', id_str='insert_id'):
		return ref_is_pos_stat[stat_str][
			ref_is_pos_stat[id_str] == x].values[0]\
				if ref_is_pos_stat[id_str].isin([x]).any() else 'NA' 
	for col in ['start1_is', 'end1_is']:
		sv_df[col+'_status'] = sv_df[col].apply(lambda x: get_status(x, ref_genome_is_pos_status))
		sv_df[col+'_status_right'] = sv_df[col].apply(lambda x: get_status(x+1, ref_genome_is_pos_status))
	for col in ['start2_is', 'end2_is']:
		for fr in ['f', 'r']:
			sv_df[col+'_status_'+fr] = sv_df[col].apply(lambda x: get_status(x, is_classified_df, stat_str='target_locus_status_'+fr, id_str='is_id_desc'))
	sv_df['new_cnt'] = sv_df.apply(lambda x: (1 if x.start1_is_status == 'new' else 0) + (1 if x.end1_is_status == 'new' else 0), axis=1)

	# fixes the boundaries of the SVs so that SV length does not include the boundary ISs
	for r in sv_df.itertuples():
		id = r[0]
		r_is_flank_df = is_flank_seq_pos_in_ref_genome_df.query('genome == "Ref"')

	# fix sv boundaries in reference genome
		if r.start1_is != -1:
			r_is_id = r_is_flank_df.query('insert_id == @r.start1_is').pos_id.values
			if len(r_is_id) > 0: # see if it was in the reference
				r_is_id = r_is_id[0]
				# check if the left flank is lost. If so it implies that the IS is lost with deletion
				if r.SVTYPE == 'DEL' and r.start1_is_status == 'lost':
					sv_df.loc[id, 'start1_fixed'] = r_is_cluster.start[r_is_id]
				else: # otherwise, set the boundary to the right flank
					sv_df.loc[id, 'start1_fixed'] = r_is_cluster.end[r_is_id] + 1
		if r.end1_is != -1:
			r_is_id = r_is_flank_df.query('insert_id == @r.end1_is').pos_id.values
			if len(r_is_id) > 0:
				r_is_id = r_is_id[0]  # see if it was in the reference
				if r.SVTYPE == 'DEL' and r.end1_is_status_right == 'lost':
					sv_df.loc[id, 'end1_fixed'] = r_is_cluster.end[r_is_id]
				else:
					sv_df.loc[id, 'end1_fixed'] = r_is_cluster.start[r_is_id] - 1
		
		# fix sv boundaries in query genome
		# the orientation is based on corresponding sequence in the reference genome
		# thus we must care about the direction
		#  for cases with <IS1> -- SV -- <IS2>
		q_is_flank_df = is_flank_seq_pos_in_ref_genome_df.query('genome == "Query"') 
		if r.start2_is < r.end2_is:
			if r.start2_is != -1:
				sv_df.loc[id, 'start2_fixed'] = q_is_cluster['end'][r.start2_is] + 1
			if r.end2_is != -1:
				sv_df.loc[id, 'end2_fixed'] = q_is_cluster['start'][r.end2_is] - 1
		else:
			if r.start2_is != -1:
				sv_df.loc[id, 'start2_fixed'] = q_is_cluster['start'][r.start2_is] - 1
			if r.end2_is != -1:	
				sv_df.loc[id, 'end2_fixed'] = q_is_cluster['end'][r.end2_is] + 1
		# if  <IS> * SV , such as TSD
		if (r.start2_is != -1) and (r.start2_is == r.end2_is):
			sv_df.loc[id, 'start2_fixed'] = sv_df.loc[id, 'start2']
			sv_df.loc[id, 'end2_fixed'] = sv_df.loc[id, 'end2']

	sv_df = sv_df.assign(
		dstart1_ = lambda x: x.start1_fixed - x.start1,
		dend1_ = lambda x: x.end1_fixed - x.end1,
		dstart2_ = lambda x: x.start2_fixed.replace('-', np.nan).astype(float) - x.start2.replace('-', np.nan).astype(float),
		dend2_ = lambda x: x.end2_fixed.replace('-', np.nan).astype(float) - x.end2.replace('-', np.nan).astype(float),
		SVLEN = lambda x: x.end1 - x.start1 + 1,
		SVLEN_fixed = lambda x: x.end1_fixed - x.start1_fixed + 1,
	)
	sv_df['SVLEN'] = sv_df.apply(lambda x: x.SVLEN if not(x.SVTYPE == 'INS' and x.SVLEN < 20) else abs(x.end2-x.start2)+1, axis=1)
	sv_df['SVLEN_fixed'] = sv_df.apply(lambda x: x.SVLEN_fixed if not(x.SVTYPE == 'INS' and x.SVLEN_fixed < 20) else x.SVLEN, axis=1)
	return sv_df

def combind_insertions_to_is_less_syri_out():
	return 0

def detect_is_loss():
	# is is lost if DEL includes the IS (DEL) or if 
	return 0

def run_syri_isrm_for_specific_id(tmpdir, prefix, rid, qid):
	tmpdir_isrm = os.path.join(tmpdir, prefix, 'is_deleted')
	q_fasta_isrm = os.path.join(tmpdir_isrm, f'{prefix}_genome.deleted_IS.{qid}.fasta')
	r_fasta_isrm = os.path.join(tmpdir_isrm, f'{prefix}_genome.deleted_IS.{rid}.fasta')
	bam_dir = os.path.join(tmpdir, prefix, 'bam', f'{prefix}_genome.deleted_IS.{rid}_{qid}.bam')
	SyRI_IS.align_two_genomes(r_fasta_isrm, q_fasta_isrm, bam_dir,
				preset = 'asm20', option = '-H -f 100 -r1k,10k --rmq=no', output_type = 'bam')

	syri_dir, syri_out_file = SyRI_IS.run_syri(bam_dir, r_fasta_isrm, q_fasta_isrm, f'{prefix}.{rid}_{qid}.', tmpdir_isrm, mapping_type = 'bam')
	#L01-2_genome.0.fw_coord_transformation.arrow
	q_coor_trans_rv = pd.read_feather(os.path.join(tmpdir, prefix, f'{prefix}_genome.{qid}.rv_coord_transformation.arrow'))
	r_coor_trans_rv = pd.read_feather(os.path.join(tmpdir, prefix, f'{prefix}_genome.{rid}.rv_coord_transformation.arrow'))
	new_syri_out_file = reverse_lift_up_syri_out(syri_out_file, r_coor_trans_rv, q_coor_trans_rv)
	SyRI_IS.convert_syri_cpg_cpl(new_syri_out_file)

def get_depth_intervals_from_bam(bam_dir, max_gap = 20):
	logging.info('Calculating depth...')
	pysamOut = pysam.depth('-a', bam_dir, split_lines=True)
	depthDF = pd.DataFrame([x.split('\t') for x in pysamOut])
	depthDF = depthDF.drop(depthDF.columns[[0]], axis=1).set_axis(['locus', 'depth'], axis =1)
	depthDF = depthDF.astype(int)

	# get intervals with same depth
	ser = depthDF.depth
	grp_ser = ser.groupby((ser.diff() !=0).cumsum()).transform('size') # each component of the group has value size of the interval

	# Identify small groups and try setting it to '1'
	mask_small_groups = grp_ser < max_gap 
	depthDF.loc[mask_small_groups, 'depth'] = 1

	# if filling gaps with 1 still leaves a small group as a gap, then set it to nan to gap fill with close values
	ser = depthDF.depth
	grp_ser = ser.groupby((ser.diff() !=0).cumsum()).transform('size')
	mask_small_groups = grp_ser < max_gap 
	depthDF.loc[mask_small_groups, 'depth'] = np.nan
	depthDF['depth'] = depthDF.depth.ffill().bfill()

	# assign group ids from 1
	ser = depthDF.depth
	depthDF['group'] = (ser.diff() !=0).cumsum()

	logging.info('Depth calculation finished')

	return depthDF.groupby('group').\
	 agg(start=('locus', 'min'), end=('locus', 'max'), depth=('depth', 'mean')).\
		reset_index().\
		assign(length = lambda x: x.end - x.start + 1)

def assign_is_id_to_interval_bound(depthDF,  r_coor_trans_rv, is_pos_conv_vec_ref):
	depthDF['start_fix'] = depthDF.start.apply(lambda x: r_coor_trans_rv.coord_in_original_genome[x])
	depthDF['end_fix'] = depthDF.end.apply(lambda x: r_coor_trans_rv.coord_in_original_genome[x])
	depthDF['is_start'] = depthDF.start_fix.apply(lambda x: is_pos_conv_vec_ref.iloc[x,1])
	depthDF['is_end'] = depthDF.end_fix.apply(lambda x: is_pos_conv_vec_ref.iloc[x,1])
	return depthDF
