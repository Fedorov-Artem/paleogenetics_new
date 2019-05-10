#This is the scrip to check bam file for known Y-SNPs

import bamnostic as bs
import numpy as np
import pandas as pd

bam_file_path = 'project_data/Botai/TU/TU45.mapped.rmdup.q30.bam'
#bam_file_path = 'project_data/Cheddar_man/SB524A4_lib.merged.markdup.bam'
target_snps_path = 'project_data/R_fixed.csv'
output_file_name = 'project_data/botai_tu45.csv'

bam = bs.AlignmentFile(bam_file_path, 'rb')
df = pd.read_csv(target_snps_path)
# это - номер хромосомы, 0 - это первая хромосома, 21 - 22-ая, 22 - X, 23 - Y
tid_number = 23
df['result_A'] = 0
df['result_T'] = 0
df['result_G'] = 0
df['result_C'] = 0


'''
for bam_reads in bam:
    if bam_reads.tid == tid_number:
        if (bam_reads.pos <= 7792789) and (bam_reads.pos + bam_reads.l_seq > 7792789):
                print(bam_reads.query_sequence)
                relative_pos = 7792789 - bam_reads.pos
                print(relative_pos)
                relative_pos1 = relative_pos + 1
                if len(bam_reads.query_sequence) > relative_pos:
                    print(bam_reads.query_sequence[relative_pos])
                if len(bam_reads.query_sequence) > relative_pos1:
                    print(bam_reads.query_sequence[relative_pos1])

'''

for bam_reads in bam:
    if bam_reads.tid == tid_number:
        for i in range(0, len(df)):
            if (bam_reads.pos < df['build_37'][i]) and (bam_reads.pos + bam_reads.l_seq >= df['build_37'][i]):
                print(bam_reads.query_sequence)
                relative_pos = df['build_37'][i] - bam_reads.pos - 1
                print(relative_pos)
                if relative_pos >= 0:
                    print(bam_reads.query_sequence[relative_pos], df['Haplogroup'][i], df['ancestral'][i], df['derived'][i])
                    if bam_reads.query_sequence[relative_pos] == 'A':
                        df['result_A'][i] = df['result_A'][i] + 1
                    if bam_reads.query_sequence[relative_pos] == 'T':
                        df['result_T'][i] = df['result_T'][i] + 1
                    if bam_reads.query_sequence[relative_pos] == 'G':
                        df['result_G'][i] = df['result_G'][i] + 1
                    if bam_reads.query_sequence[relative_pos] == 'C':
                        df['result_C'][i] = df['result_C'][i] + 1
                else:
                    print('ERROR!!!!!')

df_export1 = df[df['result_A'] + df['result_T'] + df['result_G'] + df['result_C'] > 0]
df_export = df_export1.copy()

df_export['result_text'] = '0'
df_export['result_text'] = df_export[['result_A','result_T','result_G','result_C']].idxmax(axis=1)
df_export['result_text'] = df_export['result_text'].str.strip().str[-1]

#pd.options.mode.chained_assignment = None
df_export['result_final'] = 'new'
#df_export[df_export['result_text'] == df_export['ancestral']]['result_final'] = 'ancestral'
#df_export[df_export['result_text'] == df_export['derived']]['result_final'] = 'derived'
df_export.loc[df_export['result_text'] == df_export['ancestral'], ['result_final']] = 'ancestral'
df_export.loc[df_export['result_text'] == df_export['derived'], ['result_final']] = 'derived'

df_export.to_csv(output_file_name, encoding='utf-8', index=False)