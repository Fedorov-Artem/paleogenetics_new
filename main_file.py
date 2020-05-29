#This is the script to check bam file for known Y-SNPs

import bamnostic as bs
import numpy as np
import pandas as pd
import time


bam_file_path = r'D:\Workspace\paleogenetics\project_data\med_isles\I14675.hg19.bam'
#bam_file_path = 'project_data/Cheddar_man/SB524A4_lib.merged.markdup.bam'
target_snps_path = 'project_data/trunk_and_i_new.csv'
output_file_name = 'project_data/med_isles_I14675.csv'

bam = bs.AlignmentFile(bam_file_path, 'rb')
df = pd.read_csv(target_snps_path)
# это - номер хромосомы, 0 - это первая хромосома, 21 - 22-ая, 22 - X, 23 - Y
tid_number = 23
df['result_A'] = 0
df['result_T'] = 0
df['result_G'] = 0
df['result_C'] = 0
df['result_A_end'] = 0
df['result_T_end'] = 0
df['result_G_end'] = 0
df['result_C_end'] = 0

build = 'build_37'
start = time.time()


for bam_reads in bam:
    if bam_reads.tid == tid_number:
        first = bam_reads.pos
        last = bam_reads.pos + bam_reads.l_seq
        if len(df.loc[(first < df[build]) & (last >= df[build])]) > 0:
            for i in df.loc[(first < df[build]) & (last >= df[build])].index:
                absolute_pos = df[build][i]
                print(bam_reads.query_sequence)
                relative_pos = absolute_pos - first - 1
                print(relative_pos)
                if relative_pos >= 0:
                    print(bam_reads.query_sequence[relative_pos], df['Haplogroup'][i], df['ancestral'][i], df['derived'][i])
                    if bam_reads.query_sequence[relative_pos] == 'A':
                        df.loc[i, 'result_A'] += 1
                        if relative_pos < 4 or bam_reads.l_seq - relative_pos < 4:
                            df.loc[i, 'result_A_end'] += 1
                    elif bam_reads.query_sequence[relative_pos] == 'T':
                        df.loc[i, 'result_T'] += 1
                        if relative_pos < 4 or bam_reads.l_seq - relative_pos < 4:
                            df.loc[i, 'result_T_end'] += 1
                    elif bam_reads.query_sequence[relative_pos] == 'G':
                        df.loc[i, 'result_G'] += 1
                        if relative_pos < 4 or bam_reads.l_seq - relative_pos < 4:
                            df.loc[i, 'result_G_end'] += 1
                    elif bam_reads.query_sequence[relative_pos] == 'C':
                        df.loc[i, 'result_C'] += 1
                        if relative_pos < 4 or bam_reads.l_seq - relative_pos < 4:
                            df.loc[i, 'result_C_end'] += 1
                else:
                    print('ERROR!!!!!')

end = time.time()
print(end - start)

df_export = df.loc[df['result_A'] + df['result_T'] + df['result_G'] + df['result_C'] > 0].copy()

df_export['result_text'] = '0'
df_export['result_text'] = df_export[['result_A', 'result_T', 'result_G', 'result_C']].idxmax(axis=1)
df_export['result_text'] = df_export['result_text'].str.strip().str[-1]

#pd.options.mode.chained_assignment = None
df_export['result_final'] = 'new'
#df_export[df_export['result_text'] == df_export['ancestral']]['result_final'] = 'ancestral'
#df_export[df_export['result_text'] == df_export['derived']]['result_final'] = 'derived'
df_export.loc[df_export['result_text'] == df_export['ancestral'], ['result_final']] = 'ancestral'
df_export.loc[df_export['result_text'] == df_export['derived'], ['result_final']] = 'derived'

df_export.to_csv(output_file_name, encoding='utf-8', index=False)