#This code checks 23andme result files for known Y-SNPS

import numpy as np
import pandas as pd

#bam_file_path = 'project_data/GA231_GA247_merged.bam'
thandme_file_path = 'project_data/Artem_Y_23andme.csv'
target_snps_path = 'project_data/trunk_and_i_fixed.csv'
output_file_name = 'project_data/artem_results2.csv'

df = pd.read_csv(target_snps_path)
df_known_data = pd.read_csv(thandme_file_path)
df['result'] = 'no call'

pd.options.mode.chained_assignment = None

for i in range(0, len(df)):
    print(i)
    for j in range(0, len(df_known_data)):
        if df['build_37'][i] == df_known_data['position'][j]:
            df['result'][i] = df_known_data['genotype'][j]

df_export = df[df['result'] != 'no call']
df_export.to_csv(output_file_name, encoding='utf-8', index=False)
