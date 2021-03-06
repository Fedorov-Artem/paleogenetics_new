#This code creates a file with SNP list in the required format from a file, copied fom YSOGG
import pandas as pd

path_import = 'project_data/J_new.csv'
path_export = 'project_data/J_new_fixed.csv'

df = pd.read_csv(path_import)

df = df.loc[df['Mutation info'].str.len() == 4]
df['ancestral'] = df['Mutation info'].astype(str).str[0]
df['derived'] = df['Mutation info'].str.strip().str[-1]
df = df.drop(['Mutation info'], axis=1)
df = df.rename(index=str, columns={"Build 37 #": "build_37", "Build 38 #": "build_38"})
df['build_37'] = pd.to_numeric(df['build_37'], errors='coerce', downcast='integer').fillna(0)
df['build_38'] = pd.to_numeric(df['build_38'], errors='coerce', downcast='integer').fillna(0)
df = df.loc[(df['build_38'] > 0) & (df['build_37'] > 0)]
df = df.drop_duplicates(subset={'Haplogroup', 'build_38'})
df.to_csv(path_export, encoding='utf-8', index=False)