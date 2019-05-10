#This code creates a file with SNP list in the required format from a file, copied fom YSOGG
#Note: you need to delete from YSOGG file info about all mutations that are not SNP before running this script
import pandas as pd

path_import = 'project_data/R_happl.csv'
path_export = 'project_data/R_fixed.csv'

df = pd.read_csv(path_import)

#df = df.drop(['Unnamed: 7'], axis=1)
df['ancestral'] = df['Mutation info'].astype(str).str[0]
df['derived'] = df['Mutation info'].str.strip().str[-1]
df = df.drop(['Mutation info'], axis=1)
df = df.rename(index=str, columns={"Build 37 #": "build_37", "Build 38 #": "build_38"})
df['build_37'] = df['build_37'].fillna(0)
df = df.dropna(subset={'build_38'})
df['build_37'] = df['build_37'].astype(int)
df['build_38'] = df['build_38'].astype(int)
df = df.drop_duplicates(subset={'Haplogroup', 'build_38'})
df.to_csv(path_export, encoding='utf-8', index=False)