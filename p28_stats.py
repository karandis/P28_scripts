import numpy as np
import pandas as pd
from pathlib import Path

num_corr = len(list(Path.cwd().glob('*/corrected/*_cr.fits')))

df = pd.read_csv('prog28_good_20190702.csv',header=0,sep=',')

df.columns = [x.strip() for x in df.columns]
print('NEW',df.columns)
# for i,row in df.iterrows():
#     night = int(row['filename'].split("/")[4])
#     # object = " ".join(row['  starname  '].strip()).replace('*','_')
#     object = row['hom-starname'].strip().replace(' ','_').replace('*','_')
#     unseq = int(row['unseq'])
nights = [int(x.split("/")[4]) for x in df.filename]
df=df.assign(night=nights)

#filter for nights after 2011
df = df.loc[df.night>20111231]

count=0
list_responses = list(Path('/STER/karansinghd/PhD/ResponseCorrection/responses/').glob('*.txt'))
list_responses = [str(p) for p in list_responses]
#print(list_responses)
for n in df.night.values:
	list_p2 = [x for x in list_responses if str(n) in x]
	print(list_p2)
	exit()
	if str(n) in str(j.stem):
		count+=1
		break

print(len(df),num_corr,count)
