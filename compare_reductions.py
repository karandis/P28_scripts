'''
Script to compare the reductions from April 2021 (P28_c) and Jan 2022 (P28_cr_2022)
'''
import numpy as np
import pandas as pd
from pathlib import Path
import matplotlib.pyplot as pl

old_std = Path('/STER/karansinghd/PhD/Projects/P28_c/stdinfo.csv')
new_std = Path('/STER/karansinghd/P28_cr_2022/stdinfo.csv')

old =  pd.read_csv(old_std)
new = pd.read_csv(new_std)

old_std1 = old.loc[old.night == old.stdnight]
old_std0 = old.loc[old.night != old.stdnight]

new_std1 = new.loc[new.night == new.stdnight]
new_std0 = new.loc[new.night != new.stdnight]

old_Ftypes = old.loc[old.stdname.isin(['HD185395','HD220657','HD206826'])]

modified = new.loc[new.unseq.isin(old_Ftypes.unseq)]
print(modified.stdname.value_counts())
print(modified.loc[modified.unseq.isin(old_Ftypes.loc[old_Ftypes.night == old_Ftypes.stdnight].unseq)])
# print(old_Ftypes.loc[old_Ftypes.night != old_Ftypes.stdnight])
save = modified.loc[modified.unseq.isin(old_Ftypes.loc[old_Ftypes.night == old_Ftypes.stdnight].unseq)]

save.to_csv('/STER/karansinghd/PhD/Projects/P28_scripts/modified.csv',index=False)
