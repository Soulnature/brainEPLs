import numpy as np
import os
import pandas as pd


l = []
for linkf in os.listdir('./data/eQTL_link/'):
    if not '.eQTL_link.txt' in linkf:
        continue

    links_df = pd.read_csv(os.path.join('./data/eQTL_link/', linkf), sep='\t', header=None,index_col=None)

    enhancers = links_df.iloc[:, 0].drop_duplicates()
    promoters = links_df.iloc[:, 1].drop_duplicates()

    r = len(links_df)/(len(enhancers) * len(promoters))
    l.append([linkf, r])
    print('tissue: %s, r: %.5f' % (linkf, r))

df = pd.DataFrame(l, columns=['tissue', 'sparsity']).sort_values(by='sparsity',ascending=False)

linkf = 'Kidney_Cortex.eQTL_link.txt'
links_df = pd.read_csv(os.path.join('./data/eQTL_link/', linkf), sep='\t', header=None,index_col=None)

enhancers = links_df.iloc[:, 0].drop_duplicates()
promoters = links_df.iloc[:, 1].drop_duplicates()
