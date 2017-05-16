
# coding: utf-8

# In[ ]:

import numpy as np
from collections import Counter
from Bio import SeqIO

'''Generate the word-context count'''
def collect_pairs(sentence,window,pair_counter):
    for ind_1 in range(len(sentence)):
        word_1 = sentence[ind_1]
        for ind_2 in range(ind_1+1,ind_1+window+1):
            if ind_2 < len(sentence): 
                word_2 = sentence[ind_2]
                pair_counter[(word_1,word_2)] += 1
                pair_counter[(word_2,word_1)] += 1 # symmetrize the counting
    return pair_counter

'''Given an exome creates the associated sentence'''
def exome2sentence(exome,minimo,massimo):
    sentence = []
    for start in range(len(exome)):
        word = str(exome[start:start+np.random.random_integers(minimo,massimo)]) # random word-length from 3 to 9
        if len(word) >= 3: sentence.append(word) # to avoid last short words
    return sentence


# In[ ]:

pair_counter = Counter()
window = 10
minimo = 3
massimo = 9
num_cycles = 50

with open("/home/garner1/Work/dataset/gene2vec/genedoc/896.fasta", "rU") as handle:
    for record in SeqIO.parse(handle, "fasta"): # loop over different exomes
        exome = record.seq
        for cycle in range(num_cycles): # generate different random contexts for each exome
            sentence = exome2sentence(exome,minimo,massimo)
            pair_counter = collect_pairs(sentence,window,pair_counter)
#     print pair_counter.most_common(4)


# In[ ]:

import pandas as pd

df = pd.DataFrame.from_dict(pair_counter, orient='index').reset_index()
df = df.rename(columns={'index':'pair', 0:'count'})
df.sort_values('count',ascending=False,inplace=True)

new_col_list = ['word_1','word_2']
for n,col in enumerate(new_col_list):
    df[col] = df['pair'].apply(lambda pair: pair[n])
df = df.drop('pair',axis=1)


# In[ ]:

d = {ni: indi for indi, ni in enumerate(set(df['word_1']))} # create a dictionary between unique words and numbers
df['ind_1'] = df['word_1'].map(d) # assign the dict to word_1
df['ind_2'] = df['word_2'].map(d) # assign the dict to word_2


# In[ ]:

df_grouped = df.groupby(['word_1']).sum() # group by word (row) and sum along columns
d = df_grouped.to_dict()['count'] # create a dictionary between word and their col-sum
df['ind_1_sumOver2'] = df['word_1'].map(d) # add the norm1(row) field

df_grouped = df.groupby(['word_2']).sum() # group by word (col) and sum along rows
d = df_grouped.to_dict()['count'] # create a dictionary between word and their row-sum
df['ind_2_sumOver1'] = df['word_2'].map(d) # add the norm1(col) field


# In[ ]:

import scipy.sparse as sps
card_of_D = df.shape[0]
dim = len(set(df['word_1']))
mat = sps.csr_matrix((np.log(df['count']*card_of_D*1.0/(df['ind_2_sumOver1']*df['ind_1_sumOver2'])), (df['ind_1'], df['ind_2'])),shape=(dim,dim)) # convert the pair count into a sparse matrix


# In[ ]:

from scipy.sparse.linalg import svds
U, S, UT = svds(mat, k = 2)

import matplotlib.pyplot as plt
plt.figure()

# plt.scatter(range(len(S)),S)
num_points = 1000
# plt.scatter(U[:num_points,0],U[:num_points,1])
plt.scatter(U[:,0],U[:,1])

plt.show()


# In[ ]:

from scipy.sparse.linalg import svds
U, S, UT = svds(mat, k = 3)

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

ax.scatter(U[:,0],U[:,1],U[:,2])
plt.show()


# In[ ]:



