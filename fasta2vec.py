
# coding: utf-8
import numpy as np
from collections import Counter
from Bio import SeqIO
import pandas as pd
import scipy.sparse as sps
from scipy.sparse.linalg import svds, norm
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from itertools import product
from sklearn.metrics.pairwise import cosine_similarity
import sys
import pickle

'''Generate and add to the pair-counter the new word-context counts'''
def collect_pairs(sentence,window,pair_counter):
    for ind_1 in range(len(sentence)):
        word_1 = sentence[ind_1]
        for ind_2 in range(ind_1+1,ind_1+window+1):
            if ind_2 < len(sentence): 
                word_2 = sentence[ind_2]
                pair_counter[(word_1,word_2)] += 1
                pair_counter[(word_2,word_1)] += 1 # symmetrize the counting
    return pair_counter

'''Given an exon creates the associated sentence'''
def exon2sentence(exon,minimo,massimo):
    sentence = []
    for start in range(len(exon)):
        word = str(exon[start:start+np.random.random_integers(minimo,massimo)]) # random word-length from 3 to 9
        if len(word) >= 3: sentence.append(word) # to avoid last short words
    return sentence

'''Given a counter produces a PPMI matrix'''
def counter2mat(pair_counter_in_sentence,dic,dim):
    df = pd.DataFrame.from_dict(pair_counter_in_sentence, orient='index').reset_index()
    df = df.rename(columns={'index':'pair', 0:'count'})
    df.sort_values('count',ascending=False,inplace=True)
    new_col_list = ['word_1','word_2'] # split pair into 2 columns
    for n,col in enumerate(new_col_list):
        df[col] = df['pair'].apply(lambda pair: pair[n])
    df = df.drop('pair',axis=1)
    df['ind_1'] = df['word_1'].map(d) # assign the dict to word_1
    df['ind_2'] = df['word_2'].map(d) # assign the dict to word_2
    df_grouped = df.groupby(['word_1']).sum() # group by word (row) and sum along columns
    d2 = df_grouped.to_dict()['count'] # create a dictionary between word and their col-sum
    df['ind_1_sumOver2'] = df['word_1'].map(d2) # add the norm1(row) field
    df_grouped = df.groupby(['word_2']).sum() # group by word (col) and sum along rows
    d2 = df_grouped.to_dict()['count'] # create a dictionary between word and their row-sum
    df['ind_2_sumOver1'] = df['word_2'].map(d2) # add the norm1(col) field

    card_of_D = df.shape[0]
    mat = sps.csr_matrix((np.log(df['count']*card_of_D*1.0/(df['ind_2_sumOver1']*df['ind_1_sumOver2'])), (df['ind_1'], df['ind_2'])),shape=(dim,dim)) # convert the pair count into a sparse matrix
    return mat

'''
Generate the vocabulary, i.e. the collections of possible kmers and assign them an index for the PPMI matrix
'''
bases = ['A','T','G','C']
kmers = []
for k in np.arange(3,10):
    kmers.extend([''.join(p) for p in product(bases, repeat=k)])
'''
Create a dictionary between unique words and numbers
'''
d = {ni: indi for indi, ni in enumerate(set(kmers))} # 
dim = len(set(kmers)) # dimension of the sparse matrix
'''
Set the fasta file to analyze
'''
file_name = sys.argv[1] # ex: "/home/garner1/Work/dataset/gene2vec/genedoc/27.fasta"
'''
Iniziatize the counter
'''
pair_counter_in_gene = Counter()
'''
Determine the window of the context, i.e. # of words in the neighborhood, and the min/max word-length
'''
window = int(sys.argv[2])            # ex:10
minimo = int(sys.argv[3])            # ex:3
massimo = int(sys.argv[4])           # ex:9
'''
Process the fasta file
'''
with open(file_name, "rU") as handle:
    '''Loop over different exons'''
    for record in SeqIO.parse(handle, "fasta"):
        pair_counter_in_sentence = Counter()
        exon = record.seq
        cycle = 0
        condition = True
        while condition :
            cycle += 1
            '''Generate a random sentence from the exon'''
            sentence = exon2sentence(exon,minimo,massimo)
            '''Add new word co-occurrences inside the context to the old pair_counter'''
            pair_counter_in_sentence = collect_pairs(sentence,window,pair_counter_in_sentence)
            
            # print '{0}\r'.format(str(cycle)+' of '+str(len(exon))),
            '''Check every len(exon) new sentences if the stop condition is satisfied'''
            if cycle%len(exon) == 0: 
                '''Generate the PPMI matrix'''
                mat1 = counter2mat(pair_counter_in_sentence,d,dim)
                '''Generate a perturbation to the PPMI matrix via a new random sentence'''
                sentence = exon2sentence(exon,minimo,massimo)
                pair_counter_in_sentence = collect_pairs(sentence,window,pair_counter_in_sentence)
                mat2 = counter2mat(pair_counter_in_sentence,d,dim)
                '''Estimate the relative variation in the norm of the diff of the similarity matrices'''
                pdist_1 = cosine_similarity(mat1,dense_output=False)
                pdist_2 = cosine_similarity(mat2,dense_output=False)
                variation = 100*norm(pdist_1-pdist_2,'fro')/norm(pdist_1,'fro')
                condition = variation >= 0.1
                '''If the condition is satisfied add the pair_counter to the gene pair_counter'''
                if ~condition: 
                    # print '{0}\r'.format('exon length:' +str(len(exon))+'; variation: '+str(variation)+'; loop #: '+str(cycle))
                    pair_counter_in_gene += pair_counter_in_sentence
'''
Write the pair counter of the gene in a file
'''
with open('/media/bicroserver-seq/gene2vec/'+str(file_name.split('/')[-1])+ '.pickle', 'wb') as outputfile:
    pickle.dump(pair_counter_in_gene, outputfile)

# '''
# Load the pair counter from file
# '''
# with open(file_name + '.pickle', 'rb') as inputfile:
#     print pickle.load(inputfile)
