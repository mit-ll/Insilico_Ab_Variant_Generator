'''
    __version__="1.0"
    __description__ = "Script to create 'k' mutations in the CDR regions of an antibody sequence"
    __copyright__= "© 2021 MASSACHUSETTS INSTITUTE OF TECHNOLOGY"
    __disclaimer__="THE SOFTWARE/FIRMWARE IS PROVIDED TO YOU ON AN “AS-IS” BASIS."
    __SPDX_License_Identifier__="BSD-2-Clause"

'''



import random, math
from collections import Counter, OrderedDict
from numpy.random import choice

import itertools
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from multiprocessing import Pool
from functools import partial
from scripts.parse_region import *
from copy import deepcopy
from Bio.SeqIO.FastaIO import FastaWriter

from Bio import pairwise2

def permute_seq(variable, p):

    '''Produces all combinations for input CDR region indices 
    and all permutations of amino acids for k mutations

    :param list variable: list of index positions for CDR regions of antibody sequence
    :param int p: k mutation
    :returns:
        - idx - if k=1 mutation, simply return indices (list of lists)
        - aa - if k=1 mutation, simply return amino acids (list of lists)
        - permute_idx - combinations of indices given k mutations
        - permute_aa - permutations of 20 amino acids given k mutations 

    '''

    if p == 1:
        idx = [[v] for v in variable]
        aa = [[a] for a in IUPAC_CODES.values()]
        return idx, aa
    permute_idx = [list(pi) for pi in list(itertools.combinations(variable, p))]
    permute_aa = [list(pa) for pa in list(itertools.product(IUPAC_CODES.values(), repeat=p))]


    return permute_idx, permute_aa

        

def replace_aa(sequence, p, idx_aa):
    '''Replaces amino acids at given indices in original antibody sequence to produce k "mutations"

    :param list sequence: antibody sequence in list form
    :param int p: k # mutations
    :param list idx_aa: list of index, amino acid pairs for k mutations (e.g. [[0,1],['A','C']])
    :returns:
        - chain - string sequence with amino acids replaced in index positions

    '''

    mapping = list(zip(idx_aa[0],idx_aa[1]))

    seq = deepcopy(sequence)

    for (i,j) in mapping:
        if p == 1 and seq[i]==j:
            return []
        else:
            seq[i] = j
    chain = ''.join(seq)
    return chain



def merge_lists(idx, aa):
    ''' Helper function to merge CDR index combinations with amino acid permutations if k=1 mutation
    :param list idx: list of index combinations
    :param list aa: list of amino acid permutations
    :returns:
        - list - list of combinations of indices and amino acids

    '''

    for i in idx:
        for a in aa:
            yield [i,a]
                
def merge_lists_w_sampling(aa, sample,sequence, diff, i):

    ''' Merges combinations of CDR region indices and permutations of amino acids for k mutations, 
    randomly sampled but equally distributed across CDR region indices

    :param list aa: list of amino acid permutations for k mutations
    :param int sample: max sample size
    :param list sequence: antibody sequence in list form
    :param set diff:  set containing random indices to determine if current index needs additional samples to make the required final sample size
    :returns:
        - ret - list of lists for CDR region index/indices to sampled amino acid permutations

    '''
    ret = [] 

    ###
    # Filter out amino acid permutations that have the same amino acid as the position we're trying to mutate  
    # (would result actually k-1 mutations of the current k)
    ###
    ignore = [sequence[ii] for ii in i[1]]

    reduce_aa = []

    for item in aa:

        test = list(zip(item, ignore))

        if not any(a == b for (a,b) in test):
            reduce_aa.append(item)
    mapping = list(enumerate(reduce_aa))
    
    # add on additional samples if need to get to required final sample size
    if diff and i[0] in diff:
        sample += 1

    #get random sample from reduced list of amino acid permutations
    s = set(random.sample(range(len(reduce_aa)), sample))

    #pair amino acid with CDR index position
    ret = [[i[1],a] for (c,a) in mapping if c in s]
    result_aa = [b for [a,b] in ret]


    return ret


  
   
        
def get_sample(idx, aa, sample, p, variable, sequence):
    '''Given sample size of total k mutations, gets random subset of CDR index/ amino acid combinations 
    and with equal distributions across indices. 
    
    For example, given list of indicies idx = [0,1,2] and p = 2, combinations of indices would be [0,1],[1,2] 
    and possible amino acids at those positions could be ['A','C'],['D','E']... etc. If we could only have 10 
    samples, then the script would randomly pick 5 combinations of amino acids for indices [0,1] and 5 
    combinations for indices [1,2].

    :param list idx: list of lists for all possible combinations of CDR region indices for k mutations
    :param list aa: list of lists for all possible permutations of amino acids for k mutations
    :param int sample: max sample size for k mutations
    :param int p: k mutation number
    :param list sequence: sequence in list form
    :returns:
        - res - list of index/amino acid combinations randomly sampled, but equally distributed across CDR regions

    '''

    if p > 1:    

        max_samp = int(sample/len(idx))
        proj_total = max_samp*len(idx)

        # make sure the projected total equals the required sample size after equally distributing random aa samples
        diff_samp = set()
        if proj_total < sample:
            diff = sample-proj_total
            diff_samp = set(random.sample(range(len(idx)), diff))
        enum_idx = [list(item) for item in list(enumerate(idx))]
        
        pool = Pool()
        func = partial(merge_lists_w_sampling, aa, max_samp, sequence, diff_samp)

        result_list = pool.map(func, enum_idx)
        pool.close()
        pool.join()
        res = [item for items in result_list for item in items]
        
        counts = Counter()
        for item in res:
            counts += Counter(item[0])
    
    else:
        res = list(merge_lists(idx, aa))
        counts = Counter()
        for item in res:
            counts += Counter(item[0])
        for k,v in counts.items():
            counts[k]-=1

    return res

def write_fasta(mutations, filename, heavy):
    ''' Writes sequences to fasta file.
    
    :param list mutations: list of mutated sequences (str)
    :param str filename: name of output file
    :param bool heavy: whether the input strings are for heavy chain or light chain
    :returns:
        - None

    '''
    
    if heavy:
        fname = filename+'-VH'
    else:
        fname = filename+'-VL'
        
    with open(fname+'.fa',"a") as handle:
        writer = FastaWriter(handle, wrap=None)
        c = 1
        for m in mutations:
            record1 = SeqRecord(Seq(m), id=fname+'-'+str(c))
            SeqIO.write(record1, handle, "fasta-2line")

            c+=1       

def parallelize(variable, sequence, p, filename, samples, heavy=False):
    '''Runs multiprocessing to produce k mutations and write to fast file

    :param list variable: list of indices in sequence for CDR regions
    :param list sequence: sequence in list form 
    :param int p: k mutation number
    :param str filename: name of output file
    :param int samples: number of samples from total k mutations to keep 
    :param bool heavy:  boolean to determine whether input sequence is heavy or light
    :returns:
        - None

    '''

    permute_idx, permute_aa = permute_seq(variable, p)
   
    sampling = get_sample(permute_idx, permute_aa, samples[p], p, variable, sequence)

    pool = Pool()
    func = partial(replace_aa,sequence,p)
    result_list = pool.map(func, sampling)
    res = [r for r in result_list if r ]
    pool.close()
    pool.join()
        
    
    write_fasta(res, filename, heavy)
    
    
    

def run_permute_seq(fname, heavychain, lightchain, heavy_mutations, light_mutations, save_dir):
    
    
    ''' Parse out CDR variable regions to be used in k mutations'''
    
    h, l,  listSeqHeavy, listSeqLight, regions, lightIdx = get_region_idx(heavychain, lightchain)
   

    ''' Last check to make sure n-mutations within range '''
    updated_heavy_mutations = {}
    for k,v in heavy_mutations.items():
        if k <= len(h) and v != -1:
            updated_heavy_mutations[k] = v
        elif k<= len(h) and v == -1:
            updated_heavy_mutations[k] = len(h)*20-len(h)
            print("Setting default for heavy chain: ", k, " mutations, ", len(h)*20-len(h), " samples.")
        elif k > len(h):
            print("Error: N-mutation is larger than total heavy chain CDR length. No variants generated for input: ", k)

    updated_light_mutations = {}
    for k,v in light_mutations.items():
        if k <= len(l) and v != -1:
            updated_light_mutations[k] = v
        elif k<= len(l) and v == -1:
            updated_light_mutations[k] = len(l)*20-len(l)
            print("Setting default for light chain: ", k, " mutations, ", len(h)*20-len(h), " samples.")
        elif k > len(l):
            print("Error: N-mutation is larger than total light chain CDR length. No variants generated for input:  ", k)

    ''' For each k mutation, run parallelize script to produce mutations and write to fasta files '''
    
    print("Producing sequence variants for: ", fname)

    for i in list(heavy_mutations.keys()):
        filename = save_dir+fname + '_heavy_mutation'+ str(i)
        if len(updated_heavy_mutations) > 0:
            
            parallelize(h, listSeqHeavy, i, filename, updated_heavy_mutations, heavy = True) 

        print("Variants stored in ", filename)
    for i in list(light_mutations.keys()):
        filename = save_dir + fname + '_light_mutation'+str(i)
        if len(updated_light_mutations) > 0:
            parallelize(l, listSeqLight, i, filename, updated_light_mutations, heavy = False)
        
        print("Variants stored in ", filename)
    
