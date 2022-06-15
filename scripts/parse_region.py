'''
    __version__="1.0"
    __description__ = "Script to parse CDR regions in both heavy and light chains using the 
                       Martin numbering scheme: http://www.bioinf.org.uk/abs/info.html#martinnum"
    __copyright__= "© 2021 MASSACHUSETTS INSTITUTE OF TECHNOLOGY"

    __disclaimer__="THE SOFTWARE/FIRMWARE IS PROVIDED TO YOU ON AN “AS-IS” BASIS."

    __SPDX_License_Identifier__="BSD-2-Clause"
'''

from collections import OrderedDict
import math

# Dictionary of 20 standard amino acids encoding 3-letter to 1-letter codes
IUPAC_CODES = OrderedDict([
    ('Ala', 'A'),
    ('Cys', 'C'),
    ('Asp', 'D'),
    ('Glu', 'E'),
    ('Phe', 'F'),
    ('Gly', 'G'),
    ('His', 'H'),
    ('Ile', 'I'),
    ('Lys', 'K'),
    ('Leu', 'L'),
    ('Met', 'M'),
    ('Asn', 'N'),
    ('Pro', 'P'),
    ('Gln', 'Q'),
    ('Arg', 'R'),
    ('Ser', 'S'),
    ('Thr', 'T'),
    ('Val', 'V'),
    ('Trp', 'W'),
    ('Tyr', 'Y')])

def find_aa(aa_list, aa, idx1, idx2, range=False):
    if range:
        target = [i for (i,a) in enumerate(aa_list) if a==aa and (abs(idx1-i)<10 or abs(idx2-i)<10)]
    else:
        target = [i for (i,a) in enumerate(aa_list) if a==aa and abs(idx1-i)<10]
    return target

def get_cdrl1(listSeq):
    ''' Parses indices in sequence for CDRL1 region

    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRL1 region

    '''
    start_target = find_aa(listSeq, 'C', 23, -1, range=False)
    
    start = math.inf
    for s in start_target:
        if s < start:
            start = s
    start += 1
    end_target = find_aa(listSeq, 'W', start+9, start+16, range=True)
    
    end_aa = [['W','Y','Q'],['W','L','Q'],['W','F','Q'],['W','Y','L']]
    end = -1
    for idx in end_target:
        if listSeq[idx:idx+3] in end_aa:
            end = idx
            break
    end -=1
    assert len(listSeq[start:end+1])>=10 and len(listSeq[start:end+1]) <=17
    return [start, end]
        

def get_cdrl2(cdrl1_end, listSeq):
    ''' Parses indices in sequence for CDRL2 region
    
    :param int cdrl1_end: end index of CDRL1 (i.e. beginning of CDRL2)
    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRL2 region

    '''
    start = cdrl1_end + 16
    
    before_seq = [['I','Y'],['V','Y'],['I','K'],['I','F']]
    if not listSeq[start-2:start] in before_seq:
        end_target = find_aa(listSeq, 'C', start+40, -1, range=False)
        end = -1
        for idx in end_target:
            if listSeq[idx] == 'C':
                end = idx
                break
        end -= 32
        assert listSeq[end+32] == 'C'
        return [start, end]
    else:
        end = start+6
        
        assert len(listSeq[start:end+1])==7
        assert listSeq[end+32] == 'C'
        return [start, end]

def get_cdrl3(cdrl2_end, listSeq):
    ''' Parses indices in sequence for CDRL3 region
    
    :param int cdrl1_end: end index of CDRL2 (i.e. beginning of CDRL3)
    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRL3 region

    '''


    start = cdrl2_end + 33
    assert listSeq[start-1] == 'C'

    end_target = find_aa(listSeq, 'F', start+6, start+10, range=True)
    end = -1
    for idx in end_target:
        if listSeq[idx:idx+2] == ['F','G'] and listSeq[idx+3] == 'G':
            end = idx
            break
    end -= 1
    assert len(listSeq[start:end+1])>=7 and len(listSeq[start:end+1]) <=11
    return [start,end]


def get_cdrh1(listSeq):
    ''' Parses indices in sequence for CDRH1 region
    
    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRH1 region

    '''

    start_target = find_aa(listSeq, 'C', 26, -1, range=False)
    start = -1
    for idx in start_target:
        if listSeq[idx] == 'C':
            start = idx
            break
    start += 4
    
    end_target = find_aa(listSeq, 'W', start+9, start+11, range=True)
    after_seq = [['W','V'],['W','I'],['W','A']]
    end = -1
    for idx in end_target:
        if listSeq[idx:idx+2] in after_seq:
            end = idx
            break
    end -=1
    assert len(listSeq[start:end+1])>=10 and len(listSeq[start:end+1]) <=12
    
    return [start, end]                           

def get_cdrh2(cdrh1_end, listSeq):
    ''' Parses indices in sequence for CDRH2 region
    
    :param int cdrl1_end: end index of CDRH1 (i.e. beginning of CDRH2)
    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRH1 region

    '''

    start = cdrh1_end + 15
    
    end_target = find_aa(listSeq, 'C', start+42, start+45, range=True)
    end = -1
    for idx in end_target:
        if listSeq[idx] == 'C':
            end = idx
            break
    end -= 31
    assert len(listSeq[start:end+1])>=16 and len(listSeq[start:end+1]) <=19
   
    return [start,end]

def get_cdrh3(cdrh2_end, listSeq):
    ''' Parses indices in sequence for CDRH3 region
    
    :param int cdrl1_end: end index of CDRH2 (i.e. beginning of CDRH3)
    :param list listSeq: antibody sequence in list form
    :returns:
        - list - start and end indices for CDRH3 region

    '''

    start = cdrh2_end + 34
    assert listSeq[start-3] == 'C'
    end_target = find_aa(listSeq, 'W', start+ 3, start+25, range=True)
    end = -1
    for idx in end_target:
        if listSeq[idx:idx+2] == ['W','G'] and listSeq[idx+3] == 'G':
            end = idx
            break
    end -= 1
    assert len(listSeq[start:end+1])>=3 and len(listSeq[start:end+1]) <=25
    return [start, end]

def get_region_idx(heavy, light):
    '''Parses CDR regions of antibody sequence using indices

    :param str heavy: sequence for heavy chain 
    :param str light: sequence for light chain
    :returns:
        - heavy_idxs - list of indices for heavy CDR regions
        - light_idxs - list of indices for light CDR regions
        - listSeqHeavy - sequence in list form
        - listSeqLight - sequence in list form
        - regions - dictionary with CDR regions (key) and start/end indices list (key)
        - lightIdx - index of start of light chain if you combine heavy/light chain into one string (not being used currently)



    '''
    listSeqHeavy = list(heavy)
    listSeqLight = list(light)
    lightIdx = len(listSeqHeavy)
    regions = {}
   
    regions['cdrh1'] = get_cdrh1(listSeqHeavy)
    regions['cdrh2'] = get_cdrh2(regions['cdrh1'][1], listSeqHeavy)
    regions['cdrh3'] = get_cdrh3(regions['cdrh2'][1], listSeqHeavy)
    
    regions['cdrl1'] = get_cdrl1(listSeqLight)
    regions['cdrl2'] = get_cdrl2(regions['cdrl1'][1], listSeqLight)
    regions['cdrl3'] = get_cdrl3(regions['cdrl2'][1], listSeqLight)
    
    heavy_idxs = list(range(regions['cdrh1'][0],regions['cdrh1'][1]+1))+\
               list(range(regions['cdrh2'][0],regions['cdrh2'][1]+1))+\
               list(range(regions['cdrh3'][0],regions['cdrh3'][1]+1))
    light_idxs = list(range(regions['cdrl1'][0],regions['cdrl1'][1]+1))+\
               list(range(regions['cdrl2'][0],regions['cdrl2'][1]+1))+\
               list(range(regions['cdrl3'][0],regions['cdrl3'][1]+1))
    #print('heavy: ',len(heavy))
    #print('light: ',len(light))
               
    print("cdrl1: ", ''.join(listSeqLight[regions['cdrl1'][0]:regions['cdrl1'][1]+1]))
    print("cdrl2: ", ''.join(listSeqLight[regions['cdrl2'][0]:regions['cdrl2'][1]+1]))
    print("cdrl3: ", ''.join(listSeqLight[regions['cdrl3'][0]:regions['cdrl3'][1]+1]))
    print("cdrh1: ", ''.join(listSeqHeavy[regions['cdrh1'][0]:regions['cdrh1'][1]+1]))
    print("cdrh2: ", ''.join(listSeqHeavy[regions['cdrh2'][0]:regions['cdrh2'][1]+1]))
    print("cdrh3: ", ''.join(listSeqHeavy[regions['cdrh3'][0]:regions['cdrh3'][1]+1]))
    return heavy_idxs, light_idxs, listSeqHeavy, listSeqLight, regions, lightIdx

            
