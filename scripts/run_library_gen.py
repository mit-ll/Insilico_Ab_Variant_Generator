'''
    __version__="1.0"
    __description__ = "Main script to initiate config/command line argument parsing"
    __copyright__= "© 2021 MASSACHUSETTS INSTITUTE OF TECHNOLOGY"
    __disclaimer__="THE SOFTWARE/FIRMWARE IS PROVIDED TO YOU ON AN “AS-IS” BASIS."
    __patent_ownership__="Subject to FAR 52.227-11 – Patent Rights – Ownership by the Contractor (May 2014)"
    __SPDX_License_Identifier__="BSD-2-Clause"

'''



import argparse
import json
from scripts.permute_seq import *

import os

def command_line():
    ''' List out sequence name, heavy and light chain sequences, and dictionary for k mutation : # of sequence samples to produce '''
    parser = argparse.ArgumentParser()
    parser.add_argument("config_file", type=str, help="Path to configuration file to load inputs.")
    parser.add_argument("--output_dir", type=str, help="Path to save output files.")

    args = parser.parse_args()

    inputs = json.load(open(args.config_file,'r'))
    save_dir = args.output_dir if args.output_dir[-1]=='/' else args.output_dir+'/'

    if not os.path.exists(save_dir):
        os.mkdir(save_dir)

    print("Parsing config file")
    print("###################")

    for item in inputs['inputs']:
        f = item['filename']
        heavychain = item['heavychain']
        lightchain = item['lightchain']
        heavy_mutations = {}
        len_heavy = len(heavychain)
        print("Parsing heavy chain mutation/samples: ", f)
        print("#######################################")
        print()
        for hm in item['heavy_mutations']:
            if len(hm) == 0 and heavychain != '':
                print('Error, n-mutations and variant count not specified. Empty list')
                continue
            elif len(hm) >2:
                print('Error, too many values to unpack: ', len(hm))
                continue
            elif len(hm) == 1:
                heavy_mutations[hm[0]] = -1
            elif len(hm) == 2:
                if hm[0] > 6:
                    print("Num mutations for heavy chain set too high. Will cause performance issue. Defaulting to 6")
                    heavy_mutations[6] = hm[1]
                else:
                    heavy_mutations[hm[0]] = hm[1]
        print()
        print("Parsing light chain mutation/samples: ", f)
        print("#######################################")
        print()
        
        light_mutations = {}
        len_light = len(lightchain)
        for lm in item['light_mutations']:
            if len(hm) == 0 and lightchain != '':
                print('Error, n-mutations and variant count not specified. Empty list')
                continue
            elif len(hm) >2:
                print('Error, too many values to unpack: ', len(hm))
                continue
            elif len(hm) == 1:
                light_mutations[lm[0]] = -1
            elif len(lm) == 2:
                if lm[0] > 6:
                    print("Num mutations for light chain set too high. Will cause performance issue. Defaulting to 6.")
                    light_mutations[6] = lm[1]
                else:
                    light_mutations[lm[0]] = lm[1]
        
        run_permute_seq(f, heavychain, lightchain, heavy_mutations, light_mutations, save_dir)

    print("Antibody library generation complete")

