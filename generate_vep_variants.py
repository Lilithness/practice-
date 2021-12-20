import pandas as pd
from corvus.cmd import get_cmd_output
import vcf
import io
import os
import requests
import sys
import json
import pprint

def read_vcf(path):
    """
    Read a VCF file of alligned sequences
    :param path: VCF file path
    :return: Data frame of all data    
    """
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('##')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': int, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})
