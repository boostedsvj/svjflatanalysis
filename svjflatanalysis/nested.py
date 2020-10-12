"""
Legacy module to deal with old, nested ntuples from TreeMaker
"""

import svjflatanalysis

def filter_zerojet_events(arrays, inplace=True):
    passes = (arrays[b'JetsAK15'].counts >= 1)
    return select(arrays, passes, inplace)
