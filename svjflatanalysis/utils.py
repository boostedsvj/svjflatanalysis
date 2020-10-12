import matplotlib.pyplot as plt
from time import strftime
import os.path as osp, os
import flatanalysis
logger = flatanalysis.logger

def get_ax(**kwargs):
    import matplotlib.pyplot as plt
    if not 'figsize' in kwargs: kwargs['figsize'] = (8,8)
    fig = plt.figure(**kwargs)
    ax = fig.gca()
    return ax

def is_string(string):
    """
    Checks strictly whether `string` is a string
    Python 2/3 compatibility (https://stackoverflow.com/a/22679982/9209944)
    """
    try:
        basestring
    except NameError:
        basestring = str
    return isinstance(string, basestring)

def bytes_to_human_readable(num, suffix='B'):
    """
    Convert number of bytes to a human readable string
    """
    for unit in ['','k','M','G','T','P','E','Z']:
        if abs(num) < 1024.0:
            return '{0:3.1f} {1}b'.format(num, unit)
        num /= 1024.0
    return '{0:3.1f} {1}b'.format(num, 'Y')

_i_call_default_color = 0
DEFAULT_COLORS = plt.rcParams['axes.prop_cycle'].by_key()['color']
def get_default_color():
    """
    Returns one of the default pyplot colors. Every call returns a new color,
    cycling back at some point (for most versions 10 colors)
    """
    global _i_call_default_color
    color = DEFAULT_COLORS[_i_call_default_color]
    _i_call_default_color += 1
    if _i_call_default_color == len(DEFAULT_COLORS):
        _i_call_default_color = 0
    return color

def reset_color():
    global _i_call_default_color
    _i_call_default_color = 0

def make_plotdir(path):
    path = path.replace('%d', strftime('%b%d'))
    path = osp.abspath(path)
    if not osp.isdir(path):
        logger.info('Creating %s', path)
        os.makedirs(path)
    return path