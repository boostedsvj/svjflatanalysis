from time import strftime
import os.path as osp, os, re
import svjflatanalysis
logger = svjflatanalysis.logger

def get_ax(**kwargs):
    if 'ax' in kwargs and not(kwargs['ax'] is None):
        return kwargs['ax']
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
def get_default_color():
    """
    Returns one of the default pyplot colors. Every call returns a new color,
    cycling back at some point (for most versions 10 colors)
    """
    global _i_call_default_color
    import matplotlib.pyplot as plt
    default_colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    color = default_colors[_i_call_default_color]
    _i_call_default_color += 1
    if _i_call_default_color == len(default_colors):
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

def xkcd_color_wheel(seed=44):
    import matplotlib._color_data, random
    colors = list(sorted(matplotlib._color_data.XKCD_COLORS.keys()))
    random.seed(seed)
    random.shuffle(colors)
    return colors

VARTITLES = {
    'mt' : r'$m_{T}$',
    'met' : 'MET',
    'rt' : r'$R_{T}$',
    'dphimet' : r'$\Delta\phi_{\mathrm{MET}}$',
    'deltaeta' : r'$\Delta\eta_{\mathrm{l1,l2}}$',
    'axismajor' : 'Axis major',
    'axisminor' : 'Axis minor',
    'girth' : 'Girth',
    'ptd' : r'$p_{TD}$',
    'tau21' : r'$\tau_{21}$',
    'tau31' : r'$\tau_{31}$',
    'tau32' : r'$\tau_{32}$',
    }

def ecf_title_formatter(ecf):
    match = re.search(r'ecfn(\d)b(\d)', ecf)
    if match:
        n = match.group(1)
        beta = match.group(2)
        title = r'ECF N{}$\beta${}'.format(n, beta)
        return title
    return None    

def get_title(varname):
    varname = varname.lower().strip()
    if varname in VARTITLES:
        return VARTITLES[varname]
    title = ecf_title_formatter(varname)
    if not(title is None):
        return title
    logger.info('Could not find a nice title for varname %s', varname)
    return varname

