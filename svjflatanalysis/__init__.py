import os.path as osp, logging

DEFAULT_LOGGING_LEVEL = logging.DEBUG
def setup_logger(name='datasets'):
    if name in logging.Logger.manager.loggerDict:
        logger = logging.getLogger(name)
        logger.info('Logger %s is already defined', name)
    else:
        fmt = logging.Formatter(
            fmt = (
                '\033[33m%(levelname)7s:%(asctime)s:%(module)s:%(lineno)s\033[0m'
                + ' %(message)s'
                ),
            datefmt='%Y-%m-%d %H:%M:%S'
            )
        handler = logging.StreamHandler()
        handler.setFormatter(fmt)
        logger = logging.getLogger(name)
        logger.setLevel(DEFAULT_LOGGING_LEVEL)
        logger.addHandler(handler)
    return logger
logger = setup_logger()

def is_interactive():
    import __main__ as main
    return not hasattr(main, '__file__')
IS_INTERACTIVE = is_interactive()

try:
    if IS_INTERACTIVE:
        from tqdm.notebook import tqdm
        logger.info('Using tqdm notebook')
    else:
        from tqdm import tqdm
    HAS_TQDM = True
except ImportError:
    logger.error('Could not import tqdm')
    HAS_TQDM = False

from . import utils
from . import arrayutils
from . import crosssections
from . import dataset
from .dataset import (
    Dataset,
    add_to_bytestring,
    numentries,
    iterate_events,
    iterate,
    iterate_events_datasets,
    )
from . import samples
from . import roccurve
from . import trigger

# Default plotting params
try:
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 18})
except ImportError:
    logger.error('Could not import matplotlib')
try:
    import mplhep as hep
    plt.style.use(hep.style.CMS)
except ImportError:
    logger.error('Could not import mplhep')

TRIGGER_TITLES = None
def get_trigger_titles(example_rootfile=None):
    """
    Reads triggers from any treemaker root file
    Caches result in global var for any calls after the first one
    """
    import uproot
    if example_rootfile is None:
        example_rootfile = samples.get_bkg_rootfiles()[0]
    global TRIGGER_TITLES
    if TRIGGER_TITLES is None:
        title_bstring = uproot.open(example_rootfile).get('TreeMaker2/PreSelection')[ b'TriggerPass'].title
        TRIGGER_TITLES = title_bstring.decode('utf-8').split(',')
    return TRIGGER_TITLES
