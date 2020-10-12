# import uproot4 as uproot
import uproot, re
import flatanalysis
logger = flatanalysis.logger
tqdm = flatanalysis.tqdm

# ______________________________________________________________________
# arrays-level helpers

def add_to_bytestring(bytestring, tag):
    normal_string = bytestring.decode('utf-8')
    normal_string += tag
    return normal_string.encode('utf-8')

def numentries(arrays):
    return arrays[list(arrays.keys())[0]].shape[0]
    
def iterate_events(arrays):
    """
    Iterates event by event from an arrays of events
    """
    n = numentries(arrays)
    for i in range(n):
        yield { k : v[i:i+1] for k, v in arrays.items() }

# ______________________________________________________________________
# Operators for list of datasets

def iterate(datasets, **kwargs):
    """
    Yields arrays and the dataset from a list of datasets
    """
    # Convert to iterable if not an iterable
    if hasattr(datasets, 'treename'): datasets = [datasets]
    for dataset in datasets:
        for arrays in dataset.iterate(**kwargs):
            yield arrays, dataset

def iterate_events_datasets(datasets, **kwargs):
    """
    Yields arrays and the dataset from a list of datasets
    """
    for arrays, dataset in iterate(datasets, **kwargs):
        for event in iterate_events(arrays):
            yield event, dataset

# ______________________________________________________________________
# Dataset classes

class Dataset(object):
    """
    Container for a bunch of root files with an iterate method to easily do event loops
    over an arbitrary number of files.
    Has a cache functionality to store events in memory, making relooping very fast.
    """

    treename = 'TreeMaker2/PreSelection'
    
    def __init__(self, rootfiles, make_cache=True, **kwargs):
        # logger.debug('Initializing dataset with rootfiles: %s', rootfiles)
        super().__init__()
        self.rootfiles = [rootfiles] if flatanalysis.utils.is_string(rootfiles) else rootfiles
        if len(self.rootfiles) == 0:
            logger.warning('No rootfiles for %s; iterators will be empty', self)
            make_cache = False
        self.cache = []
        if make_cache:
            self.make_cache(**kwargs)

    def __repr__(self):
        return super().__repr__().replace('Dataset', 'Dataset ({0} root files)'.format(len(self.rootfiles)))

    def iterate(self, progressbar=False, n_files=None, use_cache=True, **kwargs):
        """
        Wrapper around uproot.iterate:
        - Gets a progress bar option
        - Possibility to limit number of files
        - Can use a class cache variable
        """
        if use_cache:
            if not len(self.cache):
                logger.warning('use_cache was True but cache is empty for %s', self)
            # logger.debug('Using cache')
            iterator = iter(self.cache)
            total = len(self.cache)
        else:
            # Allow reading only the first n_files root files
            rootfiles = self.rootfiles[:]
            if n_files: rootfiles = rootfiles[:n_files]
            # rootfiles = [ r + ':' + self.treename for r in rootfiles ]
            iterator = uproot.iterate(rootfiles, self.treename, **kwargs)
            total = len(rootfiles)
            if progressbar: logger.info('Iterating over %s rootfiles for %s', total, self)
        if progressbar:
            iterator = tqdm(iterator, total=total, desc='arrays' if use_cache else 'root files')
        for arrays in iterator:
            yield arrays
            
    def iterate_events(self, **kwargs):
        """
        Like self.iterate(), but yields a single event per iteration
        """
        for arrays in self.iterate(**kwargs):
            for event in iterate_events(arrays):
                yield event
            
    def make_cache(self, max_entries=None, **kwargs):
        """
        Stores result of self.iterate in a class variable for fast reuse.
        If max_entries is set, it will fill the cache until max_entries is exceeded
        for the first time or until the iterator runs out of files
        """
        if len(self.cache): logger.info('Overwriting cache for %s', self)
        self.cache = []
        self.sizeof_cache = 0
        self.numentries_cache = 0
        branches = None
        for arrays in self.iterate(use_cache=False, **kwargs):
            if branches is None: branches = list(arrays.keys())
            numentries = flatanalysis.arrayutils.numentries(arrays)
            if max_entries and self.numentries_cache + numentries > max_entries:
                # Cut away some entries to get to max_entries
                needed_entries = max_entries - self.numentries_cache
                arrays = { k : v[:needed_entries] for k, v in arrays.items() }
            self.cache.append(arrays)
            self.sizeof_cache += sum([ v.nbytes for v in arrays.values() ])
            self.numentries_cache += flatanalysis.arrayutils.numentries(arrays)
            if max_entries and self.numentries_cache >= max_entries:
                break
        logger.info(
            'Cached ~%s (%s entries, %s branches) for %s',
            flatanalysis.utils.bytes_to_human_readable(self.sizeof_cache), self.numentries_cache, len(branches), self
            )

    def clear_cache(self):
        self.cache = []
        
    def get_event(self, i=0, **kwargs):
        i_entry_start = 0
        for arrays in self.iterate(**kwargs):
            i_entry_stop = i_entry_start + numentries(arrays) - 1
            if i > i_entry_stop:
                i_entry_start = i_entry_stop + 1
                continue
            # Cut out the one entry we're interested in in a new arrays
            return { k : v[i-i_entry_start:i-i_entry_start+1] for k, v in arrays.items() }
        else:
            raise Exception(
                'Requested entry {0} not in range; reached end of stored events at entry {1}'
                .format(i, i_entry_stop)
                )

    def numentries(self, use_cache=True):
        if use_cache:
            return self.numentries_cache
        else:
            if not hasattr(self, 'numentries'):
                logger.info(
                    'Calculating numentries for %s rootfiles in %s',
                    len(self.rootfiles), self.shortname()
                    )
                self._numentries = 0
                for rootfile in self.rootfiles:
                    self._numentries += uproot.open(rootfile).get(self.treename).numentries
            return self._numentries



# ______________________________________________________________________
# Basic analysis things

DEFAULT_BRANCHES = [
    b'JetsAK15',
    b'JetsAK15_softDropMass',
    # b'JetsAK15_subjets',
    b'TriggerPass',
    b'MET',
    b'METPhi',
    b'HT',
    ]

def basic_svj_analysis_branches(arrays):
    """
    Standard analysis array operations that we want to run basically always
    """
    flatanalysis.arrayutils.filter_zerojet_events(arrays)
    flatanalysis.arrayutils.get_leading_and_subleading_jet(arrays)
    flatanalysis.arrayutils.get_jet_closest_to_met(arrays)
    # flatanalysis.arrayutils.get_summedsoftdropsubjets(arrays)
    flatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15')
    flatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_leading')
    flatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_subleading')
    flatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_closest')
    flatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_subclosest')

class SVJDataset(Dataset):
    """
    SVJ-specific things about the dataset
    """
    def __init__(self, name, rootfiles, *args, **kwargs):
        self.name = name
        # Make sure the default branches are in there
        if 'branches' in kwargs:
            kwargs['branches'] = list(set(kwargs['branches'] + DEFAULT_BRANCHES))
        else:
            kwargs['branches'] = DEFAULT_BRANCHES[:]
        # Set a default value for the number of entries to be kept in the cache
        if not 'max_entries' in kwargs: kwargs['max_entries'] = 1000
        super().__init__(rootfiles, *args, **kwargs)
        # Apply the standard calculations on the cache
        if self.cache:
            for arrays in self.cache:
                basic_svj_analysis_branches(arrays) 

    def __repr__(self):
        return super().__repr__().replace('Dataset', 'Dataset ' + self.get_category())

    def shortname(self):
        return self.name[:20]

    def get_weight(self, use_cache=True):
        """
        Requires get_xs(self) implementation
        """
        numentries = self.numentries(use_cache)
        if numentries == 0: return 0.0
        return self.get_xs() / float(numentries)

    def get_category(self):
        return self.name


class SignalDataset(SVJDataset):
    def get_xs(self):
        if hasattr(self, 'xs'): return self.xs
        self.xs = 100.
        return self.xs


class BackgroundDataset(SVJDataset):

    titles = {
        'ttjets' : r'$t\bar{t}$',
        'qcd' : 'QCD',
        'zjets' : 'Z+jets',
        'wjets' : 'W+jets',
        }

    def get_xs(self):
        if hasattr(self, 'xs'): return self.xs
        datasetname = self.name.split('.', 1)[1]
        if '_ext' in datasetname: datasetname = datasetname.split('_ext')[0]
        self.xs = flatanalysis.crosssections.get_xs(datasetname)
        if self.xs is None:
            raise RuntimeError('No cross section for {0}'.format(self.name))
        return self.xs

    def get_category(self):
        for key in [ 'ttjets', 'qcd', 'zjets', 'wjets' ]:
            if key in self.name.lower():
                return key
        else:
            return super().get_category()
            # raise Exception(
            #     'No category could be determined from name {0}'.format(self.name)
            #     )

    def get_title(self):
        return self.titles.get(self.get_category(), self.get_category())