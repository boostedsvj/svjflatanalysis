# import uproot4 as uproot
import uproot, re, numpy as np
import svjflatanalysis
logger = svjflatanalysis.logger

# ______________________________________________________________________
# arrays-level helpers

def add_to_bytestring(bytestring, tag):
    normal_string = bytestring.decode('utf-8')
    normal_string += tag
    return normal_string.encode('utf-8')

def numentries(arrays):
    return arrays[list(arrays.keys())[0]].shape[0]
    
def iterate_events(arrays, flat=False):
    """
    Iterates event by event from an arrays of events
    """
    n = numentries(arrays)
    for i in range(n):
        if flat:
            yield { k : v[i] for k, v in arrays.items() }
        else:
            yield { k : v[i:i+1] for k, v in arrays.items() }

def group_per_category(bkgs):
    """
    Groups a flat list of datasets into sublists with the same category
    E.g. [ ttjet, ttjet, ttjet, qcd, qcd ] --> [ [ttjet, ttjet, ttjet], [qcd, qcd] ]
    """
    cats = list(set(bkg.get_category() for bkg in bkgs))
    cats.sort()
    return [ [ b for b in bkgs if b.get_category() == cat ] for cat in cats ]


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
        self.rootfiles = [rootfiles] if svjflatanalysis.utils.is_string(rootfiles) else rootfiles
        if len(self.rootfiles) == 0:
            logger.warning('No rootfiles for %s; iterators will be empty', self)
            make_cache = False
        self.cache = []
        if 'treename' in kwargs:
            self.treename = kwargs.pop('treename')
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
            if not svjflatanalysis.HAS_TQDM:
                logger.error('tqdm could not be imported, progressbars are disabled')
            else:
                iterator = svjflatanalysis.tqdm(iterator, total=total, desc='arrays' if use_cache else 'root files')
        for arrays in iterator:
            yield arrays

    def iterate_batch(self, batch_size=32, **kwargs):
        """
        Like self.iterate, but yields fixed sized chunks of events
        (except the last chunk, which will be however many events are left)
        """
        import awkward
        # Helper function to split an arrays
        def split(arrays, n):
            first  = { k : v[:n] for k, v in arrays.items() }
            second = { k : v[n:] for k, v in arrays.items() }
            return first, second
        # Helper function to merge an arrays
        def merge(batch):
            r = {}
            for key in batch[0].keys():
                r[key] = awkward.concatenate((b[key] for b in batch))
            return r
        # Initialize some variables
        arrays_iterator = self.iterate(**dict(kwargs, entrysteps=batch_size))
        n_todo = batch_size
        batch = []
        arrays = next(arrays_iterator)
        n_block = svjflatanalysis.arrayutils.numentries(arrays)
        # Build batches
        while True:
            if n_block > n_todo:
                # Take off the needed number of chunks to fill the batch
                for_batch, arrays = split(arrays, n_todo)
                batch.append(for_batch)
                n_block -= n_todo
                n_todo = 0
            elif n_block <= n_todo:
                # No need to split, just add the whole arrays to the batch and open a new arrays
                batch.append(arrays)
                n_todo -= n_block
                try:
                    arrays = next(arrays_iterator)
                except StopIteration:
                    # No more new arrays to take
                    break
                # Reset block counter
                n_block = svjflatanalysis.arrayutils.numentries(arrays)
            if n_todo == 0:
                # Batch is ready
                if len(batch) == 0:
                    raise Exception('No events in batch')
                elif len(batch) == 1:
                    # No need to merge
                    yield batch[0]
                elif len(batch) > 1:
                    yield merge(batch)
                # Reset todo counter and the batch
                n_todo = batch_size
                batch = []
        # Yield the remaining events in the batch if there are any
        if len(batch) == 1:
            # No need to merge
            yield batch[0]
        elif len(batch) > 1:
            yield merge(batch)

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
        if max_entries: kwargs['entrysteps'] = max_entries
        for arrays in self.iterate(use_cache=False, **kwargs):
            if branches is None: branches = list(arrays.keys())
            numentries = svjflatanalysis.arrayutils.numentries(arrays)
            if max_entries and self.numentries_cache + numentries > max_entries:
                # Cut away some entries to get to max_entries
                needed_entries = max_entries - self.numentries_cache
                arrays = { k : v[:needed_entries] for k, v in arrays.items() }
            self.cache.append(arrays)
            self.sizeof_cache += sum([ v.nbytes for v in arrays.values() ])
            self.numentries_cache += svjflatanalysis.arrayutils.numentries(arrays)
            if max_entries and self.numentries_cache >= max_entries:
                break
        if branches is None:
            logger.error('Did not make cache - no arrays to loop over')
        logger.info(
            'Cached ~%s (%s entries, %s branches) for %s',
            svjflatanalysis.utils.bytes_to_human_readable(self.sizeof_cache),
            self.numentries_cache,
            '?' if branches is None else len(branches),
            self
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

def basic_svj_analysis_branches(arrays):
    """
    Standard analysis array operations that we want to run basically always
    """
    # First check if the branches are already added:
    if b'JetsAK15_subleading' in arrays: return
    svjflatanalysis.arrayutils.filter_zerojet_events(arrays)
    svjflatanalysis.arrayutils.get_leading_and_subleading_jet(arrays)
    svjflatanalysis.arrayutils.get_jet_closest_to_met(arrays)
    # svjflatanalysis.arrayutils.get_summedsoftdropsubjets(arrays)
    svjflatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15')
    svjflatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_leading')
    svjflatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_subleading')
    svjflatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_closest')
    svjflatanalysis.arrayutils.calculate_mt(arrays, b'JetsAK15_subclosest')

class SVJDataset(Dataset):
    """zzz
    SVJ-specific things about the dataset
    """
    _is_signal = None
    def __init__(self, name, rootfiles, *args, **kwargs):
        self.name = name
        # Set a default value for the number of entries to be kept in the cache
        if not 'max_entries' in kwargs: kwargs['max_entries'] = 1000
        super().__init__(rootfiles, *args, **kwargs)
        # # Apply the standard calculations on the cache
        # if self.cache:
        #     for arrays in self.cache:
        #         basic_svj_analysis_branches(arrays) 

    def __repr__(self):
        return super().__repr__().replace('Dataset', 'Dataset ' + self.get_category())

    def iterate(self, *args, **kwargs):
        """
        Like the super iterate, but applies some default branches and
        adds some more calculated branches
        """
        # Insert the default branches
        if not 'branches' in kwargs: kwargs['branches'] = svjflatanalysis.samples.svj_branches
        for arrays in super().iterate(*args, **kwargs):
            basic_svj_analysis_branches(arrays) 
            yield arrays

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

    def is_signal(self):
        return self._is_signal

    def is_bkg(self):
        return not(self._is_signal)


class SignalDataset(SVJDataset):
    _is_signal = True
    def get_xs(self):
        if hasattr(self, 'xs'): return self.xs
        self.xs = 100.
        return self.xs

    def get_title(self):
        match = re.search(r'mz(\d+)', self.name)
        if not match:
            return self.name
        else:
            return r'$m_{{Z\prime}}={}$ GeV'.format(match.group(1))

class BackgroundDataset(SVJDataset):
    _is_signal = False
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
        self.xs = svjflatanalysis.crosssections.get_xs(datasetname)
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


class FeatureDataset():
    pass





def n_events_weighted_by_xs(n_target, datasets):
    xss = np.array([d.get_xs() for d in datasets])
    n_events_per_dataset = xss / np.sum(xss) * n_target
    for dataset, n_dataset in zip(datasets, n_events_per_dataset):
        logger.info('{:4d} ({:7.2f}) events for {}'.format(int(n_dataset), n_dataset, dataset.name))
    n_events_per_dataset = n_events_per_dataset.astype(np.int64)
    return n_events_per_dataset

def build_feature_array(datasets, n=100, use_cache=True):
    """
    Takes a list of datasets and returns a flat numpy array with features per event.
    Tries to weight by cross section
    """
    if hasattr(datasets, 'treename'): datasets = [datasets]
    # Determine number of events to take per dataset, depending on cross section
    n_events_per_dataset = n_events_weighted_by_xs(n, datasets)
    # Build the array
    to_flatten = []
    y = []
    for dataset, n_dataset in zip(datasets, n_events_per_dataset):
        if n_dataset == 0: continue
        logger.info('Fetching %s events for dataset %s...', n_dataset, dataset.name)
        if use_cache:
            arrays = next(dataset.iterate_batch(n_dataset))
        else:
            arrays = next(dataset.iterate_batch(
                n_dataset,
                use_cache=False,
                branches=svjflatanalysis.samples.svj_branches
                ))
            basic_svj_analysis_branches(arrays)
        feature_array = svjflatanalysis.arrayutils.to_feature_array(arrays)
        n_actual = feature_array.shape[0]
        if n_actual != n_dataset:
            logger.warning(
                'Retrieved %s values rather than target %s values for dataset %s',
                n_actual, n_dataset, dataset.name
                )
        to_flatten.append(feature_array)
        y.append(np.zeros(n_actual, dtype=np.int8) if dataset.is_bkg() else np.ones(n_actual, dtype=np.int8))
        del arrays

    # Concat per dataset and return
    X = np.concatenate(to_flatten).T
    y = np.concatenate(y)
    return X, y

