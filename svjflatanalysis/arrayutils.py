import numpy as np, awkward
from math import pi

import svjflatanalysis
logger = svjflatanalysis.logger


def safe_divide(a, b):
    """
    Returns a/b, but puts 0 where b==0
    """
    return np.divide(a, b, out=np.zeros_like(a), where=b>0.)


feature_labels = [
    'met',
    'subleading_pt',
    'subleading_eta',
    'subleading_energy',
    'subleading_mt',
    'subleading_rt',
    'subleading_rtnew',
    'subleading_msd',
    'subleading_axismajor',
    'subleading_axisminor',
    'subleading_ecfN2b1',
    'subleading_ecfN2b2',
    'subleading_girth',
    'subleading_ptd',
    'leading_pt',
    'leading_eta',
    'leading_energy',
    'leading_mt',
    'leading_rt',
    'leading_rtnew',
    'leading_msd',
    'leading_axismajor',
    'leading_axisminor',
    'leading_ecfN2b1',
    'leading_ecfN2b2',
    'leading_girth',
    'leading_ptd',
    'subleading_tau21',
    'leading_tau21',
    'deltaeta',
    'dphimet',
    ]

def to_feature_array(arrays):
    """
    Transforms a dict-like arrays object to a flat numpy array for ML purposes
    """
    # Simply direct variables from the arrays object
    plain_features = [
        b'MET',
        b'JetsAK15_subleading.fCoordinates.fPt',
        b'JetsAK15_subleading.fCoordinates.fEta',
        b'JetsAK15_subleading.fCoordinates.fE',
        b'JetsAK15_subleading_MT',
        b'JetsAK15_subleading_RT',
        b'JetsAK15_subleading_RTnew',
        b'JetsAK15_subleading_softDropMass',
        b'JetsAK15_subleading_axismajor',
        b'JetsAK15_subleading_axisminor',
        b'JetsAK15_subleading_ecfN2b1',
        b'JetsAK15_subleading_ecfN2b2',
        b'JetsAK15_subleading_girth',
        b'JetsAK15_subleading_ptD',
        b'JetsAK15_leading.fCoordinates.fPt',
        b'JetsAK15_leading.fCoordinates.fEta',
        b'JetsAK15_leading.fCoordinates.fE',
        b'JetsAK15_leading_MT',
        b'JetsAK15_leading_RT',
        b'JetsAK15_leading_RTnew',
        b'JetsAK15_leading_softDropMass',
        b'JetsAK15_leading_axismajor',
        b'JetsAK15_leading_axisminor',
        b'JetsAK15_leading_ecfN2b1',
        b'JetsAK15_leading_ecfN2b2',
        b'JetsAK15_leading_girth',
        b'JetsAK15_leading_ptD',
        ]
    # Quantities that have to be calculated from branches
    calculated_features = [
        # tau21
        lambda arrays: arrays[b'JetsAK15_subleading_NsubjettinessTau2'] / arrays[b'JetsAK15_subleading_NsubjettinessTau1'],
        lambda arrays: arrays[b'JetsAK15_leading_NsubjettinessTau2'] / arrays[b'JetsAK15_leading_NsubjettinessTau1'],
        svjflatanalysis.arrayutils.calculate_deltaeta,
        svjflatanalysis.arrayutils.calculate_dphimet,
        ]
    # Build the returnable numpy array
    plain = [ arrays[b].flatten() for b in plain_features ]
    for fn in calculated_features:
        plain.append(fn(arrays).flatten())
    plain = np.array(plain)
    return plain.T

# ------------------------
# Event level

def numentries(arrays):
    """
    Counts the number of entries in a typical arrays from a ROOT file,
    by looking at the length of the first key
    """
    return arrays[list(arrays.keys())[0]].shape[0]

def select(arrays, selection, inplace=True):
    """
    If inplace is True, modifies arrays to only contain the selection
    Otherwise it creates a new arrays object with the selection in place
    """
    if inplace:
        for key in arrays:
            arrays[key] = arrays[key][selection]
    else:
        return { k : v[selection] for k, v in arrays.items() }

_WARNED_ABOUT_TRIGGERS=[]
def get_trigger_indices(triggers):
    """
    Maps trigger names as strings to trigger indices as defined in the 'TriggerPass' key
    """
    trigger_titles = svjflatanalysis.get_trigger_titles()
    # Allow a year to be passed
    if triggers in [ 2016, 2017, 2018 ]:
        triggers = svjflatanalysis.trigger.get_triggers_for_year(triggers)
    # Indices of the triggers we want to apply
    trigger_indices = []
    for trigger in triggers:
        if not trigger in trigger_titles:
            # Only warn once per trigger to avoid tons of warning messages
            global _WARNED_ABOUT_TRIGGERS
            if not(trigger in _WARNED_ABOUT_TRIGGERS):
                logger.warning('Trigger {} not a known trigger - ignoring!'.format(trigger))
                _WARNED_ABOUT_TRIGGERS.append(trigger)
            continue
        trigger_indices.append(trigger_titles.index(trigger))
    return np.array(trigger_indices)

def passed_trigger(arrays, triggers, trigger_pass_branch=b'TriggerPass'):
    """
    Returns a boolean array, True if the event passed the given triggers or False if not
    """
    trigger_indices = get_trigger_indices(triggers)
    # Check, per event, if any of the triggers we want equals 1
    trigger_decisions = arrays[trigger_pass_branch]
    passes = (trigger_decisions[:,trigger_indices] == 1).any(axis=1)
    return passes

def apply_trigger(arrays, triggers, return_triggered=True, return_counts=True):
    """
    Given a structure `arrays` as given by uproot, applies triggers and returns a filtered arrays.
    The branch b'TriggerPass' must be in arrays
    """
    passes = passed_trigger(arrays, triggers)
    n_total = passes.shape[0]
    n_pass = passes.nonzero()[0].shape[0]
    if return_triggered:
        arrays_triggered = select(arrays, passes, inplace=False)
        if return_counts:
            return arrays_triggered, n_pass, n_total
        else:
            return arrays_triggered
    else:
        return n_pass, n_total

def passes_trigger_and_jetpt550(arrays, triggers=2018):
    jets = get_jets(arrays, b'JetsAK15')
    passes = passed_trigger(arrays, triggers) & (jets.pt > 550.).any()
    return passes

def apply_trigger_and_jetpt550(arrays, triggers=2018, return_counts=False):
    passes = passes_trigger_and_jetpt550(arrays, triggers)
    if return_counts:
        return passes.sum(), passes.counts()
    else:
        return select(arrays, passes, inplace=False)

def count_triggers(trigger_decisions, triggers):
    """
    Counts the number of events that pass a trigger
    """
    trigger_indices = get_trigger_indices(triggers)
    passes = (trigger_decisions[:,trigger_indices] == 1).any(axis=1)
    # n_total = trigger_decisions.shape[0]
    n_pass = passes.nonzero()[0].shape[0]
    return n_pass #, n_total

def filter_zerojet_events(arrays, inplace=True):
    jets = get_jets(arrays, b'JetsAK15')
    passes = (jets.counts >= 1)
    return select(arrays, passes, inplace)

def is_nested(arrays):
    """
    Checks whether an array is of the old, nested style type, or the newer non-nested
    """
    if numentries(arrays) == 0: return False
    return not(arrays[b'JetsAK15'].dtype == 'int32')

def get_jets(arrays, name=b'JetsAK15', attributes=None, has_subjets=True):
    """
    Gets a jets collection.
    Compatible with nested and non-nested arrays
    """
    if is_nested(arrays):
        jets = arrays[name]
    else:
        if attributes is None:
            attributes = [
                'axismajor',
                'axisminor',
                'ecfN2b1',
                'ecfN2b2',
                'ecfN3b1',
                'ecfN3b2',
                'girth',
                'NsubjettinessTau1',
                'NsubjettinessTau2',
                'NsubjettinessTau3',
                'ptD',
                ]
        jets = Jets(name, arrays, attributes, has_subjets=has_subjets)
    return jets

def apply_window(arrays, left, right, branch=b'JetsAK15_subleading_MT'):
    selection = (
        (arrays[b'JetsAK15_subleading_MT'] > left) & (arrays[b'JetsAK15_subleading_MT'] < right)
        ).any()
    return select(arrays, selection, inplace=False)

# ------------------------
# In-place array modifications

class Jets(object):
    """
    Fakes a TLorentzVector-style jet object in ntuples with a high splitlevel
    """
    subtructure_branches = [
        'axismajor',
        'axisminor',
        'ecfN2b1',
        'ecfN2b2',
        'ecfN3b1',
        'ecfN3b2',
        'girth',
        'NsubjettinessTau1',
        'NsubjettinessTau2',
        'NsubjettinessTau3',
        'ptD',
        ]

    @classmethod
    def from_ptetaphie(cls, name, pt, eta, phi, energy):
        instance = cls(name, attributes=[], has_subjets=False, has_substructure=False)
        instance.pt = pt
        instance.eta = eta
        instance.phi = phi
        instance.energy = energy
        instance.mass2 = instance._mag2()
        instance.mass = instance._mass()
        return instance

    def __init__(self, name, arrays=None, attributes=None, has_subjets=True, has_preprocessed_subjets=False, has_substructure=True):
        self.has_subjets = has_subjets
        self.has_substructure = has_substructure
        self.name = name
        self.name_decoded = name.decode('utf-8')
        self.attributes = [] if attributes is None else attributes

        if arrays:
            self.counts = arrays[self.name]
            self.pt = arrays[add_to_bytestring(self.name, '.fCoordinates.fPt')]
            self.eta = arrays[add_to_bytestring(self.name, '.fCoordinates.fEta')]
            self.phi = arrays[add_to_bytestring(self.name, '.fCoordinates.fPhi')]
            self.energy = arrays[add_to_bytestring(self.name, '.fCoordinates.fE')]
            self.mass2 = self._mag2()
            self.mass = self._mass()
            self.softdropmass = arrays[add_to_bytestring(self.name, '_softDropMass')]
            
            if self.has_subjets:
                if has_preprocessed_subjets:                
                    self.n_subjets = arrays[add_to_bytestring(self.name, '_subjetsCounts')]
                    self.subjet_pt = arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fPt')]
                    self.subjet_eta = arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fEta')]
                    self.subjet_phi = arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fPhi')]
                    self.subjet_energy = arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fE')]
                else:
                    for offsets_branch in [ '_subjetsOffsets', '_subjetsCounts' ]:
                        if add_to_bytestring(self.name, offsets_branch) in arrays.keys():
                            break
                    else:
                        raise Exception('Could not determine any subjets offset branch in arrays')
                    self.offsets_branch = offsets_branch
                    self.n_subjets = arrays[add_to_bytestring(self.name, offsets_branch)]
                    self.subjet_pt = self._doublejagged_from_nsubjet(arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fPt')])
                    self.subjet_eta = self._doublejagged_from_nsubjet(arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fEta')])
                    self.subjet_phi = self._doublejagged_from_nsubjet(arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fPhi')])
                    self.subjet_energy = self._doublejagged_from_nsubjet(arrays[add_to_bytestring(self.name, '_subjets.fCoordinates.fE')])

            if self.has_substructure:
                for attr in self.subtructure_branches:
                    setattr(self, attr, arrays[(self.name_decoded + '_' + attr).encode('utf-8')])

    def _mag2(self):
        # Copied from https://github.com/scikit-hep/uproot-methods/blob/master/uproot_methods/classes/TLorentzVector.py
        x = self.pt * np.cos(self.phi)
        y = self.pt * np.sin(self.phi)
        z = self.pt * np.sinh(self.eta)
        return self.energy*self.energy - (x*x + y*y + z*z)

    def _mass(self):
        # Protect sqrt from <0 errors
        mag2 = self._mag2()
        mag2[mag2 > 0.] = np.sqrt(mag2[mag2 > 0.])
        mag2[mag2 < 0.] = 0.
        return mag2
        # return np.sqrt(mag2, out=np.zeros_like(mag2), where=mag2>0.)
    
    def _doublejagged_from_nsubjet(self, array):
        """
        Returns a 'nested' JaggedArray based on n_subjets
        """
        return awkward.JaggedArray.fromoffsets(
            self.n_subjets.offsets, awkward.JaggedArray.fromcounts(self.n_subjets.content, array.content)
            )

    def _all_attrs(self):
        attributes = ['pt', 'eta', 'phi', 'energy', 'mass2', 'mass', 'softdropmass'] + self.attributes
        if self.has_subjets:
            attributes += ['n_subjets', 'subjet_pt', 'subjet_eta', 'subjet_phi', 'subjet_energy']
        if self.has_substructure:
            attributes += self.subtructure_branches
        return attributes

    def __getitem__(self, *args, **kwargs):
        """
        Calls __getitem__ on all the attributes
        """
        selection = Jets(self.name)
        for attr in self._all_attrs():
            try:
                setattr(selection, attr, getattr(self, attr).__getitem__(*args, **kwargs))            
            except IndexError:
                logger.error(
                    'IndexError on attribute \'{0}\';'
                    '\n    .{0} = {1}'
                    '\n    args = {2}, kwargs = {3}'
                    .format(attr, getattr(self, attr), args, kwargs)
                    )
                raise
        selection.counts = selection.pt.counts
        return selection

    def save(self, arrays, tag):
        """
        Puts this collection of jets into the arrays
        """
        name = add_to_bytestring(self.name, tag)
        arrays[name] = self.counts
        arrays[add_to_bytestring(name, '.fCoordinates.fPt')] = self.pt
        arrays[add_to_bytestring(name, '.fCoordinates.fEta')] = self.eta
        arrays[add_to_bytestring(name, '.fCoordinates.fPhi')] = self.phi
        arrays[add_to_bytestring(name, '.fCoordinates.fE')] = self.energy
        arrays[add_to_bytestring(name, '_softDropMass')] = self.softdropmass
        if self.has_subjets:
            arrays[add_to_bytestring(name, '_subjets.fCoordinates.fPt')] = self.subjet_pt
            arrays[add_to_bytestring(name, '_subjets.fCoordinates.fEta')] = self.subjet_eta
            arrays[add_to_bytestring(name, '_subjets.fCoordinates.fPhi')] = self.subjet_phi
            arrays[add_to_bytestring(name, '_subjets.fCoordinates.fE')] = self.subjet_energy
            arrays[add_to_bytestring(name, '_subjetsCounts')] = self.n_subjets
        if self.has_substructure:
            for attr in self.subtructure_branches:
                key = (self.name_decoded + tag + '_' + attr).encode('utf-8')
                arrays[key] = getattr(self, attr)
        for attr in self.attributes:
            key = (self.name_decoded + tag + '_' + attr).encode('utf-8')
            arrays[key] = getattr(self, attr)


def add_to_bytestring(bytestring, tag):
    normal_string = bytestring.decode('utf-8')
    normal_string += tag
    return normal_string.encode('utf-8')

def nonnested_branches(jets_branch=b'JetsAK15', add_subjets=False, old_style=False):
    tags = ['.fCoordinates.fPt', '.fCoordinates.fEta', '.fCoordinates.fPhi', '.fCoordinates.fE']
    branches = [add_to_bytestring(jets_branch, tag) for tag in tags]
    if add_subjets:
        branches.append(add_to_bytestring(jets_branch, '_subjets'))
        branches.append(add_to_bytestring(jets_branch, '_subjetsOffsets' if old_style else '_subjetsCounts'))
        branches.extend(nonnested_branches(add_to_bytestring(jets_branch, '_subjets'), add_subjets=False))
    return branches

def get_leading_jet(arrays, jets_branch=b'JetsAK15'):
    jets = get_jets(arrays, jets_branch)
    arrays[add_to_bytestring(jets_branch, '_leading')] = arrays[jets_branch][jets.pt.argmax()]

def get_leading_and_subleading_jet(arrays, jets_branch=b'JetsAK15'):
    jets = Jets(jets_branch, arrays)
    leading_pt_index = jets.pt.argmax()
    # Set the max pt to a very low value in this copy,
    # so that the subleading pt becomes the new maximum
    pt_copy = jets.pt[:,:]
    pt_copy[leading_pt_index] = -1e5
    subleading_pt_index = pt_copy.argmax()
    # Clearly the subleading shouldn't be the same one as the leading
    subleading_pt_index = subleading_pt_index[subleading_pt_index != leading_pt_index]
    if is_nested(arrays):
        arrays[add_to_bytestring(jets_branch, '_leading')] = jets[leading_pt_index]
        arrays[add_to_bytestring(jets_branch, '_subleading')] = jets[subleading_pt_index]
    else:
        jets[leading_pt_index].save(arrays, '_leading')
        jets[subleading_pt_index].save(arrays, '_subleading')

def get_jet_closest_to_met(arrays, jets_branch=b'JetsAK15'):
    metphi = arrays[b'METPhi']
    jets = get_jets(arrays, jets_branch)
    # Get absolute dphi between jet and met, find closest angle too
    dphi = np.abs(jets.phi - metphi)
    dphi[dphi > 2.*pi] = dphi[dphi > 2.*pi] - 2.*pi  # Whole circles subtracted
    dphi[dphi > pi] = 2.*pi - dphi[dphi > pi]  # Pick the smaller angle
    closest_index = dphi.argmin()
    # Also get subclosest
    dphi_copy = dphi[:,:]
    dphi_copy[closest_index] = 1e6
    subclosest_index = dphi_copy.argmin()
    # Clearly the subleading shouldn't be the same one as the leading
    subclosest_index = subclosest_index[subclosest_index != closest_index]
    if is_nested(arrays):
        arrays[add_to_bytestring(jets_branch, '_closest')] = jets[closest_index]
        arrays[add_to_bytestring(jets_branch, '_subclosest')] = jets[subclosest_index]
    else:
        jets[closest_index].save(arrays, '_closest')
        jets[subclosest_index].save(arrays, '_subclosest')


def get_summedsoftdropsubjets(arrays, subjets_branch=b'JetsAK15_subjets'):
    raise NotImplementedError(
        'uproot cannot handle ROOT.vector<vector<TLorentzVector>> at the moment'
        )
    subjets = arrays[subjets_branch]
    summedsubjets = bla

def calculate_mt(arrays, jets_branch=b'JetsAK15_leading'):
    metx = np.cos(arrays[b'METPhi']) * arrays[b'MET']
    mety = np.sin(arrays[b'METPhi']) * arrays[b'MET']
    mete = np.sqrt(metx**2 + mety**2) # Should be == arrays[b'MET'] actually
    
    jets = get_jets(arrays, jets_branch)
    jetsx = np.cos(jets.phi) * jets.pt
    jetsy = np.sin(jets.phi) * jets.pt
    jetse = np.sqrt(jets.mass2 + jets.pt**2)

    mt = np.sqrt( (jetse + mete)**2 - (jetsx + metx)**2 - (jetsy + mety)**2 )
    rt = arrays[b'MET'] / mt
    rtnew = arrays[b'MET'] / np.sqrt(mt**2 + jets.pt**2)
    arrays[add_to_bytestring(jets_branch, '_MT')] = mt
    arrays[add_to_bytestring(jets_branch, '_RT')] = rt
    arrays[add_to_bytestring(jets_branch, '_RTnew')] = rtnew
    
def calculate_deltaeta(arrays, leading_branch=b'JetsAK15_leading', subleading_branch=b'JetsAK15_subleading'):
    leading_jets = Jets(leading_branch, arrays, has_subjets=False, has_substructure=False)
    subleading_jets = Jets(subleading_branch, arrays, has_subjets=False, has_substructure=False)
    # Only get the leading jets for when there is a subleading jet too
    leading_jets = leading_jets[subleading_jets.counts > 0]
    subleading_jets = subleading_jets[subleading_jets.counts > 0]
    deta = np.abs(leading_jets.eta - subleading_jets.eta)
    return deta

def calculate_dphimet(arrays, jets_branch=b'JetsAK15_subleading'):
    jets = Jets(jets_branch, arrays, has_subjets=False, has_substructure=False)
    dphi = np.abs(jets.phi - arrays[b'METPhi'])
    dphi[dphi > 2.*pi] = dphi[dphi > 2.*pi] - 2.*pi  # Whole circles subtracted
    dphi[dphi > pi] = 2.*pi - dphi[dphi > pi]  # Pick the smaller angle
    return dphi
