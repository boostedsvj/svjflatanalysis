"""
List of datasets currently being analyzed
"""

import os.path as osp
import seutils

import svjflatanalysis
logger = svjflatanalysis.logger

# ______________________________________________________________________
# Signal

# Cross sections are estimates
SIGNAL_CROSSSECTIONS = {
    150 : 2300.,
    250 : 800.,
    450 : 200.,
    650 : 80.,
    }

HT1000_EFFICIENCIES = {
    150 : 0.00030,
    250 : 0.00102,
    450 : 0.00482,
    650 : 0.01542,
    }

# Most reliable according do Sarah
SIGNAL_CROSSSECTIONS_PYTHIA_SARA = { # pb
    50  : 144500,
    100 : 19970,
    150 : 5980,
    200 : 2386,
    250 : 1135,
    300 : 638.4,
    }

SIGNAL_CROSSSECTIONS_SVJ_SARA = { # pb
    50  : 108528.67,
    100 : 15142.18,
    150 : 4443.26,
    200 : 1812.38,
    250 : 876.96,
    300 : 487.20,
    }

# See genjet250_filter_eff.py
# mz150:
# n_pass = 23464, n_total = 6881863, eff = 0.0034
# mz250:
# n_pass = 47475, n_total = 5366755, eff = 0.0088
GENJET250_EFFICIENCIES = {
    150 : 0.0034,
    250 : 0.0088,
    }


NOCUTS_TRIGGER_EFF = {
    150 : 0.0004792296242839746,
    250 : 0.0014670666229739606,
    450 : 0.005050597137387183,
    650 : 0.012032736120327363,
    }

NOCUTS_TRIGGER_PLUS_JETPT550_EFF = {
    150 : 0.0001465878850750981,
    250 : 0.0005797685967079227,
    450 : 0.002224450724769806,
    650 : 0.005264322052643222,
    }


def init_sig_ht1000(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/klijnsma/semivis/treemakerht1000/mz{0}_*.root'
        .format(int(mz))
        )
    if not 'max_entries' in kwargs: kwargs['max_entries'] = None
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz] # * HT1000_EFFICIENCIES[mz]
    return signal

def init_sigs_ht1000(**kwargs):
    return [ init_sig_ht1000(mz, **kwargs) for mz in [150, 250, 450, 650] ]


def init_sig_2018_nohtcut(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/nohtcut_year2018/*mz{}*/*.root'
        .format(int(mz))
        )
    if not 'max_entries' in kwargs: kwargs['max_entries'] = None
    kwargs['branches'] = kwargs.get('branches', []) + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True, old_style=True)
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz]
    return signal

def init_sig_2017_nohtcut(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/nohtcut_year2017/*mz{}*/*.root'
        .format(int(mz))
        )
    kwargs['branches'] = kwargs.get('branches', []) + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True)
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz]
    return signal

def init_sig_2016_nohtcut(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/nohtcut_year2016/*mz{}*/*.root'
        .format(int(mz))
        )
    kwargs['branches'] = kwargs.get('branches', []) + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True)
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz]
    return signal

def init_sigs_2018_nohtcut(**kwargs):
    return [ init_sig_2018_nohtcut(mz, **kwargs) for mz in [150, 250, 450, 650] ]

def init_sigs_2017_nohtcut(**kwargs):
    return [ init_sig_2017_nohtcut(mz, **kwargs) for mz in [150, 250, 450, 650] ]

def init_sigs_2016_nohtcut(**kwargs):
    return [ init_sig_2016_nohtcut(mz, **kwargs) for mz in [150, 250] ]

def init_sigs_nohtcut(year, **kwargs):
    return globals()['init_sigs_{}_nohtcut'.format(year)](**kwargs)


# ______________________________________________________________________
# Triggered samples

triggered_path = 'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/triggered_and_jetpt550'
def init_sig_triggered(year, mz, **kwargs):
    rootfiles = [osp.join(triggered_path, 'year{}_mz{}.root'.format(year, mz))]
    assert seutils.isfile(rootfiles[0])
    kwargs['branches'] = kwargs.get('branches', []) + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True)
    signal = svjflatanalysis.dataset.SignalDataset('mz{}_year{}'.format(mz, year), rootfiles, treename='PreSelection', **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz] * NOCUTS_TRIGGER_PLUS_JETPT550_EFF[mz]
    return signal
    
def init_sigs_triggered(year, **kwargs):
    return [init_sig_triggered(year, mz, **kwargs) for mz in [150, 250]]



# ______________________________________________________________________
# genjet250 samples

def init_sig_jetpt250(mz, year, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/jetpt250*mz{}*year{}*/*.root'
        .format(int(mz), year)
        )
    kwargs['branches'] = kwargs.get('branches', []) + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True)
    if not 'max_entries' in kwargs: kwargs['max_entries'] = None
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS_PYTHIA_SARA[mz]
    return signal

def init_sigs_2018_jetpt250(**kwargs):
    return [ init_sig_jetpt250(mz, 2018, **kwargs) for mz in [150, 250] ]

# ______________________________________________________________________
# Background

ttjets = [
    # TTJets
    # 'Autumn18.TTJets_TuneCP5_13TeV-madgraphMLM-pythia8',  # <-- All combined prob?
    'Autumn18.TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    ]

qcd = [
    # QCD Pt
    'Autumn18.QCD_Pt_80to120_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_120to170_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_170to300_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_300to470_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_470to600_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_600to800_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_800to1000_TuneCP5_13TeV_pythia8_ext1',
    'Autumn18.QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8',
    'Autumn18.QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8',
    ]

wjets = [ 
    # WJetsToLNu
    'Autumn18.WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8',
    'Autumn18.WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8',
    ]

zjets = [
    # ZJetsToNuNu
    'Autumn18.ZJetsToNuNu_HT-100To200_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-200To400_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-400To600_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-600To800_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-800To1200_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-1200To2500_13TeV-madgraph',
    'Autumn18.ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph',
    ]

all_bkgs = ttjets + qcd + wjets + zjets


BKG_ROOTFILES = None
def get_bkg_rootfiles():
    """
    Gets rootfiles for bkg once, then just returns the cached list
    """
    global BKG_ROOTFILES
    if BKG_ROOTFILES is None:
        BKG_ROOTFILES = seutils.ls_root('root://cmseos.fnal.gov//store/user/klijnsma/semivis/treemaker_bkg_May18')
    return BKG_ROOTFILES

def init_bkg(name, **kwargs):
    """
    Gets all the rootfiles that belong to this background from the cache
    """
    belongs = lambda rootfile: osp.basename(rootfile).startswith(name)
    rootfiles = [ rootfile for rootfile in get_bkg_rootfiles() if belongs(rootfile) ]
    return svjflatanalysis.dataset.BackgroundDataset(name, rootfiles, **kwargs)

def init_ttjets(**kwargs):
    return [ init_bkg(name, **kwargs) for name in ttjets ]

def init_qcd(**kwargs):
    return [ init_bkg(name, **kwargs) for name in qcd ]

def init_wjets(**kwargs):
    return [ init_bkg(name, **kwargs) for name in wjets ]

def init_zjets(**kwargs):
    return [ init_bkg(name, **kwargs) for name in zjets ]

def init_bkgs(**kwargs):
    return [ init_bkg(name, **kwargs) for name in all_bkgs ]

def init_bkgs_limited(**kwargs):
    logger.warning('Using a limited amount of background samples; shapes will be sculpted!!')
    all_bkgs = ttjets[:2] + qcd[:2] + wjets[:2] + zjets[:2]
    return [ init_bkg(name, **kwargs) for name in all_bkgs ]

def init_ttjets_test(**kwargs):
    """
    Returns a list with a single dataset, with a single root file, and a small prepared cache
    """
    rootfiles = [ rootfile for rootfile in get_bkg_rootfiles() if osp.basename(rootfile).startswith(ttjets[0]) ][:1]
    kwargs.update(make_cache=True, max_entries=50)
    return svjflatanalysis.dataset.BackgroundDataset('ttjets_testsample', rootfiles, **kwargs)

