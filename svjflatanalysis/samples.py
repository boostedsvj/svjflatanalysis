"""
List of datasets currently being analyzed
"""

import os, os.path as osp, glob
import seutils

import svjflatanalysis
logger = svjflatanalysis.logger

# ______________________________________________________________________
basic_branches = [
    b'JetsAK15',
    b'TriggerPass',
    b'MET',
    b'METPhi',
    b'HT',
    ]

substructure_branches = [
    b'JetsAK15_softDropMass',
    b'JetsAK15_axismajor',
    b'JetsAK15_axisminor',
    b'JetsAK15_ecfN2b1',
    b'JetsAK15_ecfN2b2',
    b'JetsAK15_ecfN3b1',
    b'JetsAK15_ecfN3b2',
    b'JetsAK15_girth',
    b'JetsAK15_NsubjettinessTau1',
    b'JetsAK15_NsubjettinessTau2',
    b'JetsAK15_NsubjettinessTau3',
    b'JetsAK15_ptD',
    ]

svj_branches = (
    basic_branches
    + svjflatanalysis.arrayutils.nonnested_branches(b'JetsAK15', add_subjets=True)
    + substructure_branches
    + ['NMuons', 'NElectrons']
    )

def set_default_branches(kwargs):
    if not 'branches' in kwargs:
        kwargs['branches'] = svj_branches


# def to_feature_array(arrays):
#     """
#     Transforms a dict-like arrays object to a flat numpy array for ML purposes
#     """


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

# ______________________________________________________________________
# Trigger efficiencies
# both for just the trigger, and trigger + jetpt>550

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

NOCUTS_TRIGGER_EFF_BKG = {
    'Autumn18.TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8' : 0.04059877523247902,
    'Autumn18.TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.043112978111872646,
    'Autumn18.TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8'             : 0.035600578871201154,
    'Autumn18.TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'        : 0.7729,
    'Autumn18.TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'       : 0.9744000000000002,
    'Autumn18.TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'      : 0.9999000000000001,
    'Autumn18.TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'       : 1.0,
    'Autumn18.QCD_Pt_80to120_TuneCP5_13TeV_pythia8'                        : 0.0,
    'Autumn18.QCD_Pt_120to170_TuneCP5_13TeV_pythia8'                       : 0.0,
    'Autumn18.QCD_Pt_170to300_TuneCP5_13TeV_pythia8'                       : 0.00247269730063878,
    'Autumn18.QCD_Pt_300to470_TuneCP5_13TeV_pythia8'                       : 0.1838183818381838,
    'Autumn18.QCD_Pt_470to600_TuneCP5_13TeV_pythia8'                       : 0.9240924092409241,
    'Autumn18.QCD_Pt_600to800_TuneCP5_13TeV_pythia8'                       : 0.9981,
    'Autumn18.QCD_Pt_800to1000_TuneCP5_13TeV_pythia8_ext1'                 : 0.9996999999999999,
    'Autumn18.QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8'                     : 0.9999,
    'Autumn18.QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8'                      : 1.0,
    'Autumn18.WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8'     : 0.0,
    'Autumn18.WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.0,
    'Autumn18.WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.00021344717182497332,
    'Autumn18.WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.022404092071611252,
    'Autumn18.WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.16463292658531706,
    'Autumn18.WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8'   : 0.6910000000000001,
    'Autumn18.WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8'  : 0.999,
    'Autumn18.WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8'   : 1.0,
    'Autumn18.ZJetsToNuNu_HT-100To200_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-200To400_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-400To600_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-600To800_13TeV-madgraph'                      : 0.14755902360944378,
    'Autumn18.ZJetsToNuNu_HT-800To1200_13TeV-madgraph'                     : 0.6389,
    'Autumn18.ZJetsToNuNu_HT-1200To2500_13TeV-madgraph'                    : 0.9964999999999999,
    'Autumn18.ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph'                     : 1.0,
    }

NOCUTS_TRIGGER_PLUS_JETPT550_EFF_BKG = {
    'Autumn18.TTJets_SingleLeptFromTbar_TuneCP5_13TeV-madgraphMLM-pythia8' : 0.015649807212519844,
    'Autumn18.TTJets_SingleLeptFromT_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.015476453681185054,
    'Autumn18.TTJets_DiLept_TuneCP5_13TeV-madgraphMLM-pythia8'             : 0.01678726483357453,
    'Autumn18.TTJets_HT-600to800_TuneCP5_13TeV-madgraphMLM-pythia8'        : 0.4633,
    'Autumn18.TTJets_HT-800to1200_TuneCP5_13TeV-madgraphMLM-pythia8'       : 0.8111,
    'Autumn18.TTJets_HT-1200to2500_TuneCP5_13TeV-madgraphMLM-pythia8'      : 0.9925,
    'Autumn18.TTJets_HT-2500toInf_TuneCP5_13TeV-madgraphMLM-pythia8'       : 1.0,
    'Autumn18.QCD_Pt_80to120_TuneCP5_13TeV_pythia8'                        : 0.0,
    'Autumn18.QCD_Pt_120to170_TuneCP5_13TeV_pythia8'                       : 0.0,
    'Autumn18.QCD_Pt_170to300_TuneCP5_13TeV_pythia8'                       : 0.00010302905419328252,
    'Autumn18.QCD_Pt_300to470_TuneCP5_13TeV_pythia8'                       : 0.0288028802880288,
    'Autumn18.QCD_Pt_470to600_TuneCP5_13TeV_pythia8'                       : 0.5894589458945895,
    'Autumn18.QCD_Pt_600to800_TuneCP5_13TeV_pythia8'                       : 0.9843000000000001,
    'Autumn18.QCD_Pt_800to1000_TuneCP5_13TeV_pythia8_ext1'                 : 0.9982,
    'Autumn18.QCD_Pt_1000to1400_TuneCP5_13TeV_pythia8'                     : 0.9997,
    'Autumn18.QCD_Pt_1400to1800_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_1800to2400_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_2400to3200_TuneCP5_13TeV_pythia8'                     : 1.0,
    'Autumn18.QCD_Pt_3200toInf_TuneCP5_13TeV_pythia8'                      : 1.0,
    'Autumn18.WJetsToLNu_HT-70To100_TuneCP5_13TeV-madgraphMLM-pythia8'     : 0.0,
    'Autumn18.WJetsToLNu_HT-100To200_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.0,
    'Autumn18.WJetsToLNu_HT-200To400_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.0,
    'Autumn18.WJetsToLNu_HT-400To600_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.0034782608695652175,
    'Autumn18.WJetsToLNu_HT-600To800_TuneCP5_13TeV-madgraphMLM-pythia8'    : 0.05071014202840567,
    'Autumn18.WJetsToLNu_HT-800To1200_TuneCP5_13TeV-madgraphMLM-pythia8'   : 0.3254,
    'Autumn18.WJetsToLNu_HT-1200To2500_TuneCP5_13TeV-madgraphMLM-pythia8'  : 0.9581,
    'Autumn18.WJetsToLNu_HT-2500ToInf_TuneCP5_13TeV-madgraphMLM-pythia8'   : 1.0,
    'Autumn18.ZJetsToNuNu_HT-100To200_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-200To400_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-400To600_13TeV-madgraph'                      : 0.0,
    'Autumn18.ZJetsToNuNu_HT-600To800_13TeV-madgraph'                      : 0.057823129251700675,
    'Autumn18.ZJetsToNuNu_HT-800To1200_13TeV-madgraph'                     : 0.3083,
    'Autumn18.ZJetsToNuNu_HT-1200To2500_13TeV-madgraph'                    : 0.9479999999999998,
    'Autumn18.ZJetsToNuNu_HT-2500ToInf_13TeV-madgraph'                     : 0.9998999999999999,
    }


def init_sig_ht1000(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/klijnsma/semivis/treemakerht1000/mz{0}_*.root'
        .format(int(mz))
        )
    if not 'max_entries' in kwargs: kwargs['max_entries'] = None
    set_default_branches(kwargs)
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
    set_default_branches(kwargs)
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz]
    return signal

def init_sig_2017_nohtcut(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/nohtcut_year2017/*mz{}*/*.root'
        .format(int(mz))
        )
    set_default_branches(kwargs)
    signal = svjflatanalysis.dataset.SignalDataset(name, rootfiles, **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz]
    return signal

def init_sig_2016_nohtcut(mz, **kwargs):
    name = 'mz' + str(mz)
    rootfiles = seutils.ls_wildcard(
        'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/nohtcut_year2016/*mz{}*/*.root'
        .format(int(mz))
        )
    set_default_branches(kwargs)
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
    
triggered_path = (
    'root://cmseos.fnal.gov//store/user/lpcsusyhad/SVJ2017/boosted/treemaker/triggered_and_jetpt550_Nov10'
    if os.uname()[1].endswith('.fnal.gov') else 
    '/Users/klijnsma/work/svj/flat/data/triggered_and_jetpt550_Nov10/'
    )

def init_sig_triggered(year, mz, **kwargs):
    rootfiles = [osp.join(triggered_path, 'year{}_mz{}.root'.format(year, mz))]
    signal = svjflatanalysis.dataset.SignalDataset('mz{}_year{}'.format(mz, year), rootfiles, treename='PreSelection', **kwargs)
    signal.xs = SIGNAL_CROSSSECTIONS[mz] * NOCUTS_TRIGGER_PLUS_JETPT550_EFF[mz]
    return signal
    
def init_sigs_triggered(year, **kwargs):
    return [init_sig_triggered(year, mz, **kwargs) for mz in [150, 250]]

def init_bkgs_triggered(category=None, **kwargs):
    bkgs = []
    path = osp.join(triggered_path, 'Autumn18.*')
    bkg_rootfiles = seutils.ls_wildcard(path) if os.uname()[1].endswith('.fnal.gov') else glob.glob(path)
    kwargs['treename'] = 'PreSelection'
    for bkg_rootfile in bkg_rootfiles:
        if category:
            if not category.lower() in bkg_rootfile.lower():
                logger.info('Skipping %s', bkg_rootfile)
                continue
        logger.info('Loading %s', bkg_rootfile)
        name = osp.basename(bkg_rootfile).replace('.root', '')
        bkg = svjflatanalysis.dataset.BackgroundDataset(name, [bkg_rootfile], **kwargs)
        bkg.get_xs()
        # Take into account the trigger efficiencies
        if not name in NOCUTS_TRIGGER_PLUS_JETPT550_EFF_BKG:
            logger.error('FIXME: efficiency for %s not calculated yet!', name)
            continue
        bkg.xs *= NOCUTS_TRIGGER_PLUS_JETPT550_EFF_BKG[name]
        bkgs.append(bkg)
    return bkgs

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

