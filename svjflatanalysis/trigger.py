import svjflatanalysis
logger = svjflatanalysis.logger

import numpy as np, os, os.path as osp
import matplotlib.pyplot as plt, mplhep

# For curve fitting
from scipy.optimize import curve_fit
from scipy.optimize import fsolve

def exp_fit_fn(x, a, b, c, d):
    return 1./(1.+np.exp(a-b*x+c*x**2-d*x**3))

def inverse_exp_fit_fn(y, a, b, c, d, x_guess=0.):
    exp_fit_fn_1dof = lambda x: exp_fit_fn(x, a, b, c, d) - y
    return fsolve(exp_fit_fn_1dof, x_guess)


def trigger_efficiency(datasets, triggers, cut_function, cut_values):
    """
    Calculates the combined trigger efficiency for a list of trigger titles over
    a list of datasets.
    `cut_function` should be a function that takes `arrays`, and returns the number
    of entries that pass according to this function and the total number of events.
    `cut_values` should be a list for which `cut_function` should be evaluated.
    """
    # If a single dataset was passed, turn it into a list
    if isinstance(datasets, svjflatanalysis.Dataset): datasets = [datasets]
    # Allow `triggers` argument to be a year
    try:
        year = int(triggers)
        if year in [ 2016, 2017, 2018 ]:
            triggers = get_triggers_for_year(year)
            logger.debug(
                'Using triggers for year %s:\n  %s',
                year, '\n  '.join(triggers)
                )
    except ValueError:
        pass
    n_pass = np.zeros(len(cut_values))
    n_total = np.zeros(len(cut_values))
    for arrays, dataset in svjflatanalysis.iterate(datasets):
        for i, cut_value in enumerate(cut_values):
            n_pass_this, n_total_this = svjflatanalysis.arrayutils.count_triggers(
                cut_function(arrays, cut_value), triggers
                )
            n_pass[i] += n_pass_this * dataset.get_weight()
            n_total[i] += n_total_this * dataset.get_weight()
    # Divide with safeguard for division by zero
    eff = np.divide(n_pass, n_total, out=np.zeros_like(n_pass), where=n_total>0.)
    return eff, n_pass, n_total

def plot_trigger_efficiency(
    cut_values, eff, cut_title='cut value', ax=None, label=None, color=None,
    include_curve_fit=False,
    include_percentile_calc=False,
    include_legend=False,
    triggers=None,
    **plot_options
    ):
    if ax is None: ax = svjflatanalysis.utils.get_ax()
    if color is None: color = svjflatanalysis.utils.get_default_color()
    line = ax.plot(cut_values, eff, **dict(color=color, marker='.', linewidth=0, **plot_options))[0]
    if label: line.set_label(label)
    ax.set_xlabel(cut_title)
    ax.set_ylabel(r'$N_{pass}$/$N_{total}$')
    # Potentially draw the curve fit and 98 percentile line
    if include_curve_fit:
        # First cut away trailing zeroes
        fit_eff = np.trim_zeros(eff, 'b')
        if len(fit_eff) == 0:
            logger.error('Trimmed zeros and nothing was left! Cannot fit!')
        else:
            fit_x = cut_values[:len(fit_eff)]
            # Fit on data with zeros removed, but plot with the full x spectrum
            pars, cov = curve_fit(f=exp_fit_fn, xdata=fit_x, ydata=fit_eff, p0=[0, 0, 0, 0], bounds=(-np.inf, np.inf))
            fitline = ax.plot(cut_values, exp_fit_fn(cut_values, *pars), alpha=.5, color=color)
            if include_percentile_calc:
                x_guess = cut_values[np.argmin(np.abs(fit_eff - 0.98))]
                logger.info('Fitting for eff @ 0.98; starting at initial guess x_guess=%d', x_guess)
                cut98 = inverse_exp_fit_fn(0.98, *pars, x_guess=x_guess)[0]
                ax.plot([cut98], [0.98], marker='o', markersize=10., color=color)
                ax.text(
                    cut98, 0.97, '{:.1f}'.format(cut98),
                    verticalalignment='top', horizontalalignment='left', fontsize=16
                    )
    if include_legend:
        if triggers:
            ax.text(
                0.97, 0.03,
                '\n'.join(triggers),
                verticalalignment='bottom', horizontalalignment='right',
                transform = ax.transAxes,
                fontsize=13
                )        
        ax.legend(loc='upper left')
    return ax


# ________________________________________________________
# Trigger titles per year

triggers_2016 = [
    # AK8PFJet triggers
    'HLT_AK8PFJet450_v',
    # CaloJet
    'HLT_CaloJet500_NoJetID_v',
    # PFJet and PFHT
    'HLT_PFJet450_v',
    'HLT_PFHT800_v',
    'HLT_PFHT900_v',
    # Trim mass triggers
    'HLT_AK8PFJet360_TrimMass30',
    'HLT_AK8PFHT700_TrimR0p1PT0p03Mass50',
    'HLT_AK8DiPFJet280_200_TrimMass30BTagCSVp20',
    # MET triggers
    # 'HLT_PFHT300_PFMET110_v',
    # 'HLT_PFHT300_PFMET100_v',
    ]

triggers_2017 = [
    # AK8PFJet triggers
    'HLT_AK8PFJet500_v',
    'HLT_AK8PFJet550_v',
    # CaloJet
    'HLT_CaloJet500_NoJetID_v',
    'HLT_CaloJet550_NoJetID_v',
    # PFJet and PFHT
    # 'HLT_PFHT800_v',
    'HLT_PFHT1050_v',
    'HLT_PFJet500_v',
    'HLT_PFJet550_v',
    # Trim mass jetpt+HT
    'HLT_AK8PFHT800_TrimMass50_v',
    'HLT_AK8PFHT850_TrimMass50_v',
    'HLT_AK8PFHT900_TrimMass50_v',
    'HLT_AK8PFJet400_TrimMass30_v',
    'HLT_AK8PFJet420_TrimMass30_v',
    # MET triggers
    # 'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v',
    # 'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v',
    # 'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v',
    # 'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v',
    # 'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v',
    # 'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v',
    ]

triggers_2018 = [
    # AK8PFJet triggers
    'HLT_AK8PFJet500_v',
    'HLT_AK8PFJet550_v',
    # CaloJet
    'HLT_CaloJet500_NoJetID_v',
    'HLT_CaloJet550_NoJetID_v',
    # PFJet and PFHT
    'HLT_PFHT1050_v', # but, interestingly, not HLT_PFHT8**_v or HLT_PFHT9**_v, according to the .txt files at least
    'HLT_PFJet500_v',
    'HLT_PFJet550_v',
    # Trim mass jetpt+HT
    'HLT_AK8PFHT800_TrimMass50_v',
    'HLT_AK8PFHT850_TrimMass50_v',
    'HLT_AK8PFHT900_TrimMass50_v',
    'HLT_AK8PFJet400_TrimMass30_v',
    'HLT_AK8PFJet420_TrimMass30_v',
    # MET triggers
    # 'HLT_PFHT500_PFMET100_PFMHT100_IDTight_v',
    # 'HLT_PFHT500_PFMET110_PFMHT110_IDTight_v',
    # 'HLT_PFHT700_PFMET85_PFMHT85_IDTight_v',
    # 'HLT_PFHT700_PFMET95_PFMHT95_IDTight_v',
    # 'HLT_PFHT800_PFMET75_PFMHT75_IDTight_v',
    # 'HLT_PFHT800_PFMET85_PFMHT85_IDTight_v',
    ]

def get_triggers_for_year(year):
    return globals()['triggers_{}'.format(year)]

# ________________________________________________________
# Common cuts

def htselection(arrays, cut_value):
    return arrays[b'HT'] >= cut_value
def htcut(arrays, cut_value):
    return arrays[b'TriggerPass'][htselection(arrays, cut_value)]
htcut_values = np.linspace(0., 1250., 60)
ht_title = 'HT (GeV)'

def jetptselection(arrays, cut_value):
    if svjflatanalysis.arrayutils.is_nested(arrays):
        passes = (arrays[b'JetsAK15_leading'].pt > cut_value).any()
    else:
        passes = (arrays[b'JetsAK15_leading.fCoordinates.fPt'] > cut_value).any()
    return passes
def jetptcut(arrays, cut_value):
    return arrays[b'TriggerPass'][jetptselection(arrays, cut_value)]
jetptcut_values = np.linspace(0., 750., 60)
jetpt_title = r'$p^{jet}_{T}$ (GeV)'

def mtselection(arrays, cut_value):
    return (arrays[b'JetsAK15_MT'] > cut_value).any()
def mtcut(arrays, cut_value):
    return arrays[b'TriggerPass'][mtselection(arrays, cut_value)]
mtcut_values = np.linspace(0., 1500., 60)
mt_title = r'$M_{T}$ (GeV)'

def metselection(arrays, cut_value):
    return arrays[b'MET'] > cut_value
def metcut(arrays, cut_value):
    return arrays[b'TriggerPass'][metselection(arrays, cut_value)]
metcut_values = np.linspace(0., 700., 50)
met_title = r'$MET$ (GeV)'

def msdselection(arrays, cut_value):
    return (arrays[b'JetsAK15_softDropMass'] > cut_value).any()
def msdcut(arrays, cut_value):
    return arrays[b'TriggerPass'][msdselection(arrays, cut_value)]
msdcut_values = np.linspace(0., 700., 60)
msd_title = r'$M_{SD}$ (GeV)'

def get_cut_fn_and_vals(variable):
    """
    Returns the cut function, list of reasonable cut values, and a formatted title
    """
    cutfn = globals()['{}cut'.format(variable)]
    cutvals = globals()['{}cut_values'.format(variable)]
    cuttitle = globals()['{}_title'.format(variable)]
    return cutfn, cutvals, cuttitle

# ________________________________________________________
# Final plot production functions

def trigger_plots_for_year_signal(year, variable, signals=None, notebook=False):
    """
    """
    svjflatanalysis.utils.reset_color()
    import re
    if not signals:
        logger.info('Recaching signals')
        signals = svjflatanalysis.samples.init_sigs_nohtcut(year, max_entries=None, progressbar=True)

    cut_function, cut_values, cut_title = get_cut_fn_and_vals(variable)

    ax = svjflatanalysis.utils.get_ax()

    for i_signal, signal in enumerate(signals):
        mass = int(re.search(r'\d+', signal.name).group())
        logger.info('Processing mz = %s', mass)
        eff, _, _ = svjflatanalysis.trigger.trigger_efficiency(signal, year, cut_function, cut_values)

        svjflatanalysis.trigger.plot_trigger_efficiency(
            cut_values, eff,
            label=r'$m_{{Z\prime}}={}$ GeV'.format(mass),
            cut_title=cut_title,
            include_curve_fit=(i_signal < 2), # Do limited number of curve fits to keep plot tidy
            include_percentile_calc=(i_signal==0), # Perc calc only for the first one (lowest mass)
            include_legend=(i_signal==len(signals)-1), # Legend only on last step
            ax=ax,
            triggers=svjflatanalysis.trigger.get_triggers_for_year(year)
            )

    ax = mplhep.cms.cmslabel(data=False, paper=False, year=str(year), ax=ax, fontsize=22)

    if notebook:
        return ax
    else:
        path = svjflatanalysis.utils.make_plotdir('trigger_plots_%d')
        plt.savefig(osp.join(path, '{}_{}.pdf'.format(year, variable)))
        plt.savefig(osp.join(path, '{}_{}.png'.format(year, variable)))


def trigger_plots_for_bkg(year, variable, bkgs=None, notebook=False):
    """
    """
    svjflatanalysis.utils.reset_color()
    import re
    if not bkgs:
        logger.info('Recaching bkgs')
        bkgs = svjflatanalysis.samples.init_bkgs_limited(progressbar=True)

    cut_function, cut_values, cut_title = get_cut_fn_and_vals(variable)

    ax = svjflatanalysis.utils.get_ax()

    bkg_cats = list(set([b.get_category() for b in bkgs]))
    bkg_cats.sort()
    for i_bkg, bkg_cat in enumerate(bkg_cats):
        bkg_datasets_this_cat = [b for b in bkgs if b.get_category() == bkg_cat]

        eff, _, _ = svjflatanalysis.trigger.trigger_efficiency(bkg_datasets_this_cat, year, cut_function, cut_values)

        svjflatanalysis.trigger.plot_trigger_efficiency(
            cut_values, eff,
            label=bkg_datasets_this_cat[0].get_title(),
            cut_title=cut_title,
            include_curve_fit=(i_bkg < 2), # Do limited number of curve fits to keep plot tidy
            include_percentile_calc=(i_bkg==0), # Perc calc only for the first one (lowest mass)
            include_legend=(i_bkg==len(bkg_cats)-1), # Legend only on last step
            ax=ax,
            triggers=svjflatanalysis.trigger.get_triggers_for_year(year)
            )

    ax = mplhep.cms.cmslabel(data=False, paper=False, year=str(year), ax=ax, fontsize=22)

    if notebook:
        return ax
    else:
        path = svjflatanalysis.utils.make_plotdir('trigger_plots_%d')
        plt.savefig(osp.join(path, 'bkg_{}_{}.pdf'.format(year, variable)))
        plt.savefig(osp.join(path, 'bkg_{}_{}.png'.format(year, variable)))
