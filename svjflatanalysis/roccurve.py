import svjflatanalysis
import numpy as np, math, uuid

DEFAULT_FONTSIZE = 18

def get_eff(arrays_iterator, cut_function, cut_values):
    """
    Loops over arrays in the arrays_iterator and evaluates the cut_function
    at the cut_values.
    Returns a list of efficiences, passed events/objects, and total events/objects.
    cut_function is expected to return a tuple (n_pass, n_total) with input (arrays, cut_value).
    """
    n_cuts = len(cut_values)
    n_total = np.zeros(n_cuts)
    n_pass = np.zeros(n_cuts)
    for arrays, dataset in arrays_iterator:
        weight = dataset.get_weight()
        for i_cut, cut in enumerate(cut_values):
            this_n_pass, this_n_total = cut_function(arrays, cut)
            n_total[i_cut] += weight * this_n_total
            n_pass[i_cut] += weight * this_n_pass
    # Basically n_pass / n_total, but returns 0 if n_total has a 0 somewhere
    eff = np.divide(n_pass, n_total, out=np.zeros_like(n_pass), where=n_total!=0)
    return eff, n_pass, n_total
    
def roccurve(signals, bkgs, cut_function, cut_values):
    """
    Expects a list of signals and a list of bkgs (Dataset objects),
    and a cut_function and cut_values.
    """
    eff_sig, n_pass_sig, n_total_sig = get_eff(svjflatanalysis.iterate(signals), cut_function, cut_values)
    eff_bkg, n_pass_bkg, n_total_bkg = get_eff(svjflatanalysis.iterate(bkgs), cut_function, cut_values)
    return eff_sig, eff_bkg, n_pass_sig, n_pass_bkg, n_total_sig, n_total_bkg

def _draw_roccurve(eff_sig, eff_bkg, cut_values, ax, plot_cut_values=True):
    # Draw the actual roccurve
    line, = ax.plot(eff_bkg, eff_sig, marker='o')
    # Draw the cut_values at the specified points
    ndec = lambda number: 0 if abs(number) >= 100. else ( 1 if abs(number) >= 10. else 2 )
    if plot_cut_values:
        for i in range(len(cut_values)):
            ax.text(
                eff_bkg[i], eff_sig[i],
                '{0:.{ndec}f}'.format(cut_values[i], ndec=ndec(cut_values[i])),
                fontsize=14
                )
    return line

def plot_roccurve(signals, bkgs, cut_function, cut_values, ax):
    """
    Basic plotting style for a single roccurve, based on multiple signal and bkgs samples.
    Expects an ax object to be given, this function is not stand-alone
    """
    eff_sig, eff_bkg, n_pass_sig, n_pass_bkg, n_total_sig, n_total_bkg = roccurve(signals, bkgs, cut_function, cut_values)
    return _draw_roccurve(eff_sig, eff_bkg, cut_values, ax)

def plot_single_roccurve(signals, bkgs, cut_function, cut_values, ax=None):
    """
    Main routine for plotting a single roccurve
    """
    # Get a default ax if none is given
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()
    # Plot the base line
    ax.plot([0.0,1.0], [0.0,1.0], linestyle='--', color='xkcd:gray')
    # Plot the single roccurve
    line = plot_roccurve(signals, bkgs, cut_function, cut_values, ax=ax)
    line.set_label(bkgs[0].get_category())
    # Plot settings
    ax.set_xlim(0.0, 1.05)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel('Signal eff.', fontsize=DEFAULT_FONTSIZE)
    ax.set_xlabel('Bkg eff.', fontsize=DEFAULT_FONTSIZE)
    ax.legend(fontsize=DEFAULT_FONTSIZE)
    return ax

def plot_roccurves_per_bkg(signals, bkgs, cut_function, cut_values, ax=None):
    """
    Plots the roccurve per background category.
    Assumes signals are all datasets of the same signal.
    """
    # Get a default ax if none is given
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()
    # Get signal efficieny once
    eff_sig, n_pass_sig, n_total_sig = get_eff(svjflatanalysis.iterate(signals), cut_function, cut_values)
    # Perform some basic plotting setup
    ax.plot([0.0,1.0], [0.0,1.0], linestyle='--', color='xkcd:gray')
    ax.set_xlim(0.0, 1.05)
    ax.set_ylim(0.0, 1.05)
    ax.set_ylabel('Signal eff.', fontsize=DEFAULT_FONTSIZE)
    ax.set_xlabel('Bkg eff.', fontsize=DEFAULT_FONTSIZE)
    ax.legend(fontsize=DEFAULT_FONTSIZE)
    # Then efficiencies per bkg category (ttjets, qcd, ...)
    bkg_categories = list(set([ b.get_category() for b in bkgs ]))
    bkg_categories.sort()
    lines = {}
    for bkg_cat in bkg_categories:
        # Get Datasets that have this category
        bkgs_this_cat = [ b for b in bkgs if b.get_category() == bkg_cat ]
        # Compute efficiency in this category
        eff_bkg, n_pass_bkg, n_total_bkg = get_eff(svjflatanalysis.iterate(bkgs_this_cat), cut_function, cut_values)
        # Draw roccurve for this category
        line = _draw_roccurve(eff_sig, eff_bkg, cut_values, ax)
        line.set_label(bkg_cat)
        # Save this line in a dict for potential outputting/modifying
        lines[bkg_cat] = line
    return ax

def plot_multiple_signals_per_bkg(signals, bkgs, cut_function, cut_values, title=None):
    # Get all the categories that are going to be plotted
    sig_categories = list(set([ s.get_category() for s in signals ]))
    sig_categories.sort()
    bkg_categories = list(set([ b.get_category() for b in bkgs ]))
    bkg_categories.sort()
    # Setup the figure
    import matplotlib.pyplot as plt
    n_sigs = len(sig_categories)
    ncols = min(2, n_sigs)
    nrows = math.ceil(n_sigs/2.)
    fig = plt.figure(figsize=(16,8*nrows))
    # Compute efficiency per category once
    sig_effs = {}
    for cat in sig_categories:
        sigs_this_cat = [ s for s in signals if s.get_category() == cat]
        eff_sig, n_pass_sig, n_total_sig = get_eff(svjflatanalysis.iterate(sigs_this_cat), cut_function, cut_values)
        sig_effs[cat] = eff_sig
    bkg_effs = {}
    for cat in bkg_categories:
        bkgs_this_cat = [ b for b in bkgs if b.get_category() == cat]
        eff_bkg, n_pass_bkg, n_total_bkg = get_eff(svjflatanalysis.iterate(bkgs_this_cat), cut_function, cut_values)
        bkg_effs[cat] = eff_bkg
    # Do the plots in subplots, 1 per signal category
    lines = []
    for i, sig_cat in enumerate(sig_categories):
        sigs_this_cat = [ s for s in signals if s.get_category() == cat]
        ax = fig.add_subplot(nrows, ncols, i+1)
        ax.plot([0.0,1.0], [0.0,1.0], linestyle='--', color='xkcd:gray')
        ax.set_xlim(0.0, 1.05)
        ax.set_ylim(0.0, 1.05)
        ax.set_ylabel('Signal eff.', fontsize=DEFAULT_FONTSIZE)
        ax.set_xlabel('Bkg eff.', fontsize=DEFAULT_FONTSIZE)
        ax.set_title(sig_cat, fontsize=DEFAULT_FONTSIZE)
        for bkg_cat in bkg_categories:
            # print(sig_cat, bkg_cat, sig_effs[sig_cat], bkg_effs[bkg_cat])
            line = _draw_roccurve(sig_effs[sig_cat], bkg_effs[bkg_cat], cut_values, ax)
            line.set_label(bkg_cat)
            lines.append(line)
        ax.legend(fontsize=DEFAULT_FONTSIZE)
    if title: fig.suptitle(title, y=0.92, fontsize=DEFAULT_FONTSIZE + 3)
    return fig, lines

def plot_roc_single_sig_and_bkg_multiple_fns(signal, bkgs, cut_functions, cut_values, ax=None, titles=None):
    if ax is None: ax = svjflatanalysis.utils.get_ax()
    colors = svjflatanalysis.utils.xkcd_color_wheel(47)
    for i in range(len(cut_functions)):
        fn = cut_functions[i]
        vals = cut_values[i]
        title = get_title(fn) if titles is None else titles[i]
        # color = svjflatanalysis.utils.get_default_color()
        color = colors[i]
        eff_sig, n_pass_sig, n_total_sig = get_eff(svjflatanalysis.iterate(signal), fn, vals)
        eff_bkg, n_pass_bkg, n_total_bkg = get_eff(svjflatanalysis.iterate(bkgs), fn, vals)        
        line = _draw_roccurve(eff_sig, eff_bkg, vals, ax, plot_cut_values=False)
        line.set_label(title)
        line.set_color(color)
        i_half = len(eff_sig) // 2
    ax.set_title('Sig {} / Bkg {}'.format(signal.get_category(), bkgs[0].get_category()))
    ax.legend(fontsize=10, ncol=2)
    return ax

def hist_single_distribution(
    arrays_iterator, get_array,
    varname='somevar', vartitle=None, distrname='somedistr', distrtitle=None,
    hist=None, left=-1., right=1., nbins=50
    ):
    """
    Fills a coffea.hist.Hist for a single distribution.
    Takes a list of Dataset objects, and a function `get_array` that should return a 
    numpy-like array when given an arrays object.
    Also requires a string `name` to know in which hist to fill it
    """
    if hist is None:
        import coffea.hist
        vartitle = varname if vartitle is None else vartitle
        hist = coffea.hist.Hist(
            "Count",
            coffea.hist.Bin(varname, vartitle, nbins, left, right),
            coffea.hist.Cat('label', varname),
            )
    for arrays, dataset in arrays_iterator:
        print(dataset.get_weight(), get_array(arrays))
        hist.fill(label=distrname, weight=dataset.get_weight(), **{varname: get_array(arrays)})
    return hist

def plot_single_distribution(arrays_iterator, get_array, **kwargs):
    ax = kwargs.pop('ax', None)
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()
    plot_options = kwargs.pop('plot_options', {})
    hist = kwargs.get('hist', hist_single_distribution(arrays_iterator, get_array, **kwargs))
    import coffea.hist
    coffea.hist.plot1d(hist, ax=ax, **plot_options)
    return ax, hist


class Hist(object):
    """
    Wrapper for coffea.hist.Hist for multiple distributions of a single variable
    """
    def __init__(self, **kwargs):
        self.i_fill = 0
        self.set(**kwargs)

    def set(self, **kwargs):
        self.varname = kwargs.get('varname', 'somevar')
        self.vartitle = kwargs.get('vartitle', self.varname)
        self.ylabel = kwargs.pop('ylabel', 'Count')
        self.catname = kwargs.pop('catname', 'category')
        import coffea.hist
        self.hist = coffea.hist.Hist(
            self.ylabel,
            coffea.hist.Bin(self.varname, self.vartitle, kwargs.pop('nbins', 100), kwargs.pop('left', 0.), kwargs.pop('right', 100.)),
            coffea.hist.Cat('category', self.vartitle),
            )

    def fill(self, arrays_iterator, get_array, cat=None):
        self.i_fill += 1
        if cat is None: cat = 'auto{}'.format(self.i_fill)
        for arrays, dataset in arrays_iterator:
            self.hist.fill(category=cat, weight=dataset.get_weight(), **{self.varname: get_array(arrays)})

    def plot(self, ax=None, **plot_options):
        if ax is None: ax = svjflatanalysis.utils.get_ax()
        import coffea.hist
        coffea.hist.plot1d(self.hist, ax=ax, **plot_options)



def plot_multiple_distributions(arrays_iterators_dict, get_array, **kwargs):
    """
    Expectes arrays_iterators_dict = {'signal' : arrays_iterator, 'bkg' : ...}
    """
    import coffea.hist
    ax = kwargs.pop('ax', None)
    if ax is None:
        import matplotlib.pyplot as plt
        fig = plt.figure(figsize=(8,8))
        ax = fig.gca()

    hist = kwargs.pop('hist', None)
    if hist is None:
        varname = kwargs.get('varname', 'somevar')
        vartitle = kwargs.get('vartitle', varname)
        hist = coffea.hist.Hist(
            "Count",
            coffea.hist.Bin(varname, vartitle, kwargs.pop('nbins', 100), kwargs.pop('left', 0.), kwargs.pop('right', 100.)),
            coffea.hist.Cat('label', vartitle),
            )

    plot_options = kwargs.pop('plot_options', {})
    for distrname, arrays_iterator in arrays_iterators_dict.items():
        plot_single_distribution(
            arrays_iterator, get_array,
            distrname=distrname, ax=ax, hist=hist,
            plot_options=plot_options,
            **kwargs
            )

    return ax, hist


# ____________________________________________________________________
# Common ROC variables and cut values

def msd_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_softDropMass'] > cut).any().sum()
    n_total = (arrays[b'JetsAK15_softDropMass'] > 0.0).any().sum()
    return n_pass, n_total
msd_cut_values = [ 0.0, 10., 30., 50., 70., 100., 150., 200., 500., 1e7 ]

def msd_leading_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_leading_softDropMass'] > cut).any().sum()
    n_total = (arrays[b'JetsAK15_leading_softDropMass'] > 0.0).any().sum()
    return n_pass, n_total
msd_leading_cut_values = [ 0.0, 10., 30., 50., 70., 100., 150., 200., 500., 1e7 ]

def msd_subleading_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_subleading_softDropMass'] > cut).any().sum()
    n_total = (arrays[b'JetsAK15_subleading_softDropMass'] > 0.0).any().sum()
    return n_pass, n_total
msd_subleading_cut_values = [ 0.0, 10., 30., 50., 70., 100., 150., 200., 500., 1e7 ]

def met_cut_function(arrays, cut):
    n_pass = (arrays[b'MET'] > cut).sum()
    n_total = (arrays[b'MET'] > 0.0).sum()
    return n_pass, n_total
met_cut_values = [ 0.0, 1.0, 10., 20., 40., 60., 80., 100., 150., 200., 500., 1000, 2000., 1e7 ]

def dphimet_cut_function(arrays, cut):
    from math import pi
    try:
        phi = arrays[b'JetsAK15_subleading'].phi
    except AttributeError:
        phi = arrays[b'JetsAK15_subleading.fCoordinates.fPhi']
    dphi = np.abs(phi - arrays[b'METPhi'])
    dphi[dphi > 2.*pi] = dphi[dphi > 2.*pi] - 2.*pi  # Whole circles subtracted
    dphi[dphi > pi] = 2.*pi - dphi[dphi > pi]  # Pick the smaller angle
    n_total = dphi.flatten().shape[0]
    selection = (dphi < cut)
    n_pass = selection.flatten().sum()
    if cut == 1e6 and n_pass < n_total:
        print(dphi)
        print(dphi[dphi>=cut].flatten())
        print(dphi.flatten().shape)
        print(n_pass)
        print(n_total)
    return n_pass, n_total
dphimet_cut_function.operator = '>'
dphimet_cut_values = [ 0.0, 0.5, 0.8, 1.0, 1.5, 2.5, 1e6]

def deltaeta_cut_function(arrays, cut):
    leading_jets = svjflatanalysis.arrayutils.Jets(b'JetsAK15_leading', arrays, has_subjets=False, has_substructure=False)
    subleading_jets = svjflatanalysis.arrayutils.Jets(b'JetsAK15_subleading', arrays, has_subjets=False, has_substructure=False)
    n_total = subleading_jets.counts.sum()
    # Only get the leading jets for when there is a subleading jet too
    leading_jets = leading_jets[subleading_jets.counts > 0]
    subleading_jets = subleading_jets[subleading_jets.counts > 0]
    deta = np.abs(leading_jets.eta - subleading_jets.eta).flatten()
    n_pass = (deta < cut).sum()
    return n_pass, n_total
deltaeta_cut_function.operator = '<'
deltaeta_cut_values = [ 0.0, 0.15, 0.30, 0.50, 1.0, 1.5, 2.0, 3.0, 10.]

def rt_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_subleading_RT'] > cut).flatten().sum()
    n_total = (arrays[b'JetsAK15_subleading_RT'] > 0.0).flatten().sum()
    return n_pass, n_total
rt_cut_function.operator = '>'
rt_cut_values = [ 0.0, 0.15, 0.30, 0.50, 1.0, 1.5, 10.]

def rtnew_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_subleading_RTnew'] > cut).flatten().sum()
    n_total = (arrays[b'JetsAK15_subleading_RTnew'] > 0.0).flatten().sum()
    return n_pass, n_total
rtnew_cut_function.operator = '>'
rtnew_cut_values = [ 0.0, 0.15, 0.30, 0.50, 1.0, 1.5, 10.]

def mt_cut_function(arrays, cut):
    n_pass = (arrays[b'JetsAK15_subleading_MT'] > cut).flatten().sum()
    n_total = (arrays[b'JetsAK15_subleading_MT'] > 0.0).flatten().sum()
    return n_pass, n_total
mt_cut_function.operator = '>'
mt_cut_values = [ 0.0, 10., 30., 50., 70., 100., 150., 200., 500., 1e7 ]

def subleading_tau21_cut_function(arrays, cut):
    r = svjflatanalysis.arrayutils.safe_divide(
        arrays[b'JetsAK15_subleading_NsubjettinessTau2'].flatten(),
        arrays[b'JetsAK15_subleading_NsubjettinessTau1'].flatten()
        )
    n_pass = (r < cut).flatten().sum()
    n_total = (r > 0.0).flatten().sum()
    return n_pass, n_total
subleading_tau21_cut_function.operator = '<'

def subleading_tau32_cut_function(arrays, cut):
    r = svjflatanalysis.arrayutils.safe_divide(
        arrays[b'JetsAK15_subleading_NsubjettinessTau3'].flatten(),
        arrays[b'JetsAK15_subleading_NsubjettinessTau2'].flatten()
        )
    n_pass = (r > cut).flatten().sum()
    n_total = (r > 0.0).flatten().sum()
    return n_pass, n_total
subleading_tau32_cut_function.operator = '<'

def subleading_tau31_cut_function(arrays, cut):
    r = svjflatanalysis.arrayutils.safe_divide(
        arrays[b'JetsAK15_subleading_NsubjettinessTau3'].flatten(),
        arrays[b'JetsAK15_subleading_NsubjettinessTau1'].flatten()
        )
    n_pass = (r < cut).flatten().sum()
    n_total = (r > 0.0).flatten().sum()
    return n_pass, n_total
subleading_tau31_cut_function.operator = '<'


def generic_jetvar_cut_function_factory(varname, operator='larger'):
    def cutfn(arrays, cut):
        if operator == 'larger':
            n_pass = (arrays[varname] > cut).flatten().sum()
            n_total = (arrays[varname] > 0.0).flatten().sum()
        else:
            n_pass = (arrays[varname] < cut).flatten().sum()
            n_total = (arrays[varname] < 1e7).flatten().sum()            
        return n_pass, n_total
    cutfn.name = varname.decode('utf-8')
    return cutfn
subleading_pt_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading.fCoordinates.fPt')
subleading_pt_cut_function.name = 'pT_subleading'
subleading_axismajor_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_axismajor')
subleading_axisminor_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_axisminor')
subleading_ecfN2b1_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_ecfN2b1')
subleading_ecfN2b2_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_ecfN2b2')
subleading_ecfN3b1_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_ecfN3b1')
subleading_ecfN3b2_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_ecfN3b2')
subleading_girth_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_girth')
subleading_ptD_cut_function = generic_jetvar_cut_function_factory(b'JetsAK15_subleading_ptD', operator='smaller')


def apply_trigger_first(cut_fn):
    """
    Decorator for post-trigger cuts
    """
    def wrapped(arrays, cut):
        arrays = svjflatanalysis.arrayutils.apply_trigger_and_jetpt550(arrays, 2018)
        return cut_fn(arrays, cut)
    return wrapped

def apply_window_cutfn(cutfn, left, right):
    """
    """
    def cutfn_with_mass_window(arrays, cut):
        arrays = svjflatanalysis.arrayutils.apply_window(arrays, left, right)
        return cutfn(arrays, cut)
    if hasattr(cutfn, 'operator'): cutfn_with_mass_window.operator = cutfn.operator
    cutfn_with_mass_window.name = cutfn.name if hasattr(cutfn, 'name') else cutfn.__name__
    cutfn_with_mass_window.left = left
    cutfn_with_mass_window.right = right
    return cutfn_with_mass_window

def apply_window_onfn(fn, left, right):
    """
    """
    def fn_with_mass_window(arrays, *args, **kwargs):
        arrays = svjflatanalysis.arrayutils.apply_window(arrays, left, right)
        return fn(arrays, *args, **kwargs)
    if hasattr(fn, 'operator'): fn_with_mass_window.operator = fn.operator
    fn_with_mass_window.name = fn.__name__
    fn_with_mass_window.left = left
    fn_with_mass_window.right = right
    return fn_with_mass_window


def get_title(fn):
    """
    Takes a cut function and tries to return a title for it
    """
    title = fn.name if hasattr(fn, 'name') else fn.__name__
    title = title.replace('_cut_function','')
    suffix = []
    # if 'JetsAK15_subleading_' in title:
    #     suffix.append(r'$j^{\mathrm{AK15}}_{\mathrm{subl}}$')
    title = title.replace('JetsAK15_subleading_', '').replace('subleading_', '')
    if hasattr(fn, 'left'):
        suffix.append('({:.0f} < {} < {:.0f})'.format(fn.left, svjflatanalysis.utils.get_title('mt'), fn.right))
    # Transform variable name to title stirng
    title = svjflatanalysis.utils.get_title(title)
    if hasattr(fn, 'operator'):
        title += ' ' + fn.operator + ' cut'
    # Add the suffix
    title += ' ' + ' '.join(suffix)
    return title
