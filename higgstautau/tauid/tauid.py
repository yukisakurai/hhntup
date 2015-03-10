import os
from rootpy.io import root_open as ropen
from .p851 import selection, nvtx_to_category, \
    CATEGORIES
from .common import LEVELS, PRONGS, nprong


"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsWinterConf2012
Values are percents / 100
"""
EFFIC_UNCERT_2011 = {
    'loose': {
        1: 0.04,
        3: 0.08,
    },
    'medium': {
        1: 0.05,
        3: 0.08,
    },
    'tight': {
        1: 0.04,
        3: 0.07,
    },
}

"""
https://twiki.cern.ch/twiki/pub/AtlasProtected/TauSystematicsWinterConf2012/tauWG_wtaunu2012winter_01.27.2012_soshi.pdf
"""
EFFIC_SF_2011 = {
    'loose': {
        1: (0.9641, 0.0390),
        3: (0.9415, 0.0755),
    },
    'medium': {
        1: (0.9210, 0.0421),
        3: (0.9849, 0.0800),
    },
    'tight': {
        1: (0.9933, 0.0424),
        3: (0.9711, 0.0710),
    },
}

"""
https://twiki.cern.ch/twiki/bin/viewauth/AtlasProtected/TauSystematicsSummerConf2012#Systematics_for_Tau_identificati
"""
EFFIC_SF_2012 = {
    'loose': {
        1: {'central': (1.008, 0.030),
            'forward': (1.012, 0.038)},
        3: {'central': (1.084, 0.059),
            'forward': (1.07, 0.11)},
    },
    'medium': {
        1: {'central': (0.992, 0.030),
            'forward': (0.952, 0.040)},
        3: {'central': (1.073, 0.071),
            'forward': (0.99, 0.10)},
    },
    'tight': {
        1: {'central': (0.985, 0.045),
            'forward': (0.876, 0.045)},
        3: {'central': (1.074, 0.097),
            'forward': (0.99, 0.10)},
    },
}


def effic_sf_uncert(tau, year):

    year = year % 1000
    if tau.JetBDTSigTight:
        wp = 'tight'
    elif tau.JetBDTSigMedium:
        wp = 'medium'
    elif tau.JetBDTSigLoose:
        wp = 'loose'
    else:
        wp = 'loose'
        print ("Warning: requested efficiency for tau below BDT loose. "
               "Using BDT loose.")
    return effic_sf_uncert_exc(wp, tau, year)


def effic_sf_uncert_exc(wp, tau, year):

    year = year % 1000
    np = nprong(tau.numTrack)
    if year == 11:
        sf, sf_uncert = EFFIC_SF_2011[wp][np]
    elif year == 12:
        if abs(tau.eta) > 1.5:
            region = 'forward'
        else:
            region = 'central'
        sf, sf_uncert = EFFIC_SF_2012[wp][np][region]
    else:
        raise ValueError(
                'no efficiency scale factors defined for year %d' % year)
    return sf, sf_uncert


HERE = os.path.dirname(os.path.abspath(__file__))

# uncertainty currently only valid for 2011 MC
BDT_UNCERT = {}
with ropen(os.path.join(HERE, 'bdt_uncertainty.root')) as f:
    for level in LEVELS.keys():
        BDT_UNCERT[level] = {}
        for prong in PRONGS:
            BDT_UNCERT[level][prong] = {}
            for category in CATEGORIES.keys():
                BDT_UNCERT[level][prong][category] = {}
                for dir in ('high', 'low'):
                    BDT_UNCERT[level][prong][category][dir] = f.Get(
                            '%s_%s_%dp_%s' % (
                                level, dir, prong, category)).Clone()


def uncertainty(score, pt, prong, nvtx):

    prong = nprong(prong)
    loose = selection('loose', prong, nvtx).Eval(pt)
    medium = selection('medium', prong, nvtx).Eval(pt)
    tight = selection('tight', prong, nvtx).Eval(pt)

    """
    if score < loose:
        print score, loose
        raise ValueError(
            'No uncertainties defined for scores lower than loose')
    """

    high_loose, low_loose = selection_uncertainty('loose', pt, prong, nvtx)
    high_medium, low_medium = selection_uncertainty('medium', pt, prong, nvtx)
    high_tight, low_tight = selection_uncertainty('tight', pt, prong, nvtx)

    b_l_high = (1. - loose) + high_loose
    b_m_high = (1. - medium) + high_medium
    b_t_high = (1. - tight) + high_tight

    b_l_low = (1. - loose) - low_loose
    b_m_low = (1. - medium) - low_medium
    b_t_low = (1. - tight) - low_tight

    if b_l_low <= 0:
        raise ValueError("low BDT loose selection error too high")
    if b_m_low <= 0:
        raise ValueError("low BDT medium selection error too high")
    if b_t_low <= 0:
        raise ValueError("low BDT tight selection error too high")

    if score > tight:
        # parabolic approx
        a = high_tight / (1. - tight)**4
        x = score - tight
        dx_high = high_tight - (a * x**4)
        if score > 1. - high_tight:
            # protect against limiting case
            # so we don't end up with scores above 1.
            if 1. - score < dx_high:
                dx_high = 1. - score
    elif score > tight - high_tight:
        #dx_high = (1. - score) * high_tight / b_t_high
        dx_high = high_tight
    elif score > medium - high_medium:
        dx_high = high_medium - (high_medium - high_tight) * (score - medium +
                high_medium) / (tight - high_tight - (medium - high_medium))
    elif score > loose - high_loose:
        dx_high = high_loose - (high_loose - high_medium) * (score - loose +
                high_loose) / (medium - high_medium - (loose - high_loose))
    else:
        dx_high = high_loose

    if score > tight:
        dx_low = low_tight
    elif score > medium + low_medium:
        dx_low = low_medium - (low_medium - low_tight) * (score - (medium +
                low_medium)) / (tight - (medium + low_medium))
    elif score > loose + low_loose:
        dx_low = low_loose - (low_loose - low_medium) * (score - (loose +
            low_loose)) / (medium + low_medium - (loose + low_loose))
    else:
        dx_low = low_loose

    return score + dx_high, score - dx_low


def selection_uncertainty(level, pt, prong, nvtx):

    uncert = BDT_UNCERT[level][prong][nvtx_to_category(nvtx)]
    return -1 * uncert['high'].Eval(pt), uncert['low'].Eval(pt)
