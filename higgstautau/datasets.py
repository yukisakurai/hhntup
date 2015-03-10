"""
This module generates a database of all MC and data datasets
"""
from rootpy.io import root_open, DoesNotExist

#from multiprocessing import Pool, cpu_count

import sys
from operator import itemgetter
import logging
import re
import glob
import os
import cPickle as pickle
import atexit
import fnmatch
from collections import namedtuple

import yaml

from . import log; log = log[__name__]
from .decorators import cached_property
from .yaml_utils import Serializable
from . import xsec

USE_PYAMI = True
try:
    from pyAMI.client import AMIClient
    from pyAMI.query import get_dataset_xsec_effic, \
                            get_dataset_info, \
                            get_datasets, \
                            get_provenance, \
                            get_periods, \
                            get_runs
    from pyAMI import query
    from pyAMI.auth import AMI_CONFIG, create_auth_config
except ImportError:
    USE_PYAMI = False
    log.warning("pyAMI is not installed. "
                "Cross section retrieval will be disabled.")

# data types
DATA, MC, EMBED, MCEMBED = range(4)

TYPES = {
    'DATA':    DATA,
    'MC':      MC,
    'EMBED':   EMBED,
    'MCEMBED': MCEMBED,
}

Namedset = namedtuple('Namedset',
                      'name tags meta properties')
Dataset = namedtuple('Dataset',
                     Namedset._fields + ('datatype',))


class Fileset(namedtuple('Fileset', Dataset._fields + ('files', 'treename'))):

    def split(self, partitions):
        files = self.files[:]
        fileset_files = [[] for _ in xrange(partitions)]
        while len(files) > 0:
            for fileset in fileset_files:
                if len(files) > 0:
                    fileset.append(files.pop(0))
                else:
                    break
        mydict = self._asdict()
        filesets = []
        for fileset in fileset_files:
            mydict['files'] = fileset
            filesets.append(Fileset(**mydict))
        return filesets


class Treeset(namedtuple('Treeset', Dataset._fields + ('trees',))):

    def GetEntries(self, *args, **kwargs):
        return sum([tree.GetEntries(*args, **kwargs) for tree in self.trees])

    def Scale(self, value):
        for tree in self.trees:
            tree.Scale(value)

    def __iter__(self):
        for tree in self.trees:
            yield tree

    def Draw(self, *args, **kwargs):
        for tree in self.trees:
            tree.Draw(*args, **kwargs)

ATLASFileset = namedtuple('ATLASFileset', Fileset._fields + ('year', 'grl',))

DS_PATTERN = re.compile(
    '^(?P<prefix>\S+\.)?'
    '(?P<type>(data|mc))(?P<year>\d+)_(?P<energy>\d+)TeV'
    '\.(?P<id>(\d+|period[A-Z]))'
    '\.(?P<name>\w+)'
    '(\.PhysCont)?'
    '(\.(?P<ntup>merge\.NTUP_TAU(MEDIUM)?))?'
    '\.(?P<tag>\w+)'
    '(\.small)?'
    '(\.v(?P<version>\d+))?(_s)?'
    '\.(?P<suffix>\S+)$')

MC_TAG_PATTERN1 = re.compile(
    '^e(?P<evnt>\d+)_'
    's(?P<digi>\d+)_'
    's(?P<digimerge>\d+)_'
    'r(?P<reco>\d+)_'
    'r(?P<recomerge>\d+)_'
    'p(?P<ntup>\d+)$')

# not all valid samples have a recomerge tag:
MC_TAG_PATTERN2 = re.compile(
    '^e(?P<evnt>\d+)_'
    '[sa](?P<digi>\d+)_'
    '[sa](?P<digimerge>\d+)_'
    'r(?P<reco>\d+)_'
    'p(?P<ntup>\d+)$')

# Embedded sample pattern
EMBED_PATTERN11 = re.compile(
    '^(?P<prefix>\S+)?'
    'period(?P<period>[A-Z])'
    '\.DESD_SGLMU'
    '\.pro(?P<prod>\d+)'
    '\.embedding-(?P<embedtag>\S+)?'
    '\.Ztautau_'
    '(?P<channel>(lh)|(hh))_'
    '(?P<isol>[a-z]+)_'
    '(?P<mfs>[a-z]+)_'
    'rereco_'
    'p(?P<tag>\d+)_'
    'EXT0'
    '(\.(?P<suffix>\S+))?$')

EMBED_PATTERN12 = re.compile(
    '^(?P<prefix>\S+)?'
    'period(?P<period>[A-Z])'
    '\.DESD_ZMUMU'
    '\.pro(?P<prod>\d+)'
    '\.embedding-(?P<embedtag>\S+)?'
    '\.Ztautau_'
    '(?P<channel>(lh)|(hh))_'
    '(((high)|(low))pt_)?'
    '(?P<mfs>[a-z]+)_'
    'filter_'
    'taureco_'
    'p(?P<tag>\d+)_'
    'EXT0'
    '(\.(?P<suffix>\S+))?$')

EMBED_PATTERN12_NEW = re.compile(
    '^(?P<prefix>\S+)?'
    'data12_8TeV\.'
    'period(?P<period>[A-Z])\.'
    'physics_Muons\.PhysCont\.'
    'NTUP_EMB(?P<channel>(LH)|(HH))'
    '(?P<sys>(DN)|(IM)|(UP))\.'
    '(?P<suffix>\S+)')

MC_EMBED_PATTERN = re.compile(
    '^(?P<prefix>\S+)?'
    'Pyth8.DESD_SGLMU.pro14.embedding-01-01-10.'
    'Ztautau_MCEmbedding[\d]*_hh(?P<sys>(dn)|(up))?_p1344_EXT0'
    '(\.(?P<suffix>\S+))?$')

## Common lephad ntuple pattern
CN_MC_PATTERN12 = re.compile(
    '^(?P<prefix>\S+\.)?'
    '(?P<id>\d+)'
    '\.(?P<name>\w+)'
    '\.(?P<tag>\w+)'
    '_lhCN'
    '(v(?P<version1>\d+))?'
    '(-(?P<version2>\d+))?'
    '(-(?P<version3>\d+))?'
    '\.(?P<suffix>\S+)$')

CN_DATA_PATTERN12 = re.compile(
    '^(?P<prefix>\S+\.)?'
    'data12_8TeV\.'
    '(?P<id>\S+)'
    '\.(?P<name>\w+)'
    '((\.TunaCont.2013-March-29.v03)?)'
    '\.(?P<tag>\w+)'
    '_lhCN'
    '(v(?P<version1>\d+))?'
    '(-(?P<version2>\d+))?'
    '\.(?P<suffix>\S+)$')

CN_EMBED_PATTERN12 = re.compile(
    '^(?P<prefix>\S+\.)?'
    'data12_8TeV\.'
    '(?P<id>\S+)'
    '\.(?P<name>\w+)'
    '\.PhysCont'
    '((\.NTUP_EMB)?)'
    '(?P<channel>(LH)|(HH))'
    '(?P<mfs>(IM)|(UP)|(DN))'
    '\.grp14_v02'
    '\_(?P<tag>\w+)'
    '_lhCN'
    '(v(?P<version1>\d+))?'
    '(-(?P<version2>\d+))?'
    '(-(?P<version3>\d+))?'
    '\.(?P<suffix>\S+)$')

# MC[11|12][a|b|c|...] categories are defined here
# Each MC dataset is automatically classified
# acccording to these categories by matching the reco
# and merge tags of the dataset name.
# Order by decreasing preference:
MC_CATEGORIES = {
    'mc11a': {'reco':  (2730, 2731),
              'merge': (2780, 2700)},
    'mc11b': {'reco':  (2920, 2923),
              'merge': (3063, 2993, 2900)},
    'mc11c': {'reco':  (3043, 3060, 3108),
              'merge': (3109, 3063, 2993)},
    'mc12a': {'reco':  (3753, 3752, 3658, 3605, 3553, 3542, 3549),
              'merge': (3549,)},
    'mc12b': {'reco':  (4485, 5470,),
              'merge': (4540,)}}

HERE = os.path.dirname(os.path.abspath(__file__))

# Any datasets which don't have the provenance stored properly in AMI
# should be hardcoded here (it happens)
DS_NOPROV = {}

# Cross-sections are cached so that we don't need to keep asking AMI
# for them over and over
XSEC_CACHE_FILE = os.path.join(HERE, 'xsec', 'cache.pickle')
XSEC_CACHE_MODIFIED = False
XSEC_CACHE = {}

if USE_PYAMI:
    amiclient = AMIClient()
    if not os.path.exists(AMI_CONFIG):
        create_auth_config()
    amiclient.read_config(AMI_CONFIG)


class NoMatchingDatasetsFound(Exception):
    pass


GLOBAL_BASE = '/global/'


def find_global(path):

    if not path.startswith('/global/'):
        raise ValueError("path must be absolute and rooted at /global")

    path = re.sub('^/global/', '/cluster/data%02d/export/', path)

    for node in range(1, 13):
        if os.path.exists(path % node):
            return path % node
    raise IOError('path %s does not exist' % path)


class Database(dict):

    @classmethod
    def match_to_ds(cls, match):
        """
        Construct the original NTUP dataset name from a skim match object
        """
        if match.group('year') == '11':
            ntup = 'merge.NTUP_TAUMEDIUM'
        else:
            ntup = 'merge.NTUP_TAU'
        return '%s%s_%sTeV.%s.%s.%s.%s' % (
                match.group('type'),
                match.group('year'),
                match.group('energy'),
                match.group('id'),
                match.group('name'),
                ntup,
                match.group('tag'))

    def __init__(self, name='datasets', verbose=False, stream=None):
        super(Database, self).__init__()
        self.name = name
        self.verbose = verbose
        self.filepath = os.path.join(HERE, '%s.yml' % self.name)
        if os.path.isfile(self.filepath):
            with open(self.filepath) as db:
                log.info("Loading database '%s' ..." % self.name)
                d = yaml.load(db)
                if d:
                    self.update(d)
        self.modified = False
        if stream is None:
            self.stream = sys.stdout
        else:
            self.stream = stream

    def write(self):
        if self.modified:
            with open(self.filepath, 'w') as db:
                log.info("Saving database '%s' ..." % self.name)
                yaml.dump(dict(self), db)

    def reset(self):
        return self.clear()

    def clear(self):
        # erase all datasets in database
        log.info("Resetting database '%s' ..." % self.name)
        super(Database, self).clear()
        self.modified = True

    def validate(self,
                 pattern=None,
                 datatype=None,
                 year=None):
        ds = {}
        for name, info in self.items():
            if year is not None and info.year != year:
                continue
            if datatype is not None and info.datatype != datatype:
                continue
            if info.datatype == DATA and info.id < 0:
                # only validate data run datasets
                continue
            if pattern is None or fnmatch.fnmatch(name, pattern):
                ds[name] = info
        incomplete = []
        for name, info in sorted(ds.items(), key=lambda item: item[0]):
            log.info("Validating %s ..." % name)
            complete = validate_single((name, info), child=False)
            log.info("Complete: %s" % complete)
            log.info('-' * 50)
            if not complete:
                incomplete.append(info.ds)
        #pool = Pool(processes=cpu_count())
        #for result, complete in pool.map(
        #        validate_single, sorted(ds.items(), key=itemgetter(0))):
        #    print result
        #    print "Complete: %s" % complete
        #    print '-'*50
        #    if not complete:
        #        all_complete = False
        if not incomplete:
            log.info("ALL DATASETS ARE COMPLETE")
        else:
            log.warning("SOME DATASETS ARE NOT COMPLETE:")
            for ds in incomplete:
                print ds

    def scan(self, year,
             mc_path=None,
             mc_prefix=None,
             mc_pattern=None,
             mc_treename=None,
             mc_sampletype=None,
             data_path=None,
             data_prefix=None,
             data_pattern=None,
             data_treename=None,
             data_sampletype=None,
             data_grl=None,
             data_period_containers=False,
             embed_path=None,
             embed_prefix=None,
             embed_pattern=None,
             embed_treename=None,
             embed_sampletype=None,
             versioned=False,
             deep=False):
        """
        Update the dataset database
        """
        log.info("Updating database '%s' ..." % self.name)
        self.modified = True

        ###############################
        # MC
        ###############################
        if mc_path is not None:
            if deep:
                mc_dirs = get_all_dirs_under(mc_path, prefix=mc_prefix)
            else:
                if mc_prefix:
                    mc_dirs = glob.glob(os.path.join(mc_path, mc_prefix) + '*')
                else:
                    mc_dirs = glob.glob(os.path.join(mc_path, '*'))

            for dir in mc_dirs:
                dirname, basename = os.path.split(dir)
                if mc_sampletype == 'standard':
                    match  = re.match(DS_PATTERN, basename)
                    if match:
                        if int(match.group('year')) != (year % 1E3):
                            continue
                        if match.group('type') != 'mc':
                            continue
                        ds_name = Database.match_to_ds(match)
                        name = match.group('name')
                        tag = match.group('tag')
                        try:
                            version = int(match.group('version'))
                        except IndexError:
                            version = 0
                        except:
                            log.warning(basename)
                            raise
                        tag_match = re.match(MC_TAG_PATTERN1, tag)
                        tag_match2 = re.match(MC_TAG_PATTERN2, tag)
                        MC_TAG_PATTERN = MC_TAG_PATTERN1

                        if (tag_match2 and not tag_match) :
                            tag_match = tag_match2
                            MC_TAG_PATTERN = MC_TAG_PATTERN2

                        if not tag_match:
                            log.warning("not tag-matched: %s" % basename)
                            continue
                        cat = None
                        for cat_name, cat_params in MC_CATEGORIES.items():
                            if int(tag_match.group('reco')) in cat_params['reco']:
                                cat = cat_name
                                break
                        if cat is None:
                            log.warning(
                                "does not match a category: %s" % basename)
                            continue
                        name += '.' + cat
                        dataset = self.get(name, None)
                        if dataset is not None and version == dataset.version:
                            if tag != dataset.tag:
                                this_reco = int(tag_match.group('reco'))
                                other_reco = int(

                                    re.match(dataset.tag_pattern,
                                             dataset.tag).group('reco'))
                                use_mergetag = True
                                try:
                                    this_merge = int(tag_match.group('recomerge'))
                                    other_merge = int(
                                        re.match(dataset.tag_pattern,
                                                 dataset.tag).group('recomerge'))
                                except IndexError:
                                    use_mergetag = False

                                cat_params = MC_CATEGORIES[cat]
                                reco_tags = list(cat_params['reco'])
                                merge_tags = list(cat_params['merge'])
                                assert(this_reco in reco_tags and other_reco in reco_tags)
                                take_this = False
                                if reco_tags.index(this_reco) < reco_tags.index(other_reco):
                                    take_this = True
                                elif (use_mergetag and this_reco == other_reco and
                                      (merge_tags.index(this_merge) <
                                       merge_tags.index(other_merge))):
                                    take_this = True

                                if take_this:
                                    log.warning("taking %s over %s" % (
                                        basename, dataset.ds))
                                    self[name] = Dataset(
                                        name=name,
                                        datatype=MC,
                                        treename=mc_treename,
                                        ds=ds_name,
                                        id=int(match.group('id')),
                                        category=cat,
                                        version=version,
                                        tag_pattern=MC_TAG_PATTERN.pattern,
                                        tag=tag,
                                        dirs=[dir],
                                        file_pattern=mc_pattern,
                                        year=year)
                            elif dir not in dataset.dirs:
                                dataset.dirs.append(dir)
                        elif dataset is None or (
                                dataset is not None and version > dataset.version):
                            self[name] = Dataset(
                                name=name,
                                datatype=MC,
                                treename=mc_treename,
                                ds=ds_name,
                                id=int(match.group('id')),
                                category=cat,
                                version=version,
                                tag_pattern=MC_TAG_PATTERN.pattern,
                                tag=tag,
                                dirs=[dir],
                                file_pattern=mc_pattern,
                                year=year)
                    elif self.verbose:
                        log.warning("not a valid mc dataset name: %s" % basename)

                elif mc_sampletype == 'lhCN':
                    match  = re.match(CN_MC_PATTERN12, basename)
                    if match:

                        name = match.group('name')
                        cat = 'mc12a'
                        tag = match.group('tag')
                        year = 2012

                        ## Calculate a version int
                        version_1 = match.group('version1')
                        version_2 = match.group('version2')
                        version = int(version_1)*1000 + int(version_2)*10

                        dataset = self.get(name, None)
                        if dataset is not None and version == dataset.version:
                            if dir not in dataset.dirs:
                                dataset.dirs.append(dir)
                        else:

                            log.info('\'%s\',' % name)
                            self[name] = Dataset(
                                name=name,
                                datatype=MC,
                                treename=mc_treename,
                                ds=name,
                                id=int(match.group('id')),
                                category=cat,
                                version=version,
                                tag_pattern=None,
                                tag=tag,
                                dirs=[dir],
                                file_pattern=mc_pattern,
                                year=year)


        #####################################
        # EMBEDDING
        #####################################
        if embed_path is not None:

            if deep:
                embed_dirs = get_all_dirs_under(embed_path, prefix=embed_prefix)
            else:
                if embed_prefix:
                    embed_dirs = glob.glob(
                        os.path.join(embed_path, embed_prefix) + '*')
                else:
                    embed_dirs = glob.glob(
                        os.path.join(embed_path, '*'))

            if embed_sampletype == 'new':

                EMBED_PATTERN = EMBED_PATTERN12_NEW

                # determine what channels are available
                channels = {}
                for dir in embed_dirs:
                    if os.path.isdir(dir):
                        dirname, basename = os.path.split(dir)
                        match = re.match(EMBED_PATTERN, basename)
                        if match:
                            channel = match.group('channel')
                            if channel not in channels:
                                channels[channel] = []
                            channels[channel].append(dir)
                        elif self.verbose:
                            log.warning(
                                "not a valid embedding dataset name: %s"
                                % basename)
                    elif self.verbose:
                        log.warning("skipping file: %s" % dir)

                for channel, channel_dirs in channels.items():
                    syst = {}
                    for dir in channel_dirs:
                        dirname, basename = os.path.split(dir)
                        match = re.match(EMBED_PATTERN, basename)
                        if match:
                            isol = match.group('sys')
                            if isol not in syst:
                                syst[isol] = []
                            syst[isol].append(dir)
                        elif self.verbose:
                            log.warning(
                                "not a valid embedding dataset name: %s"
                                % basename)

                    for syst_type, dirs in syst.items():
                        name = 'embed%d-%s-%s' % (
                                year % 1000, channel, syst_type)
                        self[name] = Dataset(
                            name,
                            datatype=EMBED,
                            treename=embed_treename,
                            ds=name,
                            id=1,
                            grl=data_grl,
                            dirs=dirs,
                            file_pattern=embed_pattern,
                            year=year)

            elif embed_sampletype == 'standard':
                if year == 2011:
                    EMBED_PATTERN = EMBED_PATTERN11
                else:
                    EMBED_PATTERN = EMBED_PATTERN12

                # determine what channels are available
                channels = {}
                for dir in embed_dirs:
                    if os.path.isdir(dir):
                        dirname, basename = os.path.split(dir)
                        match = re.match(EMBED_PATTERN, basename)
                        if match:
                            channel = match.group('channel')
                            if channel not in channels:
                                channels[channel] = []
                            channels[channel].append(dir)
                        elif self.verbose:
                            log.warning(
                                "not a valid embedding dataset name: %s"
                                % basename)
                    elif self.verbose:
                        log.warning("skipping file: %s" % dir)

                for channel, channel_dirs in channels.items():
                    if year == 2011:
                        # group dirs by isolation
                        isols = {}
                        for dir in channel_dirs:
                            dirname, basename = os.path.split(dir)
                            match = re.match(EMBED_PATTERN, basename)
                            if match:
                                isol = match.group('isol')
                                if isol not in isols:
                                    isols[isol] = []
                                isols[isol].append(dir)
                            elif self.verbose:
                                log.warning(
                                    "not a valid embedding dataset name: %s"
                                    % basename)

                        for isol, isol_dirs in isols.items():
                            # group dirs by mfs
                            mfss = {}
                            for dir in isol_dirs:
                                dirname, basename = os.path.split(dir)
                                match = re.match(EMBED_PATTERN, basename)
                                if match:
                                    mfs = match.group('mfs')
                                    if mfs not in mfss:
                                        mfss[mfs] = []
                                    mfss[mfs].append(dir)
                                elif self.verbose:
                                    log.warning(
                                        "not a valid embedding dataset name: %s"
                                        % basename)

                            for mfs, mfs_dirs in mfss.items():
                                name = 'embed%d-%s-%s-%s' % (
                                    year % 1000, channel, isol, mfs)

                                self[name] = Dataset(
                                    name,
                                    datatype=EMBED,
                                    treename=embed_treename,
                                    ds=name,
                                    id=1,
                                    grl=data_grl,
                                    dirs=mfs_dirs,
                                    file_pattern=embed_pattern,
                                    year=year)

                                periods = {}
                                for dir in mfs_dirs:
                                    dirname, basename = os.path.split(dir)
                                    match = re.match(EMBED_PATTERN, basename)
                                    if match:
                                        period = match.group('period')
                                        tag = match.group('tag')
                                        if period not in periods:
                                            periods[period] = {'tag': tag, 'dirs': [dir]}
                                        else:
                                            periods[period]['dirs'].append(dir)
                                            if tag != periods[period]['tag']:
                                                log.warning(
                                                    'multiple copies of run with '
                                                    'different tags: %s' %
                                                    periods[period]['dirs'])
                                    elif self.verbose:
                                        log.warning(
                                            "not a valid embeding dataset name: %s"
                                            % basename)

                                for period, info in periods.items():
                                    period_name = '%s-%s' % (name, period)

                                    self[period_name] = Dataset(
                                        name=period_name,
                                        datatype=EMBED,
                                        treename=embed_treename,
                                        ds=period_name,
                                        id=1,
                                        grl=data_grl,
                                        dirs=info['dirs'],
                                        file_pattern=embed_pattern,
                                        year=year)

                    else:
                        # group dirs by mfs
                        mfss = {}
                        for dir in channel_dirs:
                            dirname, basename = os.path.split(dir)
                            match = re.match(EMBED_PATTERN, basename)
                            if match:
                                mfs = match.group('mfs')
                                if mfs not in mfss:
                                    mfss[mfs] = []
                                mfss[mfs].append(dir)
                            elif self.verbose:
                                log.warning(
                                    "not a valid embedding dataset name: %s"
                                    % basename)

                        for mfs, mfs_dirs in mfss.items():
                            name = 'embed%d-%s-%s' % (
                                year % 1000, channel, mfs)
                            self[name] = Dataset(
                                name,
                                datatype=EMBED,
                                treename=embed_treename,
                                ds=name,
                                id=1,
                                grl=data_grl,
                                dirs=mfs_dirs,
                                file_pattern=embed_pattern,
                                year=year)

                            periods = {}
                            for dir in mfs_dirs:
                                dirname, basename = os.path.split(dir)
                                match = re.match(EMBED_PATTERN, basename)
                                if match:
                                    period = match.group('period')
                                    tag = match.group('tag')
                                    if period not in periods:
                                        periods[period] = {'tag': tag, 'dirs': [dir]}
                                    else:
                                        periods[period]['dirs'].append(dir)
                                        if tag != periods[period]['tag']:
                                            log.warning(
                                                'multiple copies of run with '
                                                'different tags: %s' %
                                                periods[period]['dirs'])
                                elif self.verbose:
                                    log.warning(
                                        "not a valid embedding dataset name: %s"
                                        % basename)

                            for period, info in periods.items():
                                period_name = '%s-%s' % (name, period)
                                self[period_name] = Dataset(
                                    name=period_name,
                                    datatype=EMBED,
                                    treename=embed_treename,
                                    ds=period_name,
                                    id=1,
                                    grl=data_grl,
                                    dirs=info['dirs'],
                                    file_pattern=embed_pattern,
                                    year=year)

            elif embed_sampletype == 'lhCN':
                year = 2012
                channels = {}
                for dir in embed_dirs:
                    if os.path.isdir(dir):
                        dirname, basename = os.path.split(dir)
                        match = re.match(CN_EMBED_PATTERN12, basename)
                        if match:
                            channel = match.group('channel')
                            if channel not in channels:
                                channels[channel] = []
                            channels[channel].append(dir)
                        elif self.verbose:
                            log.warning(
                                "not a valid embedding dataset name: %s"
                                % basename)
                    elif self.verbose:
                        log.warning("skipping file: %s" % dir)

                for channel, channel_dirs in channels.items():
                    # group dirs by mfs
                    mfss = {}
                    for dir in channel_dirs:
                        dirname, basename = os.path.split(dir)
                        match = re.match(CN_EMBED_PATTERN12, basename)
                        if match:
                            mfs = match.group('mfs')
                            if mfs not in mfss:
                                mfss[mfs] = []
                            mfss[mfs].append(dir)
                        elif self.verbose:
                            log.warning(
                                "not a valid embedding dataset name: %s"
                                % basename)

                    for mfs, mfs_dirs in mfss.items():
                        name = 'embed%d-%s-%s' % (
                            year % 1000, channel, mfs)
                        self[name] = Dataset(
                            name,
                            datatype=EMBED,
                            treename=embed_treename,
                            ds=name,
                            id=1,
                            grl=data_grl,
                            dirs=mfs_dirs,
                            file_pattern=embed_pattern,
                            year=year)

            # MC EMBEDDING
            variations = {}
            for dir in embed_dirs:
                dirname, basename = os.path.split(dir)
                match = re.match(MC_EMBED_PATTERN, basename)
                if not match:
                    continue
                syst = match.group('sys') or ''
                variations.setdefault(syst, []).append(dir)
            for variation, dirs in variations.items():
                name = 'mcembed12-hh%s' % variation
                self[name] = Dataset(
                    name,
                    datatype=MCEMBED,
                    treename=embed_treename,
                    ds=name,
                    id=1,
                    dirs=dirs,
                    file_pattern=embed_pattern,
                    year=2012)

        ##############################
        # DATA
        ##############################
        if data_path is not None:

            if deep:
                data_dirs = get_all_dirs_under(data_path, prefix=data_prefix)
            else:
                if data_prefix:
                    data_dirs = glob.glob(
                        os.path.join(data_path, data_prefix) + '*')
                else:
                    data_dirs = glob.glob(
                        os.path.join(data_path, '*'))

            if data_sampletype == 'standard':
                # classify dir by stream
                streams = {}
                for dir in data_dirs:
                    dirname, basename = os.path.split(dir)
                    match = re.match(DS_PATTERN, basename)
                    if match:

                        # pass embed
                        if re.match(EMBED_PATTERN12_NEW, basename) or \
                                re.match(EMBED_PATTERN12, basename) or \
                                re.match(EMBED_PATTERN11, basename) :
                            continue

                        if int(match.group('year')) != (year % 1E3):
                            continue
                        if match.group('type') != 'data':
                            continue
                        stream = match.group('name').split('_')[-1]
                        if stream not in streams:
                            streams[stream] = []
                        streams[stream].append(dir)
                    elif self.verbose:
                        log.warning(
                            "not a valid data dataset name: %s" % basename)

                for stream, dirs in streams.items():
                    name = 'data%d-%s' % (year % 1000, stream)
                    self[name] = Dataset(
                        name=name,
                        datatype=DATA,
                        treename=data_treename,
                        ds=name,
                        id=-1,
                        # The GRL is the same for both lephad and hadhad analyses
                        grl=data_grl,
                        dirs=dirs,
                        stream=stream,
                        file_pattern=data_pattern,
                        year=year)

                    if data_period_containers:
                        # in each stream create a separate dataset for each run
                        periods = {}
                        for dir in dirs:
                            dirname, basename = os.path.split(dir)
                            match = re.match(DS_PATTERN, basename)
                            if match:
                                period = match.group('id')
                                if not period.startswith('period'):
                                    continue
                                tag = match.group('tag')
                                if period not in periods:
                                    periods[period] = {
                                        'tag': tag,
                                        'dirs': [dir],
                                        'ds': Database.match_to_ds(match)}
                                else:
                                    periods[period]['dirs'].append(dir)
                                    if tag != periods[period]['tag']:
                                        log.warning(
                                            'multiple copies of period with different '
                                            'tags: \n%s' %
                                            ('\n'.join(periods[period]['dirs'])))
                            elif self.verbose:
                                log.warning(
                                    "not a valid data dataset name: %s" % basename)
                        # need to use the actual ds name for ds for validation
                        for period, info in periods.items():
                            name = 'data%d-%s-%s' % (year % 1000, stream, period[-1])
                            self[name] = Dataset(
                                name=name,
                                datatype=DATA,
                                treename=data_treename,
                                ds=name,
                                id=period,
                                grl=data_grl,
                                dirs=info['dirs'],
                                stream=stream,
                                file_pattern=data_pattern,
                                year=year)
                    else:
                        # in each stream create a separate dataset for each run
                        runs = {}
                        for dir in dirs:
                            dirname, basename = os.path.split(dir)
                            match = re.match(DS_PATTERN, basename)
                            if match:
                                run = int(match.group('id'))
                                tag = match.group('tag')
                                if run not in runs:
                                    runs[run] = {
                                        'tag': tag,
                                        'dirs': [dir],
                                        'ds': Database.match_to_ds(match)}
                                else:
                                    runs[run]['dirs'].append(dir)
                                    if tag != runs[run]['tag']:
                                        log.warning(
                                            'multiple copies of run with different '
                                            'tags: %s' % runs[run]['dirs'])
                            elif self.verbose:
                                log.warning(
                                    "not a valid data dataset name: %s" % basename)
                        # need to use the actual ds name for ds for validation
                        for run, info in runs.items():
                            name = 'data%d-%s-%d' % (year % 1000, stream, run)
                            self[name] = Dataset(
                                name=name,
                                datatype=DATA,
                                treename=data_treename,
                                ds=name,
                                id=run,
                                grl=data_grl,
                                dirs=info['dirs'],
                                stream=stream,
                                file_pattern=data_pattern,
                                year=year)

                        if USE_PYAMI:
                            # in each stream create a separate dataset for each period
                            run_periods = get_periods(amiclient, year=year, level=2)
                            # ignore subset periods like Ba in 2012
                            run_periods = [
                                p.name for p in run_periods if len(p.name) == 1]
                            period_runs = {}
                            for period in run_periods:
                                if period == 'VdM':
                                    continue
                                _runs = get_runs(amiclient, periods=period, year=year)
                                for run in _runs:
                                    period_runs[run] = period
                            periods = {}
                            for run, info in runs.items():
                                if run in period_runs:
                                    _period = period_runs[run]
                                else:
                                    # ignore spurious runs
                                    continue
                                if _period in periods:
                                    periods[_period] += info['dirs']
                                else:
                                    periods[_period] = info['dirs'][:]
                            for period, dirs in periods.items():
                                name = 'data%d-%s-%s' % (year % 1000, stream, period)
                                self[name] = Dataset(
                                    name=name,
                                    datatype=DATA,
                                    treename=data_treename,
                                    ds=name,
                                    id=-1,
                                    grl=data_grl,
                                    dirs=dirs,
                                    stream=stream,
                                    file_pattern=data_pattern,
                                    year=year)

            elif data_sampletype == 'lhCN':
                year = 2012
                streams = {}
                for dir in data_dirs:
                    match = re.match(CN_DATA_PATTERN12, dir)
                    if match:
                        stream = match.group('name')
                        if stream not in streams:
                            streams[stream] = []
                        streams[stream].append(dir)
                    elif self.verbose:
                        log.warning("not a valid data dataset name: %s" % dir)

                for stream, dirs in streams.items():
                    name = 'data%d-%s' % (year % 1000, stream)
                    log.info('\'%s\',' % name)
                    self[name] = Dataset(
                        name=name,
                        datatype=DATA,
                        treename=data_treename,
                        ds=name,
                        id=-1,
                        # The GRL is the same for both lephad and hadhad analyses
                        grl=data_grl,
                        dirs=dirs,
                        stream=stream,
                        file_pattern=data_pattern,
                        year=year)

                    # in each stream create a separate dataset for each period
                    periods = {}
                    for dir in dirs:
                        match = re.match(CN_DATA_PATTERN12, dir)
                        if match:
                            period = match.group('id')
                            tag = match.group('tag')
                            if period not in periods:
                                periods[period] = {
                                    'tag': tag,
                                    'dirs': [dir],
                                    'ds': -1}
                            else:
                                periods[period]['dirs'].append(dir)
                                if tag != periods[period]['tag']:
                                    log.warning(
                                        'multiple copies of period with different '
                                        'tags: %s' % periods[period]['dirs'])
                        elif self.verbose:
                            log.warning(
                                "not a valid data dataset name: %s" % dir)
                    # need to use the actual ds name for ds for validation
                    for period, info in periods.items():
                        name = 'data%d-%s-%s' % (year % 1000, stream, period)
                        log.info('\'%s\',' % name)
                        self[name] = Dataset(
                            name=name,
                            datatype=DATA,
                            treename=data_treename,
                            ds=info['ds'],
                            id=period,
                            grl=data_grl,
                            dirs=info['dirs'],
                            stream=stream,
                            file_pattern=data_pattern,
                            year=year)

    def __setitem__(self, name, ds):
        if self.verbose:
            print >> self.stream, str(ds)
        super(Database, self).__setitem__(name, ds)

    def search(self, pattern):
        data = []
        patterns = pattern
        if not isinstance(pattern, (list, tuple)):
            patterns = [pattern]
        for name, ds in self.items():
            for pattern in patterns:
                if fnmatch.fnmatch(name, pattern):
                    data.append(ds)
                    continue
                if not pattern.startswith('^'):
                    pattern = '^' + pattern
                if not pattern.endswith('$'):
                    pattern = pattern + '$'
                if re.match(pattern, name):
                    data.append(ds)
                    continue
        return data


class Dataset(Serializable):

    yaml_tag = u'!Dataset'

    def __init__(self, name, datatype, treename, ds, dirs,
                 file_pattern='*.root*',
                 id=None,
                 category=None,
                 version=None,
                 tag_pattern=None,
                 tag=None,
                 grl=None,
                 year=None,
                 stream=None):
        self.name = name
        self.datatype = datatype
        self.treename = treename
        self.id = id
        self.ds = ds
        self.category = category
        self.version = version
        self.tag_pattern = tag_pattern
        self.tag = tag
        self.dirs = dirs
        self.file_pattern = file_pattern
        self.grl = grl
        self.year = year
        self.stream = stream

    def __repr__(self):
        return ("%s(name=%r, datatype=%r, treename=%r, "
                "id=%r, ds=%r, category=%r, version=%r, "
                "tag_pattern=%r, tag=%r, dirs=%r, "
                "file_pattern=%r, grl=%r, year=%r, stream=%r)") % (
                        self.__class__.__name__,
                        self.name, self.datatype, self.treename,
                        self.id, self.ds, self.category, self.version,
                        self.tag_pattern, self.tag, self.dirs,
                        self.file_pattern, self.grl, self.year, self.stream)

    @cached_property
    def xsec_kfact_effic(self):
        global XSEC_CACHE_MODIFIED
        year = self.year % 1E3
        if self.datatype == DATA:
            return 1., 1., 1.
        if year in XSEC_CACHE and self.name in XSEC_CACHE[year]:
            log.warning("using cached cross section for dataset %s" % self.ds)
            return XSEC_CACHE[year][self.name]
        try:
            return xsec.xsec_kfact_effic(self.year, self.id)
        except KeyError:
            log.warning("cross section of dataset %s not available locally."
                        "Looking it up in AMI instead. AMI cross sections can be very"
                        "wrong! You have been warned!"
                        % self.ds)
        if USE_PYAMI:
            if self.ds in DS_NOPROV:
                xs, effic = get_dataset_xsec_effic(amiclient, DS_NOPROV[self.ds])
            else:
                xs, effic = get_dataset_xsec_effic(amiclient, self.ds)
            if year not in XSEC_CACHE:
                XSEC_CACHE[year] = {}
            XSEC_CACHE[year][self.name] = (xs, 1., effic)
            XSEC_CACHE_MODIFIED = True
            return xs, 1., effic
        raise Exception("cross section of dataset %s is not known!" % self.ds)

    @cached_property
    def files(self):
        if not self.dirs:
            log.warning(
                "files requested from dataset %s "
                "with an empty list of directories" % self.name)
        _files = []
        for dir in self.dirs:
            if not os.path.exists(dir):
                raise IOError("%s is not readable" % dir)
            for path, dirs, files in os.walk(dir):
                _files += [os.path.join(path, f) for f in
                           fnmatch.filter(files, self.file_pattern)]
        return _files

    def __str__(self):
        return "%s (%d files):\n\t%s" % (
                self.name,
                len(self.files),
                self.ds)

def dataset_constructor(loader, node):
    kwargs = loader.construct_mapping(node)
    try:
        return Dataset(**kwargs)
    except:
        fields = '\n'.join('%s = %s' % item for item in kwargs.items())
        log.error("unable to load dataset %s with these fields:\n\n%s\n" %
                  (kwargs['name'], fields))
        raise

yaml.add_constructor(u'!Dataset', dataset_constructor)

if os.path.isfile(XSEC_CACHE_FILE):
    with open(XSEC_CACHE_FILE) as cache:
        log.info("Loading cross section cache in %s ..." % XSEC_CACHE_FILE)
        XSEC_CACHE = pickle.load(cache)


@atexit.register
def write_cache():
    if XSEC_CACHE_MODIFIED:
        with open(XSEC_CACHE_FILE, 'w') as cache:
            log.info("Saving cross-section cache to disk...")
            pickle.dump(XSEC_CACHE, cache)


def validate_single(args, child=True):
    if child:
        from cStringIO import StringIO
        sys.stdout = out = StringIO()
        sys.stderr = out
    name = args[0]
    info = args[1]
    complete = True
    try:
        dirs = info.dirs
        root_files = []
        for dir in dirs:
            root_files += glob.glob(os.path.join(dir, info.file_pattern))
        events = 0
        for fname in root_files:
            try:
                with root_open(fname) as rfile:
                    try: # skimmed dataset
                        events += int(rfile.cutflow_event[0])
                    except DoesNotExist: # unskimmed dataset
                        tree = rfile.tau
                        events += tree.GetEntries()
            except IOError:
                log.warning("Currupt file: %s" % fname)
                pass
        # determine events in original ntuples
        # use first dir
        ds_name = info.ds
        log.info('NTUP: ' + ds_name)
        ds_info = get_dataset_info(amiclient, ds_name)
        ntuple_events = int(ds_info.info['totalEvents'])
        try:
            # determine events in AODs
            prov = get_provenance(amiclient, ds_name, type='AOD')
            AOD_ds = prov.values()[0][0].replace('recon', 'merge')
            log.info('AOD: ' + AOD_ds)
            AOD_events = int(get_datasets(amiclient, AOD_ds, fields='events',
                    flatten=True)[0][0])
        except IndexError:
            log.info('AOD: UNKNOWN')
            AOD_events = ntuple_events
        log.info(name)
        log.info("\tevts\tNTUP\tAOD")
        log.info("\t%i\t%i\t%i" % (events, ntuple_events, AOD_events))
        if events != ntuple_events:
            log.warning("NTUP MISMATCH")
        if events != AOD_events:
            log.warning("AOD MISMATCH")
        if events != ntuple_events and (events != AOD_events or AOD_events == 0):
            log.warning("MISSING EVENTS")
            complete = False
        if child:
            return out.getvalue(), complete
        return complete
    except Exception, e:
        import traceback
        log.warning("dataset %s exception" % name)
        traceback.print_exception(*sys.exc_info())
        if child:
            return out.getvalue(), False
        return False


def get_all_dirs_under(path, prefix=None):
    """
    Get list of all directories under path
    """
    dirs = []
    for dirpath, dirnames, filenames in os.walk(path):
        _dirnames = []
        for dirname in dirnames:
            fullpath = os.path.join(dirpath, dirname)
            # check if this dir contains other dirs
            subdirs_exist = False
            subdirs = os.listdir(fullpath)
            for subdir in subdirs:
                if os.path.isdir(os.path.join(fullpath, subdir)):
                    subdirs_exist = True
                    break
            if subdirs_exist:
                _dirnames.append(dirname)
            else:
                # this must be a dataset, don't walk into this dir
                if prefix is not None:
                    if not dirname.startswith(prefix):
                        continue
                dirs.append(fullpath)
        # only recurse on directories containing subdirectories
        dirnames = _dirnames
    return dirs
