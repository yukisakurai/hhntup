#!/usr/bin/env python

from rootpy.io import root_open as ropen
from rootpy.tree import Cut
from rootpy.plotting import Graph
from ..common import PRONGS, LEVELS, nprong
import os


HERE = os.path.dirname(os.path.abspath(__file__))

CATEGORIES = {
    '_3': Cut('tau_numberOfVertices<=3'),
    '3_5': Cut('3<tau_numberOfVertices<=5'),
    '5_7': Cut('5<tau_numberOfVertices<=7'),
    '7_': Cut('tau_numberOfVertices>7'),
}


if __name__ == '__main__':

    with ropen(os.path.join(HERE, 'bdt_selection.root'), 'recreate') as f:

        for prong in PRONGS:
            for cat_str, category in CATEGORIES.items():
                for level_name, level in LEVELS.items():
                    fname = 'sig-bits-%dp-%s-perfB--%d.txt' % (
                            prong, category.safe(parentheses=False), level)
                    with open(fname) as fin:
                        lines = fin.readlines()[1:]
                        graph = Graph(len(lines), name='%s_%dp_%s' % (
                            level_name, prong, cat_str))
                        for i, line in enumerate(lines):
                            pt, bdt = map(float, line.strip().split())
                            graph[i] = (pt, bdt)
                        graph.Write()
else:

    P851_SELECTION = {}
    with ropen(os.path.join(HERE, 'bdt_selection.root')) as f:
        for level in LEVELS.keys():
            P851_SELECTION[level] = {}
            for prong in PRONGS:
                P851_SELECTION[level][prong] = {}
                for category in CATEGORIES.keys():
                    P851_SELECTION[level][prong][category] = f.Get(
                            '%s_%dp_%s' % (level, prong, category)).Clone()


    def nvtx_to_category(nvtx):

        if isinstance(nvtx, str):
            return nvtx
        if nvtx <= 3:
            category = '_3'
        elif nvtx <= 5:
            category = '3_5'
        elif nvtx <= 7:
            category = '5_7'
        else:
            category = '7_'
        return category


    def selection(level, prong, nvtx):

        prong = nprong(prong)
        return P851_SELECTION[level][prong][nvtx_to_category(nvtx)]
