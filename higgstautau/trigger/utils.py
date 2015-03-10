import os
import ROOT

from externaltools import TrigRootAnalysis
from ROOT import D3PD


__all__ = [
    'get_trigger_config',
    'update_trigger_config',
    'get_tau_trigger_obj_idx',
]


def get_trigger_config():

    return D3PD.TrigConfigSvcD3PD()


def update_trigger_config(tool, name, file, tree):

    configTree = file.Get("%sMeta/TrigConfTree" % name)
    tool.SetConfigTree(configTree)


def get_tau_trigger_obj_idx(config, event, trigger):
    """
    Return list of trigger tau indicies associated with
    RoIs for this trigger
    """
    chainID = config.GetChainId(trigger)

    trigger_idx = []
    for i, id in enumerate(event.trig_Nav_chain_ChainId):
        if id == chainID + 10000:
            for roi_idx in event.trig_Nav_chain_RoIIndex[i]:
                if len(event.trig_RoI_EF_tau_Analysis__TauJetContainer) == 0:
                    continue
                try:
                    for k, status in enumerate(event.trig_RoI_EF_tau_Analysis__TauJetContainerStatus[roi_idx]):
                        if status != 1:
                            continue
                        idx = event.trig_RoI_EF_tau_Analysis__TauJetContainer[roi_idx][k]
                        trigger_idx.append(idx)
                except IndexError:
                    # saw this happen once...
                    # ...Analysis__TauJetContainerStatus[roi_idx]
                    pass
    # remove duplicates
    return list(set(trigger_idx))
