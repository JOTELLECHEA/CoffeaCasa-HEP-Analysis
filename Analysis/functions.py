import matplotlib.pyplot as plt
import mplhep as hep
import awkward as ak
import numpy as np
import uproot
import hist

def ChargeOrdering(LEPTON):
    ''' Sort charge from each event negative first then postive. '''
    index = ak.argsort(LEPTON.charge)
    return LEPTON[index]

def LeptonSelectionForZBoson(DATA, PT, ETA, FLAVOR):
    ''' Event selection for Z Boson canidates: lepton antilepton pair, lepton.pt > PT, lepton.eta > ETA,
    ensures charge conservation and allows for 2 or 4 lepton pairs. '''
    if FLAVOR == 'Electron': x = DATA.Electron
    elif FLAVOR == 'Muon': x = DATA.Muon
    mask_pT_eta = ak.all(x.pt > PT, axis=-1) & ak.all(x.eta > ETA, axis=-1)
    x = x[mask_pT_eta]
    mask_charge = ak.sum(x.charge, axis=-1) == 0
    mask_2e = ak.num(x, axis=-1) == 2
    mask_4e = ak.num(x, axis=-1) == 4
    mask = mask_charge & (mask_2e | mask_4e)
    return x[mask]

def LeptonSelectionForWBoson(DATA, PT, ETA, FLAVOR, CHARGE):
    ''' Event selection for W Boson canidates: lepton antilepton pair, lepton.pt > PT, lepton.eta > ETA and only 1 lepton.'''
    x = DATA
    mask_pT_eta = ak.all(x[FLAVOR].pt > PT, axis=-1) & ak.all(x[FLAVOR].eta > ETA, axis=-1)
    x = x[mask_pT_eta]
    mask_lepCount = ak.num(x[FLAVOR], axis=-1) == 1
    x = x[mask_lepCount]
    mask_charge = ak.all(x[FLAVOR].charge == CHARGE, axis=-1)
    return x[mask_charge]

def TwoLeptonIM(LEPTON):
    ''' Calculates Invariant Mass: Sum of two lorentz vector and then calculates the new vector's mass. '''
    twoLepVec = LEPTON[:,0].add(LEPTON[:,1])
    return twoLepVec.mass

def FourLeptonIM(LEPTON):
    ''' Calculates Invariant Mass: Sum of two lorentz vector and then calculates the new vector's mass.
     This is done for all possible combination of 4 leptons and a list of the mass is return. ''' 
    MASS = []
    for i in range(2):
        for j in range(2):
            MASS.append((LEPTON[:,i].add(LEPTON[:,j+2])).mass)
    return MASS

def InvariantMassHist(LEPTON, BINS):
    ''' Creates a Histogram that ranges from 0 to 180 GeV. Event selection that makes sure that 2and 4 lepton events are         used.Then fills the histogram. '''
    h = hist.Hist(hist.axis.Regular(BINS, 0, 180, name='Invariant Mass[GeV/C^2]'))
    mask_2l = ak.num(LEPTON, axis=-1) == 2
    mask_4l = ak.num(LEPTON, axis=-1) == 4
    lep2 = LEPTON[mask_2l]
    lep4 = LEPTON[mask_4l]
    h.fill(TwoLeptonIM(lep2))
    for i in range(4):
        h.fill(FourLeptonIM(lep4)[i])
    return h 

def addStats(H):
    ''' This adds a Statistics box for the histogram that is passed. It shows count, mean, and STD. '''
    stats = (np.atleast_1d(H.profile(axis=0).view())[0])
    count , mean , sumDeltaSquared = int(stats['count']), stats['value'], stats['_sum_of_deltas_squared']
    std = np.sqrt(sumDeltaSquared/count)
    statbox = f'{"Entries":<{10}}{count:d}' + '\n' + f'{"Mean":<{10}}{mean:.3f}' + '\n' + f'{"Std Dev":<{10}}{std:.3f}'
    return statbox

def cmsPlot(H):
    ''' Plots histogram and decorates the plot for presentaio. '''
    fig, ax = plt.subplots(figsize=(10,5))
    ax.set_ylabel('Events', fontsize=24)
    ax.set_xlabel('$[GeV/C^{2}$]', fontsize=15)
    hep.histplot(H,label='$M_{ee}$',histtype='fill')
    # hep.cms.label(rlabel="pT > 20 GeV/C")
    plt.tight_layout()
    stats = addStats(H)
    plt.text(0.77, 0.90, stats, ha='left', va='top', transform=ax.transAxes, fontsize = 15, bbox = dict(alpha = 0.15))
    # plt.savefig('InvariantMassPlot.jpg')
    plt.show()