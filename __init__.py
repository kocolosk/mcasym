# -*- coding: utf-8 -*-
from _aac06 import aacpdfe
from _dns import polfit, dnsini
from _dssv import dssvfit, dssvini
from _polnlo import polnlo, nloini
import _grsv
import _grsv2000
import _lss2006
import _bb

from _ctq5par import ctq5pd

## initialization routines that only need to be called once
dnsini()
dssvini()
nloini()

_flavor_lut = {1:'d', 2:'u', -1:'dbar', -2:'ubar', 3:'s', -3:'sbar', 4:'c', -4:'cbar',
    5:'b', -5:'bbar', 6:'t', -6:'tbar', 21:'gluon'}

_grsv_iset = {'NLO':1, 'STD':1, 'MIN':2, 'ZERO':3, 'MAX':4, 'M015':5, 'M030':6,
    'M045':7, 'M060':8, 'M075':9, 'M090':10, 'M105':11, 'P030':12, 'P045':13,
    'P060':14, 'P070':15}

_polnlo_iset = {'GS_NLOA':0, 'GS_NLOB':1, 'GS_NLOC':2}

_lss2006_iset = {'LSS1':1, 'LSS2':2, 'LSS3':3}

_aac06_iset = {'AAC1':1, 'AAC2':2, 'AAC3':3}

_bb_iset = {'BB1':3, 'BB2':4}

_dns_iset = {'DNS1':1, 'DNS2':2}

## Pythia => u_valence + u_sea
## so try to match it in polpdf

def num(iset, event):
    flavors = [event.flavor(i+1) for i in range(4)]
    df1 = polpdf(flavors[0], iset, event.x1(), event.Q2())
    df2 = polpdf(flavors[1], iset, event.x2(), event.Q2())
    asym = partonicAsymmetry(flavors, event.processId(), event.s(), event.t(), event.u())
    return df1*df2*asym


def denom(iset, event):
    f1 = pdf(event.flavor(1), iset, event.x1(), event.Q2())
    f2 = pdf(event.flavor(2), iset, event.x2(), event.Q2())
    return f1*f2


def polpdf(flavor, iset, x, Q2):
    '''Δf(x) according to one of several PDF parameterizations
    
    flavor: PDG code for a quark or gluon
    iset: GRSV, GS, LSS, AAC, BB, DNS, and DSSV sets are available
    
    >>> import mcasym
    
    These are just regression tests, they really shouldn't be here in the
    docstring:
    
    >>> print round( mcasym.polpdf(2, 'LO'     , 0.5223, 302.81), 6 )
    0.198335
    
    >>> print round( mcasym.polpdf(2, 'NLO'    , 0.5223, 302.81), 6 )
    0.179623
    
    >>> print round( mcasym.polpdf(2, 'MAX'    , 0.5223, 302.81), 6 )
    0.178866
    
    >>> print round( mcasym.polpdf(2, 'MIN'    , 0.5223, 302.81), 6 )
    0.173255
    
    >>> print round( mcasym.polpdf(2, 'ZERO'   , 0.5223, 302.81), 6 )
    0.177149
    
    >>> print round( mcasym.polpdf(2, 'M015'   , 0.5223, 302.81), 6 )
    0.178162
    
    >>> print round( mcasym.polpdf(2, 'M030'   , 0.5223, 302.81), 6 )
    0.177678
    
    >>> print round( mcasym.polpdf(2, 'M045'   , 0.5223, 302.81), 6 )
    0.17579
    
    >>> print round( mcasym.polpdf(2, 'M060'   , 0.5223, 302.81), 6 )
    0.179609
    
    >>> print round( mcasym.polpdf(2, 'M075'   , 0.5223, 302.81), 6 )
    0.179027
    
    >>> print round( mcasym.polpdf(2, 'M090'   , 0.5223, 302.81), 6 )
    0.178198
    
    >>> print round( mcasym.polpdf(2, 'M105'   , 0.5223, 302.81), 6 )
    0.176486
    
    >>> print round( mcasym.polpdf(2, 'P030'   , 0.5223, 302.81), 6 )
    0.179384
    
    >>> print round( mcasym.polpdf(2, 'P045'   , 0.5223, 302.81), 6 )
    0.179025
    
    >>> print round( mcasym.polpdf(2, 'P060'   , 0.5223, 302.81), 6 )
    0.178174
    
    >>> print round( mcasym.polpdf(2, 'P070'   , 0.5223, 302.81), 6 )
    0.178352
    
    >>> print round( mcasym.polpdf(2, 'GS_NLOA', 0.5223, 302.81), 6 )
    0.185317
    
    >>> print round( mcasym.polpdf(2, 'GS_NLOB', 0.5223, 302.81), 6 )
    0.184054
    
    >>> print round( mcasym.polpdf(2, 'GS_NLOC', 0.5223, 302.81), 6 )
    0.18058
    
    >>> print round( mcasym.polpdf(2, 'LSS1'   , 0.5223, 302.81), 6 )
    0.177305
    
    >>> print round( mcasym.polpdf(2, 'LSS2'   , 0.5223, 302.81), 6 )
    0.179846
    
    >>> print round( mcasym.polpdf(2, 'LSS3'   , 0.5223, 302.81), 6 )
    0.178731
    
    >>> print round( mcasym.polpdf(2, 'AAC1'   , 0.5223, 302.81), 6 )
    0.181222
    
    >>> print round( mcasym.polpdf(2, 'AAC2'   , 0.5223, 302.81), 6 )
    0.180826
    
    >>> print round( mcasym.polpdf(2, 'AAC3'   , 0.5223, 302.81), 6 )
    0.181822
    
    >>> print round( mcasym.polpdf(2, 'BB1'    , 0.5223, 302.81), 6 )
    0.208463
    
    >>> print round( mcasym.polpdf(2, 'BB2'    , 0.5223, 302.81), 6 )
    0.208008
    
    >>> print round( mcasym.polpdf(2, 'DNS1'   , 0.5223, 302.81), 6 )
    0.170999
    
    >>> print round( mcasym.polpdf(2, 'DNS2'   , 0.5223, 302.81), 6 )
    0.172109
    
    >>> print round( mcasym.polpdf(2, 'DSSV'   , 0.5223, 302.81), 6 )
    0.183847
    '''
    if iset in _grsv_iset.keys():
        _grsv.intini.iini = 0
        u,d,ubar,dbar,s,gluon,g1p,g1n = _grsv.parpol2(_grsv_iset[iset], x, Q2)
        sbar = s
    elif iset == 'LO':
        _grsv2000.intini.iini = 0
        u,d,ubar,dbar,s,gluon,g1p,g1n = _grsv2000.parpol(3, x, Q2)
        sbar = s
    elif iset == 'DSSV':
        uv,dv,ubar,dbar,s,gluon = dssvfit(x, Q2)
        u = uv + ubar
        d = dv + dbar
        sbar = s
    elif iset in _polnlo_iset.keys():
        uv,dv,gluon,ubar,dbar,s = polnlo(_polnlo_iset[iset],x,Q2)
        u = uv + ubar
        d = dv + dbar
        sbar = s
    elif iset in _lss2006_iset.keys():
        _lss2006.intini.iini = 0
        uub,ddb,ssb,gluon,uv,dv,ubar,dbar,sbar,g1plt,g1p,g1nlt,g1n = \
            _lss2006.lss2006(_lss2006_iset[iset],x,Q2)
        u = uub - ubar
        d = ddb - dbar
        s = ssb - sbar
    elif iset in _aac06_iset.keys():
        sbar,dbar,ubar,gluon,u,d,s = aacpdfe(Q2,x,_aac06_iset[iset])
    elif iset in _bb_iset.keys():
        _bb.intini.iini = 0
        uv,duv,dv,ddv,gluon,dgl,qbar,dqb,g1p,dg1p,g1n,dg1n = \
            _bb.ppdf(_bb_iset[iset],x,Q2)
        u = uv + qbar
        d = dv + qbar
        ubar = dbar = sbar = s = qbar
    elif iset in _dns_iset.keys():
        uv,dv,ubar,dbar,s,gluon,g1p,g1n = polfit(_dns_iset[iset],x,Q2)
        u = uv + ubar
        d = dv + dbar
        sbar = s
    else:
        raise ValueError, '%s is not a supported PDF set' % iset
    
    t = tbar = b = bbar = c = cbar = sbar
    return locals()[_flavor_lut[flavor]]/x


def pdf(flavor, iset, x, Q2):
    '''Unpolarized PDFs from CTEQ5L (LO) and CTEQ5M1 (NLO).
    
    flavor: PDG code for a quark or gluon
    iset: 'LO' or 'NLO'
    
    >>> import mcasym
    
    >>> print round(mcasym.pdf(21, 'LO', 0.001, 30.0), 3)
    20567.16

    >>> print round(mcasym.pdf(21, 'NLO', 0.001, 30.0), 3)
    14004.477
    '''
    from math import sqrt
    cteq_iset = {'LO':3, 'NLO':1}
    if flavor == 21: flavor = 0
    answer, error = ctq5pd(cteq_iset[iset], flavor, x, sqrt(Q2))
    if error != 0: 
        raise ValueError, "invalid x, Q2 for CTEQ parameterization: %f,%f" % (x,Q2)
    return answer


def _weights(s,t,u,index):
    if index == 1:
        w = [ (8*s*s)/(9*t*t),
        (8*u*u)/(9*t*t) ]
    elif index == 2:
        w = [ 8*(s*s/(t*t)+s*s/(u*u)-2.0/3*s*s/(t*u))/9,
        8*(u*u/(t*t)+t*t/(u*u))/9 ]
    elif index == 3:
        w = [ 0,
        8*(t*t+u*u)/(s*s*9) ]
    elif index == 4:
        w = [ 8*s*s/(9*t*t),
        8*(u*u/(t*t)+(t*t+u*u)/(s*s)-2*u*u/(s*t*3))/9 ]
    elif index == 5:
        w = [ 0,
        64*(t*t+u*u)/(27*u*t)-16*(t*t+u*u)/(3*s*s) ]
    elif index == 6:
        w = [ 0,
        (u*u+t*t)/(3*u*t)-3*(t*t+u*u)/(4*s*s) ]
    elif index == 7:
        w = [ 2*s*s/(t*t)-8*s*s/(9*u*s),
        2*u*u/(t*t)-8*u*u/(9*u*s) ]
    elif index == 8:
        w = [ 4.5*(2*s*s/(u*t)-s*u/(t*t)-s*t/(u*u)),
        4.5*(6-2*s*s/(u*t)-s*u/(t*t)-s*t/(u*u)-2*u*t/(s*s)) ]
    else:
        raise KeyError, "unsupported index for Werner's weights: %d" % index 
    
    return w


def partonicAsymmetry(flavors, processId, s, t, u):
    '''pQCD asymmetry for various QCD 2->2 processes.
    
    flavors: PDG codes for partons 1,2,3,4
    processId: one of the QCD 2->2 processes 11,12,13,28,53,68
    s,t,u: Mandelstam variables
    
    >>> from mcasym import partonicAsymmetry as pasym
    
    Here are some tests for qiqj -> qiqj scattering.  Different matrix elements are
    required depending on whether i=j, i=-j, or other
    
    >>> print round(pasym(( 2, 2, 2, 2), 11,  3226.3932, -2887.6768, -0338.7164),6)
    0.077358
    >>> print round(pasym(( 2, 1, 2, 1), 11,  2454.9072, -0275.5169, -2179.3903),6)
    0.118485
    >>> print round(pasym((-2, 2,-2, 2), 11,  1836.2634, -0566.4859, -1269.7775),6)
    0.226179
    
    That last one was wrong in StMCAsymMaker for a long time.  Next we have 
    qqbar -> q'qbar'.  The first one is the same result for pid 11:
    
    >>> print round(pasym((-1, 1,-1, 1), 12,  1836.2634, -0566.4859, -1269.7775),6)
    0.226179
    >>> print round(pasym((-2, 2,-3, 3), 12,  0945.4256, -0569.3693, -0376.0563),6)
    -1.0
    
    qqbar -> gg
    
    >>> print round(pasym((-2, 2,21,21), 13,  1016.7810, -0644.5507, -0372.2303),6)
    -1.0
    
    qg -> qg
    
    >>> print round(pasym(( 1,21, 1,21), 28,  1669.9881, -0493.8904, -1176.0977),6)
    0.33692
    
    gg -> qqbar
    
    >>> print round(pasym((21,21,-4, 4), 53,  1281.6422, -0432.5604, -0849.0818),6)
    -1.0
    
    gg -> gg
    
    >>> print round(pasym((21,21,21,21), 68,  1078.3770, -0409.0749, -0669.3021),6)
    0.71072
    '''
    index = {
        28 : 7, ## qg->qg
        53 : 6, ## gg->qqbar
        68 : 8, ## gg->gg
        12 : 3, ## qqbar->q'qbar'
        13 : 5, ## qqbar->gg
        11 : 1  ## qq'->qq'
    }[processId]
    
    if processId == 12 and abs(flavors[0]) == abs(flavors[2]):
        index = 4
    if processId == 11 and abs(flavors[0]) == abs(flavors[1]):
        index = (flavors[0]==flavors[1]) and 2 or 4
        #index = 2  ## was a bug in StMCAsymMaker
    
    t = _weights(s,t,u,index)
    return (t[0]-t[1])/(t[0]+t[1])

