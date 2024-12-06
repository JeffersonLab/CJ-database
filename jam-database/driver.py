#!/usr/bin/env python
import sys,os
import numpy as np
import time

#--matplotlib
import matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['text.latex.preamble']=[r"\usepackage{amsmath}"]
matplotlib.rc('text',usetex=True)
from mpl_toolkits.axes_grid1 import make_axes_locatable
import pylab  as py
from matplotlib.lines       import Line2D
from tools.tools            import checkdir,lprint,load
from obslib.idis.residuals  import RESIDUALS
from obslib.idis.reader     import READER
from obslib.idis.theory     import THEORY
from qcdlib                 import aux, eweak, pdf, alphaS, mellin
from nuclib                 import deuterium,helium
from fitlib                 import resman,parman
from tools.config           import conf,load_config
from obslib.idis.theory     import HT,OFFSHELL,THEORY


def get_tabs():
    conf['datasets']={}
    conf['datasets']['idis']={}
    #------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['filters']=[]
    #------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['filters'].append("Q2>1.27**2") 
    conf['datasets']['idis']['filters'].append("W2>3.0") 
    #------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx']={}
    #------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx'][10010]='idis/expdata/10010.xlsx' # proton   | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10016]='idis/expdata/10016.xlsx' # proton   | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10020]='idis/expdata/10020.xlsx' # proton   | F2            | NMC                   
    conf['datasets']['idis']['xlsx'][10026]='idis/expdata/10026.xlsx' # proton   | sigma red     | HERA II NC e+ (1)     
    conf['datasets']['idis']['xlsx'][10027]='idis/expdata/10027.xlsx' # proton   | sigma red     | HERA II NC e+ (2)     
    conf['datasets']['idis']['xlsx'][10028]='idis/expdata/10028.xlsx' # proton   | sigma red     | HERA II NC e+ (3)     
    conf['datasets']['idis']['xlsx'][10029]='idis/expdata/10029.xlsx' # proton   | sigma red     | HERA II NC e+ (4)     
    conf['datasets']['idis']['xlsx'][10030]='idis/expdata/10030.xlsx' # proton   | sigma red     | HERA II NC e-         
    conf['datasets']['idis']['xlsx'][10031]='idis/expdata/10031.xlsx' # proton   | sigma red     | HERA II CC e+         
    conf['datasets']['idis']['xlsx'][10032]='idis/expdata/10032.xlsx' # proton   | sigma red     | HERA II CC e-         
    conf['datasets']['idis']['xlsx'][10003]='idis/expdata/10003.xlsx' # proton   | sigma red     | JLab Hall C (E00-106) 
    conf['datasets']['idis']['xlsx'][10007]='idis/expdata/10007.xlsx' # proton   | sigma red     | HERMES                
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx'][10011]='idis/expdata/10011.xlsx' # deuteron | F2            | SLAC                  
    conf['datasets']['idis']['xlsx'][10017]='idis/expdata/10017.xlsx' # deuteron | F2            | BCDMS                 
    conf['datasets']['idis']['xlsx'][10021]='idis/expdata/10021.xlsx' # d/p      | F2d/F2p       | NMC                   
    conf['datasets']['idis']['xlsx'][10006]='idis/expdata/10006.xlsx' # deuteron | F2            | HERMES               
    conf['datasets']['idis']['xlsx'][10002]='idis/expdata/10002.xlsx' # deuteron | F2            | JLab Hall C (E00-106) 
    conf['datasets']['idis']['xlsx'][10033]='idis/expdata/10033.xlsx' # n/d      | F2n/F2d       | BONUS                 
    #-------------------------------------------------------------------------------------------------------------------
    #conf['datasets']['idis']['xlsx'][10037]='idis/expdata/10037.xlsx' # proton   | sigmacc red   | HERA I and II charm rean. 
    #conf['datasets']['idis']['xlsx'][10038]='idis/expdata/10038.xlsx' # proton   | sigmabb red  | HERA I and II bottom rean.
    #conf['datasets']['idis']['xlsx'][10039]='idis/expdata/10039.xlsx' # proton   | F2cc          | EMC charm
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['xlsx'][10041]='idis/expdata/10041.xlsx' # h/d      | F2h/F2d       | Seely
    conf['datasets']['idis']['xlsx'][10050]='idis/expdata/10050.xlsx' # d/p      | F2d/F2p       | MARATHON
    conf['datasets']['idis']['xlsx'][10051]='idis/expdata/10051.xlsx' # h/t      | F2h/F2t       | MARATHON


    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm']={}
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm'][10010]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10016]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10020]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10003]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10007]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    #-------------------------------------------------------------------------------------------------------------------
    conf['datasets']['idis']['norm'][10011]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10017]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10006]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10002]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    conf['datasets']['idis']['norm'][10033]={'value':    1, 'min': 8.00000e-01, 'max': 1.20000e+00, 'fixed': False}
    #-------------------------------------------------------------------------------------------------------------------
    tabs=READER().load_data_sets('idis')
    return tabs

def test_interpolation(wdir):

    conf['order']   = 'NLO'
    conf['Q20']     = 1.27**2
    conf['aux']     = aux.AUX()
    conf['mellin']  = mellin.MELLIN(npts=4)
    conf['alphaS']  = alphaS.ALPHAS()
    conf['eweak']   = eweak.EWEAK()
    conf['pdf']     = pdf.PDF()
    conf['pdf-mom'] = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']    = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])
    conf['hsmf']    = helium.HELIUM     (path2nuctab='%s/nuclib/grids/helium'%os.environ['FITPACK'])

    conf['ht4']  = HT()
    conf['off']  = OFFSHELL()

    conf['tmc']      = 'AOT'
    conf['nuc']      = True
    conf['ht']       = True 
    conf['offshell'] = True 

    tabs=get_tabs()
    conf['idis tabs'] = tabs

    stfuncs=THEORY()
    conf['idis'] = stfuncs

    #--get random replica
    jar   = load('../../../analysis-hx/results/upol/step21/data/jar-21.dat')
    load_config('../../../analysis-hx/results/upol/step21/input.py')
    jar   = load('%s/data/jar-21.dat'%wdir)
    load_config('%s/input.py'%wdir)
    rman = resman.RESMAN(parallel=False,datasets=False)
    parman = rman.parman
    parman.order = jar['order']
    replicas = jar['replicas']
    par = replicas[0]
    parman.set_new_params(par,initial=True)
    residuals=RESIDUALS()
 
    nrows,ncols=1,1
    fig = py.figure(figsize=(ncols*12,nrows*7))
    ax11=py.subplot(nrows,ncols,1)

    divider = make_axes_locatable(ax11)
    axL = divider.append_axes("right",size=5,pad=0,sharey=ax11)
    axL.set_xlim(0.1,0.9)
    axL.spines['left'].set_visible(False)
    axL.yaxis.set_ticks_position('right')
    py.setp(axL.get_xticklabels(),visible=True)

    ax11.spines['right'].set_visible(False)


    #--get interpolated observables
    conf['idis grid'] = {}
    stfuncs.setup_interpolation()
    stfuncs._update()
    residuals.get_theory()
    approx = {}
    for idx in tabs:
        approx[idx] = residuals.tabs[idx]['thy']
    #--get exact observables
    conf['idis grid']['xlsx'] = []
    for idx in tabs:
        if idx==10027: continue
        conf['idis grid']['xlsx'].append('%s/database/idis/expdata/%s.xlsx'%(os.environ['FITPACK'],idx))
    stfuncs.setup_interpolation()
    stfuncs._update()
    residuals.get_theory()
    exact = {}
    for idx in tabs:
        exact[idx] = residuals.tabs[idx]['thy']

    hand = {}
    for idx in tabs:
        X      = tabs[idx]['X']
        Q2     = tabs[idx]['Q2']
        target = tabs[idx]['target'][0]
        obs    = tabs[idx]['obs'][0].strip()
        alpha  = tabs[idx]['alpha']
        if idx==10016:   marker,color,label  = '^','black',    'BCDMS'
        elif idx==10017: marker,color        = '^','black'
        elif idx==10020: marker,color,label  = '+','goldenrod','NMC'
        elif idx==10021: marker,color        = '+','goldenrod'
        elif idx==10010: marker,color,label  = 'v','blue',     'SLAC'
        elif idx==10011: marker,color        = 'v','blue'
        elif idx==10026: marker,color,label  = 'o','green',    'HERA'
        elif idx==10027: marker,color        = 'o','green'
        elif idx==10028: marker,color        = 'o','green'
        elif idx==10029: marker,color        = 'o','green'
        elif idx==10030: marker,color        = 'o','green'
        elif idx==10031: marker,color        = 'o','green'
        elif idx==10032: marker,color        = 'o','green'
        elif idx==10033: marker,color,label  = 's','orange',   'JLab BONuS'
        elif idx==10002: marker,color,label  = 'x','red',      'JLab Hall C'
        elif idx==10003: marker,color        = 'x','red'
        #elif idx==10050: marker,color,label  = '+','cyan',     'MARATHON'
        #elif idx==10051: marker,color        = '+','cyan'
        #elif idx==10052: marker,color        = '+','cyan'
        else: continue

        rel_thy  = np.abs((exact[idx]-approx[idx]))/exact[idx]*100
        rel_exp  = np.abs(alpha)/tabs[idx]['value']*100
        rat = rel_thy/rel_exp

        X1, XL = [],[]
        rat1, ratL = [],[]
        for i in range(len(X)):
            if X[i] < 0.1:
                X1.append(X[i])
                rat1.append(rat[i])
            if X[i] > 0.1:
                XL.append(X[i])
                ratL.append(rat[i])
        
        hand[label] = ax11.scatter(X1,rat1,marker=marker,color=color,s=10)
        axL.scatter(XL,ratL,marker=marker,color=color,s=25)

    ax11.set_xscale('log')
    ax11.set_yscale('log')

    ax11.axhline(1,0,1,alpha=0.5,ls='--',color='black')
    axL. axhline(1,0,1,alpha=0.5,ls='--',color='black')
    ax11.tick_params(axis='both',which='both',top=True,right=False,direction='in',labelsize=30)
    axL.tick_params(axis='both',which='both',top=True,right=True,labelright=False,direction='in',labelsize=30)
    ax11.set_xlim(2e-5,0.1)
    ax11.set_xticks([1e-4,1e-3,1e-2])
    ax11.set_xticklabels([r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'])
    axL.set_xticks([0.1,0.3,0.5,0.7])
    ax11.set_ylim(1e-6,5)

    ax11.axvline(0.1,color='black',ls=':',alpha=0.5)

    axL.set_xlabel(r'\boldmath$x$',size=40)
    axL.xaxis.set_label_coords(0.95,0.00)
    ax11.set_ylabel(r'\boldmath$ \delta_{\rm{rel}}^{\rm{thy}} / \delta_{\rm{rel}}^{\rm{exp}}$',size=30)

    handles,labels = [], []
    handles.append(hand['BCDMS'])
    handles.append(hand['NMC'])
    handles.append(hand['SLAC'])
    handles.append(hand['JLab BONuS'])
    handles.append(hand['JLab Hall C'])
    handles.append(hand['HERA'])
    #handles.append(hand['MARATHON'])
    labels.append(r'\textbf{\textrm{BCDMS}}')
    labels.append(r'\textbf{\textrm{NMC}}')
    labels.append(r'\textbf{\textrm{SLAC}}')
    labels.append(r'\textbf{\textrm{JLab BONuS}}')
    labels.append(r'\textbf{\textrm{JLab Hall C}}')
    labels.append(r'\textbf{\textrm{HERA}}')
    #labels.append(r'\textbf{\textrm{MARATHON}}')
    ax11.legend(handles,labels,loc='lower left',fontsize=20,frameon=False, handlelength = 1.0, handletextpad = 0.1)

    py.tight_layout()
    checkdir('gallery')
    py.savefig('gallery/interp.png')
    print('Saving image to gallery/interp.png')

def test_theory():
    
    conf['order']      = 'NLO'
    conf['Q20']        = 1.27**2
    conf['dglap mode'] ='truncated'
    conf['aux']        = aux.AUX()
    conf['mellin']     = mellin.MELLIN(npts=4)
    conf['alphaS']     = alphaS.ALPHAS()
    conf['eweak']      = eweak.EWEAK()
    conf['pdf']        = pdf.PDF()
    conf['pdf-mom']    = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']       = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']        = HT()
    conf['off']        = OFFSHELL()

    conf['tmc']        = 'AOT'
    conf['nuc']        = True
    conf['ht']         = True 
    conf['offshell']   = True 

    X=np.linspace(0.01,0.1,4000)
    Q2=np.ones(X.size)*10.0

    thy=THEORY()
    thy.setup_interpolation()

    #--no parallel (only option)
    t1=time.time()
    for i in range(10):
        lprint('%d/%d'%(i+1,10))
        thy._update()
        thy.get_stf(X,Q2,stf='F2',tar='p')
        thy.get_stf(X,Q2,stf='F2',tar='n')
        thy.get_stf(X,Q2,stf='F2',tar='d')
    t2=time.time()
    print('\nno parallel: t=%0.2f'%(t2-t1))

def test_smearing(dsmf_group = 'paris', hsmf_group = 'kpsv'):
   
    #--print out integrated smearing functions from 0 to MA/MN
    conf['order']      = 'NLO'
    conf['Q20']        = 1.27**2
    conf['dglap mode'] ='truncated'
    conf['aux']        = aux.AUX()
    conf['alphaS']     = alphaS.ALPHAS()
    conf['eweak']      = eweak.EWEAK()
    conf['pdf-mom']    = pdf.PDF(mellin.IMELLIN())
    conf['mellin']     = mellin.MELLIN(npts=4)
    conf['dsmf']       = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'],dsmf_type = dsmf_group)
    conf['hsmf']       = helium.HELIUM     (path2nuctab = '%s/nuclib/grids/helium'%os.environ['FITPACK'],group = hsmf_group)
    X   = np.linspace(0,0.99,1)
    idis=THEORY()

    Q2 = np.array([10000000000])
    gX = idis.gX
    gW = idis.gW
    XM,  gXM = np.meshgrid(X,gX)
    Q2M, gWM = np.meshgrid(Q2,gW)
    a = np.zeros(np.shape(XM))

    dsmf = idis.dsmf
    hsmf = idis.hsmf

    shape = np.ones(np.shape(XM))
    #--deuterium
    b = idis.ymaxD
    YM = 0.5*(b-a)*gXM + 0.5*(a+b)
    JM = 0.5*(b-a)
    dsmfon   = dsmf.get_fXX2('f22' ,'onshell' ,XM,Q2M,YM)*shape
    dsmfoff  = dsmf.get_fXX2('f22' ,'offshell',XM,Q2M,YM)*shape

    dsmfon   = np.einsum('ij,ij,ij->j',gWM,JM,dsmfon)  [0]
    dsmfoff  = np.einsum('ij,ij,ij->j',gWM,JM,dsmfoff) [0]

    #--helium
    b = idis.ymaxH
    YM = 0.5*(b-a)*gXM + 0.5*(a+b)
    JM = 0.5*(b-a)
    hsmfpon  = hsmf.get_fXX2('f22p','onshell' ,XM,Q2M,YM)*shape
    hsmfpoff = hsmf.get_fXX2('f22p','offshell',XM,Q2M,YM)*shape
    hsmfnon  = hsmf.get_fXX2('f22n','onshell' ,XM,Q2M,YM)*shape
    hsmfnoff = hsmf.get_fXX2('f22n','offshell',XM,Q2M,YM)*shape

    hsmfpon  = np.einsum('ij,ij,ij->j',gWM,JM,hsmfpon) [0]
    hsmfpoff = np.einsum('ij,ij,ij->j',gWM,JM,hsmfpoff)[0]
    hsmfnon  = np.einsum('ij,ij,ij->j',gWM,JM,hsmfnon) [0]
    hsmfnoff = np.einsum('ij,ij,ij->j',gWM,JM,hsmfnoff)[0]

    print('Smearing Functions:')
    print('Deuterium: %s'%dsmf_group)
    print('Deuterium on-shell: %4.3f'      %dsmfon)
    print('Deuterium off-shell: %4.3f'     %dsmfoff)
    print('Helium: %s'   %hsmf_group)
    print('Helium proton on-shell: %4.3f'  %hsmfpon)
    print('Helium neutron on-shell: %4.3f' %hsmfnon)
    print('Helium proton off-shell: %4.3f' %hsmfpoff)
    print('Helium neutron off-shell: %4.3f'%hsmfnoff)

def plot_grid():

    conf['aux']     = aux.AUX()
    conf['order']   = 'NLO'
    conf['mellin']  = mellin.MELLIN(npts=4)

    tabs=get_tabs()
    conf['idis tabs'] = tabs

    stfuncs=THEORY()
    conf['idis'] = stfuncs

    nrows,ncols=1,3
    fig = py.figure(figsize=(ncols*12,nrows*7))
    ax11=py.subplot(nrows,ncols,1)
    ax12=py.subplot(nrows,ncols,2)
    ax13=py.subplot(nrows,ncols,3)

    divider = make_axes_locatable(ax11)
    axL11 = divider.append_axes("right",size=5,pad=0,sharey=ax11)
    axL11.set_xlim(0.1,0.9)
    axL11.spines['left'].set_visible(False)
    axL11.yaxis.set_ticks_position('right')
    py.setp(axL11.get_xticklabels(),visible=True)

    divider = make_axes_locatable(ax12)
    axL12 = divider.append_axes("right",size=5,pad=0,sharey=ax12)
    axL12.set_xlim(0.1,0.9)
    axL12.spines['left'].set_visible(False)
    axL12.yaxis.set_ticks_position('right')
    py.setp(axL12.get_xticklabels(),visible=True)

    divider = make_axes_locatable(ax13)
    axL13 = divider.append_axes("right",size=5,pad=0,sharey=ax13)
    axL13.set_xlim(0.1,0.9)
    axL13.spines['left'].set_visible(False)
    axL13.yaxis.set_ticks_position('right')
    py.setp(axL13.get_xticklabels(),visible=True)

    ax11.spines['right'].set_visible(False)
    ax12.spines['right'].set_visible(False)
    ax13.spines['right'].set_visible(False)

    #--get interpolated observables
    conf['idis grid'] = {}
    stfuncs.setup_interpolation()

    hand = {}
    #--plot data
    for idx in tabs:
        X      = tabs[idx]['X']
        Q2     = tabs[idx]['Q2']
        target = tabs[idx]['target'][0]

        if target in ['p']:
            hand['data'] = ax11 .scatter(X,Q2,color='red'   ,s=15)
            axL11               .scatter(X,Q2,color='red'   ,s=15)
        if target in ['d']:
            hand['data'] = ax12 .scatter(X,Q2,color='blue'  ,s=15)
            axL12               .scatter(X,Q2,color='blue'  ,s=15)
        if target in ['d/p','n/d']:
            hand['data'] = ax11 .scatter(X,Q2,color='red'   ,s=15)
            axL11               .scatter(X,Q2,color='red'   ,s=15)
            hand['data'] = ax12 .scatter(X,Q2,color='blue'  ,s=15)
            axL12               .scatter(X,Q2,color='blue'  ,s=15)
        if target in ['h/t']:
            hand['data'] = ax13 .scatter(X,Q2,color='green' ,s=15)
            axL13               .scatter(X,Q2,color='green' ,s=15)
        if target in ['h/d']:
            hand['data'] = ax12 .scatter(X,Q2,color='blue'  ,s=15)
            axL12               .scatter(X,Q2,color='blue'  ,s=15)
            hand['data'] = ax13 .scatter(X,Q2,color='green' ,s=15)
            axL13               .scatter(X,Q2,color='green' ,s=15)

    #--plot grid
    for tar in ['p','d','h']:

        Xgrid  = conf['idis'].X[tar]
        Q2grid = conf['idis'].Q2[tar]

        if tar=='p':
            hand['grid'] = ax11.scatter(Xgrid,Q2grid,color='black',s=20)
            axL11              .scatter(Xgrid,Q2grid,color='black',s=20)
        if tar=='d':
            hand['grid'] = ax12.scatter(Xgrid,Q2grid,color='black',s=20)
            axL12              .scatter(Xgrid,Q2grid,color='black',s=20)
        if tar=='h':
            hand['grid'] = ax13.scatter(Xgrid,Q2grid,color='black',s=20)
            axL13              .scatter(Xgrid,Q2grid,color='black',s=20)

    for ax in [ax11, ax12, ax13]:
        ax.set_xscale('log')
        ax.set_yscale('log')
        ax.tick_params(axis='both',which='both',top=True,right=False,direction='in',labelsize=30)
        ax.set_xlim(1e-5,0.1)
        ax.set_xticks([1e-4,1e-3,1e-2])
        ax.set_xticklabels([r'$10^{-4}$',r'$10^{-3}$',r'$10^{-2}$'])
        ax.set_ylim(0.8,3e5)
        ax.axvline(0.1,color='black',ls=':',alpha=0.5)
        ax.set_ylabel(r'\boldmath$ \delta_{\rm{rel}}^{\rm{thy}} / \delta_{\rm{rel}}^{\rm{exp}}$',size=30)

    for axL in [axL11, axL12, axL13]:
        axL.set_xlabel(r'\boldmath$x$',size=40)
        axL.xaxis.set_label_coords(0.95,0.00)
        axL.set_xticks([0.1,0.3,0.5,0.7])
        axL.tick_params(axis='both',which='both',top=True,right=True,labelright=False,direction='in',labelsize=30)
        axL.axvline(0.1,color='black',ls=':',alpha=0.5)

    ax11.text(0.05,0.50,r'\textrm{\textbf{p/n}}',transform=ax11.transAxes,size=40)
    ax12.text(0.05,0.50,r'\textrm{\textbf{d}}'  ,transform=ax12.transAxes,size=40)
    ax13.text(0.05,0.50,r'\textrm{\textbf{h/t}}',transform=ax13.transAxes,size=40)

    #handles,labels = [], []
    #ax11.legend(handles,labels,loc='lower left',fontsize=20,frameon=False, handlelength = 1.0, handletextpad = 0.1)

    py.tight_layout()
    checkdir('gallery')
    py.savefig('gallery/grid.png')
    print('Saving image to gallery/grid.png')

#--heavy quark tests (outdated)
def main03():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=4)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 

    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()
    X=0.005
    Q2=1.2
    #x=0.083
    #Q2=2.39
    #X=np.linspace(0.01,0.1,3)
    #Q2=np.ones(X.size)*10.0
    
    #print(thy.get_stf(X,Q2,stf='FLb',tar='p'))
    #print(thy.get_stf(X,Q2,stf='F2',tar='p'))
    #print(thy.get_stf(X,Q2,stf='F2',tar='n'))
    #print(thy.get_stf(X,Q2,stf='F2',tar='d'))
    #print(thy.get_stf(X,Q2,stf='F2b',tar='p'))

    #tmc='AOT'
    #print(thy._get_HQFXN(X,Q2,stf='F2',hq='b',tmc=False))
    #print(thy._get_HQFXN(X,Q2,stf='FL',hq='b',tmc=False))

    #print(thy.get_stf(X,Q2,stf='F2',tar='p'))
    #print(thy.get_stf(X,Q2,stf='F2',tar='n'))
    #print(thy.get_stf(X,Q2,stf='F2',tar='d'))

    #print(thy.get_stf(X,Q2,stf='F2c',tar='n'))
    #print(thy.get_stf(X,Q2,stf='F2c',tar='d'))

    #print(thy.get_stf(X,Q2,stf='F2b',tar='p'))
    #print(thy.get_stf(X,Q2,stf='F2b',tar='n'))
    #print(thy.get_stf(X,Q2,stf='F2b',tar='d'))

def main04():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=4)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 

    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()


    tabs=get_tabs()

    
    ax=py.subplot(111)
    for idx in tabs:
        X =tabs[idx]['X']
        Q2=tabs[idx]['Q2']
        IF2=thy.get_stf(X,Q2,stf='F2c',tar='p')
        #IFL=thy.get_stf(X,Q2,stf='FLc',tar='p')
        EF2=np.array([thy._get_HQFXN(X[i],Q2[i],stf='F2',hq='c',tmc=False) for i in range(X.size)])
        #EFL=np.array([thy._get_HQFXN(X[i],Q2[i],stf='FL',hq='c',tmc=False) for i in range(X.size)])
        relF2=np.abs((EF2-IF2)/EF2)*0.2
        for i in range(X.size):
        #    if np.isinf(relF2[i]): continue
        #    print(relF2[i])
            if EF2[i]> 0: c='r'
            if EF2[i]==0: c='k'
            if EF2[i]< 0: c='b'
            ax.plot(X[i],Q2[i],marker='o',color=c,markersize=relF2[i])
        #ax.plot(X,relF2,'k.')
        #ax.plot(X,(EF2),'kx')
        #ax.plot(X,(IF2),'r.')
        xmin=np.amin(X)
        xmax=np.amax(X)
        Q2min=np.amin(Q2)
        Q2max=np.amax(Q2)

    #ax.axhline(1)
    #ax.set_ylabel(r'$\rm rel-err~(\%)$',size=20)
    #ax.set_ylim(None,0.00001)
    ax.semilogy()
    ax.semilogx()
    ax.set_xlabel(r'$x_{\rm bj}$',size=20)
    ax.set_ylabel(r'$Q^2$',size=20)

    #ax=py.subplot(212)
    #X =thy.X
    #Q2=thy.Q2
    #IF2=thy.get_stf(X,Q2,stf='F2c',tar='p')
    ##IFL=thy.get_stf(X,Q2,stf='FLc',tar='p')
    #EF2=np.array([thy._get_HQFXN(X[i],Q2[i],stf='F2',hq='c',tmc=False) for i in range(X.size)])
    ##EFL=np.array([thy._get_HQFXN(X[i],Q2[i],stf='FL',hq='c',tmc=False) for i in range(X.size)])
    #relF2=np.abs((EF2-IF2)/EF2)*100
    #for i in range(X.size):
    #    if X[i]<xmin: continue
    #    if X[i]>xmax: continue
    #    if Q2[i]<Q2min: continue
    #    if Q2[i]>Q2max: continue
    #    ax.plot(X[i],np.abs(EF2[i]),'kx')
    #    ax.plot(X[i],np.abs(IF2[i]),'r.')
    ##    ax.plot(X[i],Q2[i],'ko')
    ##ax.set_xlabel(r'$x_{\rm bj}$',size=20)
    ##ax.set_ylabel(r'$\rm rel-err~(\%)$',size=20)
    #ax.semilogy()
    #ax.semilogx()


    py.tight_layout()
    py.savefig('test.pdf')

    return 

def main05():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=4)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 

    
    thy=THEORY()
    thy.setup_interpolation()

    tabs=get_tabs()
    
    ax=py.subplot(111)

    X =thy.X
    Q2=thy.Q2
    ax.plot(X,Q2,'ko')

    for idx in tabs:
        X =tabs[idx]['X']
        Q2=tabs[idx]['Q2']
        ax.plot(X,Q2,'ro')

    ax.semilogy()
    ax.semilogx()



    py.tight_layout()
    py.savefig('grid.png')

    return 

def main0X():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=16)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 
    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()

    #tabs=get_tabs()
    
    ax=py.subplot(111)

    X=10**np.linspace(-4,np.log10(0.8),100)
    for Q2 in [1.27**2,4,10,1000]:
        F2=np.array([X[i]*thy._get_FXN(X[i],Q2,stf='F2',nucleon='p') for i in range(X.size)])
        p,=ax.plot(X,F2,label=r'$Q^2=%0.2f$'%Q2)

        F2=np.array([X[i]*thy._get_HQFXN(X[i],Q2,stf='F2',hq='c',tmc=False) for i in range(X.size)])
        ax.plot(X,F2,ls=':',color=p.get_color())

        F2=np.array([X[i]*thy._get_HQFXN(X[i],Q2,stf='F2',hq='b',tmc=False) for i in range(X.size)])
        ax.plot(X,F2,ls='--',color=p.get_color())

    ax.text(0.1,0.3,r'$-~{\rm Inclusive}$',transform=ax.transAxes,size=20)
    ax.text(0.1,0.2,r'$...~{\rm Charm}$',transform=ax.transAxes,size=20)
    ax.text(0.1,0.1,r'$---~{\rm Bottom}$',transform=ax.transAxes,size=20)


    ax.legend(loc=2)
    ax.semilogy()
    ax.semilogx()
    ax.set_ylabel(r'$xF_2(x)$',size=20)
    ax.set_xlabel(r'$x$',size=20)

    py.tight_layout()
    py.savefig('test.pdf')

    return 

def main0X():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=16)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 
    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()

    #tabs=get_tabs()
    
    ax=py.subplot(111)

    X=10**np.linspace(-4,np.log10(0.8),100)
    for Q2 in [1.27**2,4,10,1000]:

        EF2=np.array([thy._get_HQFXN(X[i],Q2,stf='F2',hq='c',tmc=False) for i in range(X.size)])
        #p,=ax.plot(X,F2,ls='--',label=r'$Q^2=%0.2f$'%Q2)

        IF2=thy.get_stf(X,Q2,stf='F2c',tar='p')
        #ax.plot(X,X*F2,'k:')#ls='--')#,color=p.get_color())
        

        ax.plot(X,EF2/IF2,ls='',marker='.')#ls='--')#,color=p.get_color())


    ax.legend(loc=4)
    ax.semilogy()
    ax.semilogx()
    ax.set_ylabel(r'$\rm ratio$',size=20)
    ax.set_xlabel(r'$x$',size=20)

    py.tight_layout()
    py.savefig('test.pdf')

    return 

def main0X():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=4)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 
    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()

    tabs=get_tabs()
    
    ax=py.subplot(111)

    for idx in tabs:
        X =tabs[idx]['X']
        Q2=tabs[idx]['Q2']
        IF2=thy.get_stf(X,Q2,stf='F2c',tar='p')
        EF2=np.array([thy._get_HQFXN(X[i],Q2[i],stf='F2',hq='c',tmc=False) for i in range(X.size)])
        relF2=np.abs((EF2-IF2)/EF2)*0.2
        rat=EF2/IF2
        for i in range(X.size):
        #    if np.isinf(relF2[i]): continue
        #    print(relF2[i])
            #if EF2[i]> 0: c='r'
            #if EF2[i]==0: c='k'
            #if EF2[i]< 0: c='b'
            ax.plot(X[i],Q2[i],marker='o',markersize=rat[i])
        #ax.plot(X,EF2/IF2,'k.')
        #ax.plot(X,(EF2),'kx')
        #ax.plot(X,(IF2),'r.')
        xmin=np.amin(X)
        xmax=np.amax(X)
        Q2min=np.amin(Q2)
        Q2max=np.amax(Q2)

    X,Q2=thy.get_grid()

    #X =thy.X
    #Q2=thy.Q2
    ax.plot(X,Q2,'ko')

    ax.semilogy()
    ax.semilogx()

    ax.set_xlim(xmin,xmax)
    ax.set_ylim(Q2min,Q2max)

    py.tight_layout()
    py.savefig('test.pdf')

    return 

def main0X():

    conf['order']    = 'NLO'
    conf['Q20']      = 1.27**2
    conf['aux']      = aux.AUX()
    conf['mellin']   = mellin.MELLIN(npts=4)
    conf['alphaS']   = alphaS.ALPHAS()
    conf['eweak']    = eweak.EWEAK()
    conf['pdf']      = pdf.PDF()
    conf['pdf-mom']  = pdf.PDF(mellin.IMELLIN())
    conf['dsmf']     = deuterium.DEUTERON('%s/nuclib/grids/deuteron'%os.environ['FITPACK'])

    conf['ht4']      = HT()
    conf['off']      = OFFSHELL()

    conf['tmc']      = False#'AOT'
    conf['nuc']      = False#True
    conf['ht']       = False#True 
    conf['offshell'] = False#True 
    conf['hq']       = True 
    
    thy=THEORY()
    thy.setup_interpolation()
    thy.update()


    tabs=get_tabs()
    conf['idis tabs']=tabs
    conf['idis']=thy

    RESIDUALS()
   
    ax=py.subplot(111)

    for idx in tabs:
        X =tabs[idx]['X']
        Q2=tabs[idx]['Q2']

        IF2=thy.get_stf(X,Q2,stf='F2c',tar='p')
        EF2=np.array([thy._get_HQFXN(X[i],Q2[i],stf='F2',hq='c',tmc=False) for i in range(X.size)])
        F2=thy.get_stf(X,Q2,stf='F2c',tar='p')

        rel_thy=np.abs((EF2-IF2)/F2)
        rel_exp=tabs[idx]['alpha']/tabs[idx]['value']

        ax.plot(X,rel_thy/rel_exp,'r.')


    ax.set_ylabel(r'$(\delta F_{2,c}/F_2)/(\delta F_{2,c}^{\rm exp}/F_{2,c}^{\rm exp})$')

    
    ax.semilogy()
    ax.semilogx()
    ax.set_ylim(None,10)
    ax.axhline(1)

    ax.set_xlabel(r'$x$',size=20)


    py.tight_layout()
    py.savefig('test.pdf')

    return 

if __name__=="__main__":

    #dsmf_group = 'av18'
    #hsmf_group = 'ss'
    #test_smearing(dsmf_group,hsmf_group)

    plot_grid()


