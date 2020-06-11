#!/usr/bin/env python

import sys
sys.argv.append( '-b-' )
import ROOT
from ROOT import std
ROOT.gROOT.SetBatch(True)
sys.argv.remove( '-b-' )

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from uncertainties import unumpy
import numpy as np 
import os 
import datetime

from math import hypot, sqrt, ceil

import argparse
import htcondor

m_list = ['1500_1000','1500_1200','1600_1100','1700_1200','1800_1300','1900_100','1900_1000','1900_800','2200_100','2200_800']
b_list = ['DiLepTT','DY','QCD','SemiLepTT','SingleT','TTV','VV','WJ','Data']
others_bkg = ['DiLepTT','SemiLepTT','DY','SingleT','TTV','VV'] #keep the QCD out 'QCD',

CRs = ['SR','CR2','CR3','QCDCR']

bkg_dires = []

for m in m_list : 
    small_list = []
    for b in b_list:
        small_list.append(b+'/'+m)
    bkg_dires.append(small_list)

def hadd1ds(histList):
    '''  A functon to hadd background and set it's style '''
    sumbkg = ROOT.TH1F(histList[0].Clone())
    sumbkg.Reset()
    for bkghist in histList :
        sumbkg.Add(bkghist)
    #sumbkg.Draw('goff')
    sumbkg.SetTitle('Total BKG')
    sumbkg.SetName('sumbkg')
    return sumbkg

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--infile', help='List of datasets to process', metavar='infile')
    parser.add_argument('--outdir', help='output dir', metavar='outdir')
    parser.add_argument('--syst', help='systematics',default="nom", metavar='syst')
    
    
    args = parser.parse_args()
    infile = args.infile
    syst = args.syst
    outdir = os.path.join(args.outdir,syst)
     

    if not os.path.exists(outdir):os.makedirs(outdir)
    if os.path.exists(outdir+'/alphabetagammaTable.tex'):
        os.remove(outdir+'/alphabetagammaTable.tex')
    if os.path.exists(outdir+'/alphabetagammaTable.txt'):
            os.remove(outdir+'/alphabetagammaTable.txt')

    
    textalphabetagamma = open(outdir+'/alphabetagammaTable.txt', "a+")
    textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format('mass','alpha','alphaE','beta' ,'betaE', 'delta','deltaE')+" \n")

    f = ROOT.TFile.Open(infile, "read")        
    mass_dir = outdir
    if not os.path.exists(mass_dir):os.makedirs(mass_dir)
    
    otherSR = ROOT.TH1F('otherSR','otherSR',10000,0.0,1.0)
    otherCR1 = ROOT.TH1F('otherCR1','otherCR1',10000,0.0,1.0)
    otherCR2 = ROOT.TH1F('otherCR2','otherCR2',10000,0.0,1.0)
    otherCR3 = ROOT.TH1F('otherCR3','otherCR3',10000,0.0,1.0)
    otherQCDCR = ROOT.TH1F('otherQCDCR','otherQCDCR',10000,0.0,1.0)

    txtSR  = open(mass_dir+'/SR.txt', "w+")
    txtCR1 = open(mass_dir+'/CR1.txt', "w+")
    txtCR2 = open(mass_dir+'/CR2.txt', "w+")
    txtCR3 = open(mass_dir+'/CR3.txt', "w+")
    txtQCDCR = open(mass_dir+'/QCDCR.txt', "w+")


    WJ_1 = 0   ; WJ_2 = 0   ; WJ_3 = 0   ; WJ_QCDCR = 0
    TTJ_1 = 0  ; TTJ_2 = 0  ; TTJ_3 = 0  ; TTJ_QCDCR = 0
    data_3 = 0 ; data_2 = 0 ; data_3 = 0 ; data_QCDCR = 0 
    QCD_1 =0   ; QCD_2 = 0  ; QCD_3 = 0  ; QCD_QCDCR = 0 

    err   = ROOT.Double(0.) ; err_1 = ROOT.Double(0.) ; err_3 = ROOT.Double(0.)  ;err_4 = ROOT.Double(0.)  ;err_5 = ROOT.Double(0.)  ;err_6 = ROOT.Double(0.) 
    err_7 = ROOT.Double(0.) ; err_8 = ROOT.Double(0.) ; err_9 = ROOT.Double(0.)  ;err_10 = ROOT.Double(0.) ;err_11 = ROOT.Double(0.) ;err_12 = ROOT.Double(0.) 
    err_13 = ROOT.Double(0.) ; err_14 = ROOT.Double(0.) ; err_15 = ROOT.Double(0.)  ;err_16 = ROOT.Double(0.) ;err_17 = ROOT.Double(0.) ;err_18 = ROOT.Double(0.)
    err_19 = ROOT.Double(0.) ; err_20 = ROOT.Double(0.) ; err_21 = ROOT.Double(0.) ; err_22 = ROOT.Double(0.); err_23 = ROOT.Double(0.) ; err_24 = ROOT.Double(0.) 
    err_25 = ROOT.Double(0.) ; err_2 = ROOT.Double(0.)  ;err_0   = ROOT.Double(0.) ;err_00   = ROOT.Double(0.);err_000   = ROOT.Double(0.)
    bkg_list = b_list
    df_dict = {}
    dftex_dict = {}
    for bkg in bkg_list : 
        #print(bkg)
        tdir = f.Get(bkg)
        hList = tdir.GetListOfKeys()
        for i in range(0,len(hList)):
            hname = (hList)[i].GetName()
            if "Data" in bkg : 
                if "Jec" in syst : 
                    if not '_'+syst in hname: continue
                elif not '_nom' in hname: continue
            else : 
                if not '_'+syst in hname: continue
            hist = tdir.Get(hname)
            if not 'sig_0b' in hname : continue
            for obkg in others_bkg : 
                if (obkg in bkg and '_SR' in hname ): otherSR.Add(hist)
                if (obkg in bkg and '_CR1' in hname ): otherCR1.Add(hist)
                if (obkg in bkg and '_CR2' in hname ): otherCR2.Add(hist)
                if (obkg in bkg and '_CR3' in hname ): otherCR3.Add(hist)
                if (obkg in bkg and '_QCDCR' in hname ): otherQCDCR.Add(hist) 
            
            if ('_SR' in hname and not 'Data' in bkg): 
                txtSR.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                if 'WJ' in bkg   : WJ_1  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2)
                if 'QCD' in bkg  : QCD_1 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_15), 2)
            elif '_CR1' in hname : 
                txtCR1.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
            elif '_CR2' in hname : 
                txtCR2.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                if 'WJ' in bkg   : WJ_2  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_1), 2)
                if 'QCD' in bkg  : QCD_2 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_16), 2)
                if 'Data' in bkg : Data_2 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_3), 2)
            elif '_CR3' in hname : 
                txtCR3.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                if 'WJ' in bkg      : WJ_3  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_4), 2)
                if 'QCD' in bkg     : QCD_3 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_17), 2)
                if 'Data' in bkg    : Data_3 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_5), 2)
            elif '_QCDCR' in hname  : 
                txtQCDCR.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                if 'WJ' in bkg         : WJ_QCDCR  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_19), 2)
                if 'QCD' in bkg        : QCD_QCDCR = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_21), 2)
                if 'Data' in bkg       : Data_QCDCR = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_22), 2) 


    txtSR.write((60 *('='))+'\n')
    txtSR.write("{:<20}{:<20}{:<20}".format(otherSR.GetTitle(),round(otherSR.IntegralAndError(0, otherSR.GetNbinsX()+1, err_6), 2), round(err_6, 2))+"\n")
    TTJ_1 = round(otherSR.IntegralAndError(0, otherSR.GetNbinsX()+1, err_6),2)
    txtCR1.write((60 *('='))+'\n')
    txtCR1.write("{:<20}{:<20}{:<20}".format(otherCR1.GetTitle(),round(otherCR1.IntegralAndError(0, otherCR1.GetNbinsX()+1, err_7), 2), round(err_7, 2))+"\n")
    txtCR2.write((60 *('='))+'\n')
    txtCR2.write("{:<20}{:<20}{:<20}".format(otherCR2.GetTitle(),round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_8), 2), round(err_8, 2))+"\n")
    TTJ_2 = round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_8),2)
    txtCR3.write((60 *('='))+'\n')
    txtCR3.write("{:<20}{:<20}{:<20}".format(otherCR3.GetTitle(),round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_9), 2), round(err_9, 2))+"\n")
    TTJ_3 = round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_9),2)

    txtQCDCR.write((60 *('='))+'\n')
    txtQCDCR.write("{:<20}{:<20}{:<20}".format(otherQCDCR.GetTitle(),round(otherQCDCR.IntegralAndError(0, otherQCDCR.GetNbinsX()+1, err_24), 2), round(err_24, 2))+"\n")
    WJ_QCDCR = round(otherQCDCR.IntegralAndError(0, otherQCDCR.GetNbinsX()+1, err_24),2)
    
    del otherSR,otherCR1,otherCR2,otherCR3


    df_dict['Category'] = ["T1tttt Cat.", "TTJets Cat.","WJets Cat.","QCDCS"]
    df_dict['TTJets'] = [str(TTJ_1)+" +/- "+str(round(err_6,2)),str(TTJ_2)+" +/- "+str(round(err_8,2)),str(TTJ_3)+" +/- "+str(round(err_9,2)),str(TTJ_QCDCR)+" +/- "+str(round(err_24,2))]
    df_dict['Others'] = [str(WJ_1)+" +/- "+str(round(err_0,2)),str(WJ_2)+" +/- "+str(round(err_1,2)),str(WJ_3)+" +/- "+str(round(err_4,2)),str(WJ_QCDCR)+" +/- "+str(round(err_19,2))]
    df_dict['QCD'] = [str(QCD_1)+" +/- "+str(round(err_15,2)),str(QCD_2)+" +/- "+str(round(err_16,2)),str(QCD_3)+" +/- "+str(round(err_17,2)),str(QCD_QCDCR)+" +/- "+str(round(err_21,2))]
    df_dict['Data'] =   ["-- +/- --", str(Data_2)+" +/- "+str(round(err_3,2)), str(Data_3)+" +/- "+str(round(err_5,2)), str(Data_QCDCR)+" +/- "+str(round(err_22,2))]


    dftex_dict['Category &'] = ["T1tttt Cat. &", "TTJets Cat.&","WJets Cat.&","QCDCS &"]
    dftex_dict['TTJets &'] = [str(TTJ_1)+" $\pm$ "+str(round(err_6,2))+" &",str(TTJ_2)+" $\pm$ "+str(round(err_8,2))+" &",str(TTJ_3)+" $\pm$ "+str(round(err_9,2))+" &",str(TTJ_QCDCR)+" $\pm$ "+str(round(err_24,2))+" &"]
    dftex_dict['Others &'] = [str(WJ_1)+" $\pm$ "+str(round(err_0,2))+" &",str(WJ_2)+" $\pm$ "+str(round(err_1,2))+" &",str(WJ_3)+" $\pm$ "+str(round(err_4,2))+" &",str(WJ_QCDCR)+" $\pm$ "+str(round(err_19,2))+" &"]
    dftex_dict['QCD &'] = [str(QCD_1)+" $\pm$ "+str(round(err_15,2))+" &",str(QCD_2)+" $\pm$ "+str(round(err_16,2))+" &",str(QCD_3)+" $\pm$ "+str(round(err_17,2))+" &",str(QCD_QCDCR)+" $\pm$ "+str(round(err_21,2))+" &"]
    dftex_dict['Data \\\\ \\hline'] =   ["-- $\pm$ --"+" \\\\ \\hline", str(Data_2)+" $\pm$ "+str(round(err_3,2))+"  \\\\ \\hline", str(Data_3)+" $\pm$ "+str(round(err_5,2))+"  \\\\ \\hline", str(Data_QCDCR)+" $\pm$ "+str(round(err_22,2))+" \\\\"]


    import scipy


    a = unumpy.umatrix([[TTJ_2,WJ_2,QCD_2],[TTJ_3,WJ_3,QCD_3],[TTJ_QCDCR,WJ_QCDCR,QCD_QCDCR]],[[err_8,err_1,err_16],[err_9,err_4,err_17],[err_24,err_19,err_21]])
    
    b = unumpy.umatrix([[Data_2],[Data_3],[Data_QCDCR]],[[err_3],[err_5],[err_22]])

    #fun = lambda x: np.linalg.norm(np.dot(a,x)-b)
    #Y = minimize(fun, np.zeros(n), method='L-BFGS-B', bounds=[(0.,None) for x in range(n)])
    #Y = Y['x']           
    #Y = np.linalg.solve(a,b)
    Y_ = a.I * b 
    Y = np.squeeze(np.asarray(unumpy.nominal_values(Y_)))
    Yerr = np.squeeze(np.asarray(unumpy.std_devs(Y_)))
    #print (np.allclose(np.dot(a, Y), b))
    alpha, beta , delta = round(Y[0],2) , round(Y[1],2), round(Y[2],2)
    alphaerr, betaerr , deltaerr  = round(Yerr[0],2) , round(Yerr[1],2),round(Yerr[2],2)


    import pandas as pd 
    df = pd.DataFrame(df_dict)
    dftex = pd.DataFrame(dftex_dict)

    txtalphabetagamma_ = mass_dir+'/alphabetagamma.txt'
    with open(txtalphabetagamma_, 'w') as f: df.to_string(f, col_space=10,index=False)
    txtalphabetagamma = open(txtalphabetagamma_, "a")
    txtalphabetagamma.write("\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ','alpha',' ','beta' ,' ', 'delta',' ')+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ',alpha, '+/-'+str(alphaerr),beta ,'+/-'+str(betaerr),delta,'+/-'+str(deltaerr) )+"\n")
    
    textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format("",alpha, str(alphaerr) ,beta ,str(betaerr),delta ,str(deltaerr))+"\n")

    latexalphabetagamma_ = mass_dir+'/alphabetagammaTable.tex'
    latexalphabetagamma = open(latexalphabetagamma_, "a+")
    latexalphabetagamma.write("\\documentclass[a4paper,11pt]{article} \n")
    latexalphabetagamma.write("\\usepackage{longtable} \n")
    latexalphabetagamma.write("\\begin{document} \n")
    latexalphabetagamma.write("\\begin{longtable}{|c || c | c | c| c|} \n")
    latexalphabetagamma.write("\\hline \n")
    dftex.to_string(latexalphabetagamma, col_space=10,index=False)
    latexalphabetagamma.write('\\hline\\hline \n')

    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('   & ','$\\alpha$',' & ','$\\beta$' ,'&' ,'$\\delta$' , '& ')+"\\\\ \n")
    latexalphabetagamma.write('\\hline\\hline \n')
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('   & ',alpha, ' $\pm$ '+str(alphaerr)+' & ',beta ,' $\pm$ '+str(betaerr)+' & ',delta ,' $\pm$ '+str(deltaerr)+ ' & \\\\' )+"\n")
    latexalphabetagamma.write('\\hline\\hline\\hline \n')
    
    latexalphabetagamma.write("\\caption{ Normalizaton parameters for backgrounds, calculated from control categories.}\n")
    latexalphabetagamma.write("\\label{normalize}\n")
    latexalphabetagamma.write("\\end{longtable}\n")
    latexalphabetagamma.write("\\end{document} \n")
