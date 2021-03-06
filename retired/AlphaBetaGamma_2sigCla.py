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

#m_list = ['1500_1000','1500_1200','1600_1100','1700_1200','1800_1300','1900_100','1900_1000','1900_800','2200_100','2200_800']
b_list = ['DiLepTT','DY','QCD','SemiLepTT','SingleT','TTV','VV','WJ','Data']
others_bkg = ['DY','QCD','SingleT','TTV','VV','WJ']

CRs = ['SRLDM','SRHDM','CR1','CR2','CR3']

bkg_dires = []

bkg_dires = b_list

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
    
    args = parser.parse_args()
    infile = args.infile
    outdir = args.outdir

    if not os.path.exists(outdir):os.makedirs(outdir)
    if os.path.exists(outdir+'/alphabetagammaTable.tex'):
        os.remove(outdir+'/alphabetagammaTable.tex')
    if os.path.exists(outdir+'/alphabetagammaTable.txt'):
            os.remove(outdir+'/alphabetagammaTable.txt')

    latexalphabetagamma = open(outdir+'/alphabetagammaTable.tex', "a+")
    textalphabetagamma = open(outdir+'/alphabetagammaTable.txt', "a+")
    textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format('mass','alpha','alphaE','beta' ,'betaE', 'gamma','gammaE')+" \n")

    latexalphabetagamma.write("\\documentclass[a4paper,11pt]{article} \n")
    latexalphabetagamma.write("\\usepackage{longtable} \n")
    latexalphabetagamma.write("\\begin{document} \n")
    latexalphabetagamma.write("\\begin{longtable}{|c || c | c | c| c|} \n")
    latexalphabetagamma.write("\\hline \n")

    f = ROOT.TFile.Open(infile, "read")        
    
    bkg_list = bkg_dires
    
    otherSRLDM = ROOT.TH1F('otherSRLDM','otherSRLDM',50,0.0,1.0)
    otherSRHDM = ROOT.TH1F('otherSRHDM','otherSRHDM',50,0.0,1.0)
    otherCR1 = ROOT.TH1F('otherCR1','otherCR1',50,0.0,1.0)
    otherCR2 = ROOT.TH1F('otherCR2','otherCR2',50,0.0,1.0)
    otherCR3 = ROOT.TH1F('otherCR3','otherCR3',50,0.0,1.0)
    #otherCR4 = ROOT.TH1F('otherCR4','otherCR4',100,0.0,1.0)

    txtSRLDM  = open(outdir+'/SRLDM.txt', "w+")
    txtSRHDM  = open(outdir+'/SRHDM.txt', "w+")
    txtCR1 = open(outdir+'/CR1.txt', "w+")
    txtCR2 = open(outdir+'/CR2.txt', "w+")
    txtCR3 = open(outdir+'/CR3.txt', "w+")
    #txtCR4 = open(mass_dir+'/'+bkg_list[0].split('/')[1]+'_CR4.txt', "w+")
    txtalphabetagamma = open(outdir+'/alphabetagamma.txt', "w+")

    WJ_1   = 0 ; WJ_2   = 0 ; WJ_3   = 0 ; WJ_4   = 0 ; WJ_5   = 0
    TTl_1  = 0 ; TTl_2  = 0 ; TTl_3  = 0 ; TTl_4  = 0 ; TTl_5  = 0
    TTll_1 = 0 ; TTll_2 = 0 ; TTll_3 = 0 ; TTll_4 = 0 ; TTll_5 = 0 
    Data_1 = 0 ; Data_2 = 0 ; Data_3 = 0 ; Data_4 = 0 ; Data_5 = 0
    
    err_0   = ROOT.Double(0.) ; err_1   = ROOT.Double(0.) ; err_2   = ROOT.Double(0.) ; err_3   = ROOT.Double(0.) ; 
    err_4   = ROOT.Double(0.) ; err_5   = ROOT.Double(0.) ; err_6   = ROOT.Double(0.) ; err_7   = ROOT.Double(0.) ; 
    err_8   = ROOT.Double(0.) ; err_9   = ROOT.Double(0.) ; err_10  = ROOT.Double(0.) ; err_11  = ROOT.Double(0.) ; 
    err_12  = ROOT.Double(0.) ; err_13  = ROOT.Double(0.) ; err_14  = ROOT.Double(0.) ; err_15  = ROOT.Double(0.) ; 
    err_16  = ROOT.Double(0.) ; err_17  = ROOT.Double(0.) ; err_18  = ROOT.Double(0.) ; err_19  = ROOT.Double(0.) ; 
    err_20  = ROOT.Double(0.) ; err_21  = ROOT.Double(0.) ; err_22  = ROOT.Double(0.) ; err_23  = ROOT.Double(0.) ; 
    for bkg in bkg_list : 
        tdir = f.Get(bkg)
        hList = tdir.GetListOfKeys()
        for i in range(0,len(hList)):
            hname = (hList)[i].GetName()
            hist = tdir.Get(hname)
            if not 'signal_LDM' in hname : continue 
            for obkg in others_bkg : 
                if (obkg in bkg and 'SRLDM' in hname ): otherSRLDM.Add(hist)
                if (obkg in bkg and 'SRHDM' in hname ): otherSRHDM.Add(hist)
                if (obkg in bkg and 'CR1' in hname ): otherCR1.Add(hist)
                if (obkg in bkg and 'CR2' in hname ): otherCR2.Add(hist)
                if (obkg in bkg and 'CR3' in hname ): otherCR3.Add(hist)
                #if (obkg in bkg and '_CR4' in hname ): otherCR4.Add(hist)
            #print(hname, hist.Integral())
            
            if ('_SRLDM' in hname and not 'Data' in bkg): 
                txtSRLDM.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2), round(err_0, 2))+"\n")
                if 'SemiLepTT' in bkg  : TTl_4  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_10), 2)
                if 'DiLepTT' in bkg    : TTll_4 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_11), 2)
                
            elif ('_SRHDM' in hname and not 'Data' in bkg): 
                txtSRHDM.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2), round(err_0, 2))+"\n")
                if 'SemiLepTT' in bkg  : TTl_5  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_12), 2)
                if 'DiLepTT' in bkg    : TTll_5 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_13), 2)

            elif '_CR1' in hname : 
                txtCR1.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2), round(err_0, 2))+"\n")
                if 'SemiLepTT' in bkg  : TTl_1  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_1), 2)
                if 'DiLepTT' in bkg    : TTll_1 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_2), 2)
                if 'Data' in bkg       : Data_1 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_3), 2)
            elif '_CR2' in hname : 
                txtCR2.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2), round(err_0, 2))+"\n")
                if 'SemiLepTT' in bkg : TTl_2  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_4), 2)
                if 'DiLepTT' in bkg   : TTll_2 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_5), 2)
                if 'Data' in bkg      : Data_2 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_6), 2)
            elif '_CR3' in hname : 
                txtCR3.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2), round(err_0, 2))+"\n")
                if 'SemiLepTT' in bkg  : TTl_3  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_7), 2)
                if 'DiLepTT' in bkg    : TTll_3 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_8), 2)
                if 'Data' in bkg       : Data_3 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_9), 2)
        
    txtSRLDM.write((60 *('='))+'\n')
    txtSRLDM.write("{:<20}{:<20}{:<20}".format(otherSRLDM.GetTitle(),round(otherSRLDM.IntegralAndError(0, otherSRLDM.GetNbinsX()+1, err_14), 2), round(err_14, 2))+"\n")
    WJ_4 = round(otherSRHDM.IntegralAndError(0, otherSRLDM.GetNbinsX()+1, err_14),2)
        
    txtSRHDM.write((60 *('='))+'\n')
    txtSRHDM.write("{:<20}{:<20}{:<20}".format(otherSRHDM.GetTitle(),round(otherSRHDM.IntegralAndError(0, otherSRHDM.GetNbinsX()+1, err_15), 2), round(err_15, 2))+"\n")
    WJ_5 = round(otherSRHDM.IntegralAndError(0, otherSRHDM.GetNbinsX()+1, err_15),2)

    txtCR1.write((60 *('='))+'\n')
    txtCR1.write("{:<20}{:<20}{:<20}".format(otherCR1.GetTitle(),round(otherCR1.IntegralAndError(0, otherCR1.GetNbinsX()+1, err_16), 2), round(err_16, 2))+"\n")
    WJ_1 = round(otherCR1.IntegralAndError(0, otherCR1.GetNbinsX()+1, err_16),2)
        
    txtCR2.write((60 *('='))+'\n')
    txtCR2.write("{:<20}{:<20}{:<20}".format(otherCR2.GetTitle(),round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_17), 2), round(err_17, 2))+"\n")
    WJ_2 = round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_17),2)

    txtCR3.write((60 *('='))+'\n')
    txtCR3.write("{:<20}{:<20}{:<20}".format(otherCR3.GetTitle(),round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_18), 2), round(err_18, 2))+"\n")
    WJ_3 = round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_18),2)
        
    del otherSRLDM,otherSRHDM,otherCR1,otherCR2,otherCR3

    import scipy

    a = unumpy.umatrix([[TTl_1,TTll_1,WJ_1],[TTl_2,TTll_2,WJ_2],[TTl_3,TTll_3,WJ_3]],[[err_1,err_2,err_16],[err_4,err_5,err_17],[err_7,err_8,err_18]])
        
    b = unumpy.umatrix([[Data_1],[Data_2],[Data_3]],[[err_3],[err_6],[err_9]])

    #fun = lambda x: np.linalg.norm(np.dot(a,x)-b)
    #Y = minimize(fun, np.zeros(n), method='L-BFGS-B', bounds=[(0.,None) for x in range(n)])
    #Y = Y['x']           
    #Y = np.linalg.solve(a,b)
    Y_ = a.I * b 
    Y = np.squeeze(np.asarray(unumpy.nominal_values(Y_)))
    Yerr = np.squeeze(np.asarray(unumpy.std_devs(Y_)))

        #print (np.allclose(np.dot(a, Y), b))
    alpha, beta , gamma  = round(Y[0],2) , round(Y[1],2),round( Y[2],2)
    alphaerr, betaerr , gammaerr  = round(Yerr[0],2) , round(Yerr[1],2),round( Yerr[2],2)
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' '             ,'SemiLepTT' ,' '                       ,'DiLepTT' ,' '                       ,'others' ,' '                       ,'Data' ,' ')+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('SemiLepTTCat'  ,TTl_1       ,'+/-'+str(round(err_1,2)) ,TTll_1    ,'+/-'+str(round(err_2,2)) ,WJ_1     ,'+/-'+str(round(err_16,2)) ,Data_1 ,'+/-'+str(round(err_3,2)))+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('DiLepTTCat'    ,TTl_2       ,'+/-'+str(round(err_4,2)) ,TTll_2    ,'+/-'+str(round(err_5,2)) ,WJ_2     ,'+/-'+str(round(err_17,2)) ,Data_2 ,'+/-'+str(round(err_6,2 )))+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('othersCat'     ,TTl_3       ,'+/-'+str(round(err_7,2)) ,TTll_3    ,'+/-'+str(round(err_8,2)) ,WJ_3     ,'+/-'+str(round(err_18,2)) ,Data_3 ,'+/-'+str(round(err_9,2)))+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('SigCatLDM'     ,TTl_4       ,'+/-'+str(round(err_10,2)),TTll_4    ,'+/-'+str(round(err_11,2)),WJ_4     ,'+/-'+str(round(err_14,2)) ,'---' ,'+/-'+'---')+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('SigCatHDM'     ,TTl_5       ,'+/-'+str(round(err_12,2)),TTll_5    ,'+/-'+str(round(err_13,2)),WJ_5     ,'+/-'+str(round(err_15,2)) ,'---' ,'+/-'+'---')+"\n")


    txtalphabetagamma.write((60 *('='))+'\n')

    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ','alpha',' ','beta' ,' ', 'gamma',' ')+"\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ',alpha, '+/-'+str(alphaerr),beta ,'+/-'+str(betaerr), gamma,'+/-'+str(gammaerr) )+"\n")

        
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('      &'   ,'SemiLepTT' ,' &'                       ,'DiLepTT' ,' &'                       ,'others' ,' &'                       ,'Data' ,' ')+" \\\\ \n")
    latexalphabetagamma.write('\\hline\\hline \n')

    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}".format('SemiLepTTCat & '  ,TTl_1       ,' $\pm$ '+str(round(err_1,2)) +'& ' ,TTll_1    ,' $\pm$ '+str(round(err_2,2)) +'& ' ,WJ_1     ,' $\pm$ '+str(round(err_16,2))+'& ' ,Data_1)+"\\\\ \n")
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}".format('DiLepTTCat & '    ,TTl_2       ,' $\pm$ '+str(round(err_4,2)) +'& ' ,TTll_2    ,' $\pm$ '+str(round(err_5,2)) +'& ' ,WJ_2     ,' $\pm$ '+str(round(err_17,2)) +'& ',Data_2)+"\\\\ \n")
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}".format('othersCat &'      ,TTl_3       ,' $\pm$ '+str(round(err_7,2)) +'& ' ,TTll_3    ,' $\pm$ '+str(round(err_8,2)) +'& ' ,WJ_3     ,' $\pm$ '+str(round(err_18,2)) +'& ',Data_3)+"\\\\ \n")
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}".format('SigCat &'         ,TTl_4       ,' $\pm$ '+str(round(err_10,2))+'& ' ,TTll_4    ,' $\pm$ '+str(round(err_11,2))+'& ' ,WJ_4     ,' $\pm$ '+str(round(err_14,2)) +'& ','---')+"\\\\ \n")
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}{:<12}".format('SigCat &'         ,TTl_5       ,' $\pm$ '+str(round(err_12,2))+'& ' ,TTll_5    ,' $\pm$ '+str(round(err_13,2))+'& ' ,WJ_5     ,' $\pm$ '+str(round(err_15,2)) +'& ','---')+"\\\\ \n")

    latexalphabetagamma.write('\\hline\\hline \n')

    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('   & ','$\\alpha$',' & ','$\\beta$' ,' & ', '$\\gamma$','   & ')+"\\\\ \n")
    latexalphabetagamma.write('\\hline\\hline \n')
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('   & ',alpha, ' $\pm$ '+str(alphaerr)+' & ',beta ,' $\pm$ '+str(betaerr)+' & ', gamma,' $\pm$ '+str(gammaerr)+' &   \\\\' )+"\n")
    latexalphabetagamma.write('\\hline\\hline\\hline \n')
    textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format('   ' ,alpha, str(alphaerr) ,beta ,str(betaerr), gamma,str(gammaerr) )+"\n")
    
    latexalphabetagamma.write("\\caption{ Normalizaton parameters for backgrounds, calculated from control categories.}")
    latexalphabetagamma.write("\\label{normalize}")
    latexalphabetagamma.write("\\end{longtable}\n")
    latexalphabetagamma.write("\\end{document} \n")