#!/usr/bin/env python

import sys
sys.argv.append( '-b-' )
import ROOT
from ROOT import std
ROOT.gROOT.SetBatch(True)
sys.argv.remove( '-b-' )

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)

import os ,glob
import datetime

from math import hypot, sqrt, ceil

currentDT = datetime.datetime.now()

import argparse

channels = [['T1tttt_Scan','background','QCD']]#,['T1tttt_Scan_','background']]
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
    
def simpleAsimov(s,b):
    sign = s/sqrt(s+b+(0.2*b)**2)
    return sign

def getSFs(sffile,which='alpha') : 
    token = open(sffile,'r')
    linestoken=token.readlines()
    if which == 'alpha' :  tokens_column_number = 0
    if which == 'alphaE' :  tokens_column_number = 1
    if which == 'beta' :  tokens_column_number = 2
    if which == 'betaE' :  tokens_column_number = 3
    if which == 'gamma' :  tokens_column_number = 4
    if which == 'gammaE' :  tokens_column_number = 5
    if which == 'delta' :  tokens_column_number = 6
    if which == 'deltaE' :  tokens_column_number = 7
    resulttoken=[]
    #masstoken=[]
    for x in linestoken:
        resulttoken.append(x.split()[tokens_column_number])
        #masstoken.append(x.split()[0])
    token.close()
    #idx = masstoken.index(mass)
    return resulttoken[1]


if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='Runs a NAF batch system for nanoAOD', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--indir', help='List of Signals to process', metavar='indirS')
    #prser.add_argument('--indirB', help='List of BKGs to process', metavar='indirB')
    parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
    parser.add_argument('--Limit','-L', help='calculate the limit after making the datacards',default=False, action='store_true')
    parser.add_argument('--cmsswdir', help='cmssw directory',default='/nfs/dust/cms/user/amohamed/susy-desy/deepAK8/CMSSW_9_4_11/src/', metavar='cmsswdir')
    parser.add_argument('--execu', help="wight directory",default='./batch/Limit_exec.sh' ,metavar='execu')
    parser.add_argument('--sfs', help="the text file that has alpha beta gamma",default=None ,metavar='sfs')
    parser.add_argument('--year','--Y', help="Year of the data taken",default=None ,metavar="year")
    parser.add_argument('--oneBin','--ob',help="do one bin analysis or multibins",default=False , action='store_true')

    args = parser.parse_args()
    indirS = args.indir+'/scan'
    indirB = args.indir+'/grid'
    Limit = args.Limit
    cmsswdir = args.cmsswdir
    outdir = args.outdir
    execu = args.execu   
    sfs = args.sfs
    oneBin= args.oneBin

    if not os.path.exists(outdir): os.makedirs(outdir)
    
    datacardsdir = os.path.join(outdir,'datacards/combinedCards')
    if not os.path.exists(datacardsdir):os.makedirs(datacardsdir)       
    signal_List  = os.listdir(indirS)
    #background_List = os.listdir(indirB)
    text = open(outdir+"/background_comp.txt", "w+")

    for signal in signal_List : 
        systList = { 
                "Jec_Up"            :    []  ,
                "Jec_Down"          :    []  ,
                "btagSF_b_Up"       :    []  ,
                "btagSF_b_Down"     :    []  ,
                "btagSF_l_Up"       :    []  ,
                "btagSF_l_Down"     :    []  ,
                "ISR_Up"            :    []  ,
                "ISR_Down"          :    []  ,
                "lepSF_Up"          :    []  ,
                "lepSF_Down"        :    []  ,
                "PU_Up"             :    []  ,
                "PU_Down"           :    []  ,
                #"TTxsec_Up"         :    []  ,
                #"TTxsec_Down"       :    []  ,
                #"TTVxsec_Up"        :    []  ,
                #"TTVxsec_Down"      :    []  ,
                "Wpol_Up"           :    []  ,
                "Wpol_Down"         :    []  ,
                #"Wxsec_Up"          :    []  ,
                #"Wxsec_Down"        :    []
                }
        systListQCD = { 
                "Jec_Up"            :    []  ,
                "Jec_Down"          :    []  ,
                "btagSF_b_Up"       :    []  ,
                "btagSF_b_Down"     :    []  ,
                "btagSF_l_Up"       :    []  ,
                "btagSF_l_Down"     :    []  ,
                "ISR_Up"            :    []  ,
                "ISR_Down"          :    []  ,
                "lepSF_Up"          :    []  ,
                "lepSF_Down"        :    []  ,
                "PU_Up"             :    []  ,
                "PU_Down"           :    []  ,
                #"TTxsec_Up"         :    []  ,
                #"TTxsec_Down"       :    []  ,
                #"TTVxsec_Up"        :    []  ,
                #"TTVxsec_Down"      :    []  ,
                "Wpol_Up"           :    []  ,
                "Wpol_Down"         :    []  ,
                #"Wxsec_Up"          :    []  ,
                #"Wxsec_Down"        :    []
                }
        SigsystList = { 
                "Jec_Up"            :    []  ,
                "Jec_Down"          :    []  ,
                "btagSF_b_Up"       :    []  ,
                "btagSF_b_Down"     :    []  ,
                "btagSF_l_Up"       :    []  ,
                "btagSF_l_Down"     :    []  ,
                "ISR_Up"            :    []  ,
                "ISR_Down"          :    []  ,
                "lepSF_Up"          :    []  ,
                "lepSF_Down"        :    []  ,
                "Wpol_Up"           :    []  ,
                "Wpol_Down"         :    []  ,
                }

        if int(args.year) != 2016 : 
            systList.pop('ISR_Up', None)
            systList.pop('ISR_Down', None)
            SigsystList.pop('ISR_Up', None)
            SigsystList.pop('ISR_Down', None)
        intersection_syst = set(systList).intersection(set(SigsystList))
        
        #print(intersection_syst)

        lists = []
        signal_files = glob.glob(os.path.join(indirS,signal+'/*SR.root'))
        if not signal_files : 
            print("cannot find this file for ", signal,"will escape it, please check")
            continue
        sf = ROOT.TFile.Open(signal_files[0], "read")
        shist = sf.Get('Signal_1/'+signal+'/sig_SR_nom')
        for skey in SigsystList : 
            SigsystList[skey].append(sf.Get('Signal_1/'+signal+'/sig_SR_'+skey))

        lists.append(shist)
        mgo = float(signal.split('_')[0])
        mlsp = float(signal.split('_')[1])
        #bkgToUse_ = bkgToUse(mgo,mlsp)
        print (signal_files)
        bkg_files_SR = glob.glob(os.path.join(indirB+'/*SR.root'))
        bkg_files_CR1 = glob.glob(os.path.join(indirB+'/*CR1.root'))
        #print (bkg_files_SR)
        hist = ROOT.TH1F('background','background',10000,0.0,1.0)
        for key in systList : 
            systList[key].append(ROOT.TH1F('background_'+key,'background_'+key,10000,0.0,1.0))
        #print(systList)
        
        for bkgf in bkg_files_SR : 
            bkg_name = bkgf.split('/')[-1].split('_')[1]
            if bkg_name == 'Data' : continue 
            if bkg_name == 'QCD' : continue 
            bf = ROOT.TFile.Open(bkgf, "read")
            #print(bkg_name+'/'+bkg+'/'+bkg_name+bkg+'sig_SR_nominal')
            bhist = bf.Get(bkg_name+'/sig_SR_nom')
            for key in systList : 
                systList[key].append(bf.Get(bkg_name+'/sig_SR_'+key))
            #print(bhist.Integral())
            which = ''
            scalefactor = 1.0
            if bkg_name == 'DiLepTT' : which = 'beta'
            elif bkg_name == 'SemiLepTT' :which = 'alpha'
            elif bkg_name != 'Data' : which = 'gamma'
            scalefactor = float(getSFs(sfs,which=which))
            #print(bkg,bkgf)
            #print(bkg)
            bhist.Scale(scalefactor)
            for key in systList : 
                for i,syst in enumerate(systList[key]) : 
                    scalefactor_syst = float(getSFs(sfs.replace("nom",key),which=which))
                    if i == 0 : continue
                    syst.Scale(scalefactor_syst)
                    systList[key][0].Add(syst)
                del systList[key][1:]
            hist.Add(bhist)
            lists.append(hist)
            #blisthist.append(bhist)
        #print('totalBKG ',hist.Integral())
        NBins = hist.GetNbinsX()
        prevSigni = 0.0
        factor = 1.0 
        if mgo < 1400 : factor = 1.0
        bestBin = 0.0
        # getting the QCD number from antiSel SR which is name as CR1
        '''dHist = None
        for bkgf in bkg_files_CR1 : 
            bkg_name = bkgf.split('/')[-1].split('_')[1]
            bf = ROOT.TFile.Open(bkgf, "read")
            if bkg_name == 'Data' : dHist = bf.Get(bkg_name+'/sig_CR1_nom')
            #print(bkg_name+'/'+bkg+'/'+bkg_name+bkg+'sig_SR_nominal')
            bhist = bf.Get(bkg_name+'/sig_CR1_nom')
            for key in systListQCD : 
                systListQCD[key].append(bf.Get(bkg_name+'/sig_CR1_'+key))
            #print(bhist.Integral())
            which = ''
            scalefactor = 1.0
            if bkg_name == 'DiLepTT' : which = 'beta'
            elif bkg_name == 'SemiLepTT' :which = 'alpha'
            elif bkg_name = 'QCD' : which = 'delta'
            elif bkg_name != 'Data' : which = 'gamma'
            scalefactor = float(getSFs(sfs,which=which))
            #print(bkg,bkgf)
            #print(bkg)
            bhist.Scale(scalefactor)
            for key in systList : 
                for i,syst in enumerate(systList[key]) : 
                    scalefactor_syst = float(getSFs(sfs.replace("nom",key),which=which))
                    if i == 0 : continue
                    syst.Scale(scalefactor_syst)
                    systList[key][0].Add(syst)
                del systList[key][1:]
            hist.Add(bhist)
            lists.append(hist)
            #blisthist.append(bhist)
        #print('totalBKG ',hist.Integral())'''
        NBins = hist.GetNbinsX()
        prevSigni = 0.0
        factor = 1.0 
        if mgo < 1400 : factor = 1.0
        bestBin = 0.0
           
        '''for i in range(NBins,9000,-1):
            s = shist.Integral(i,NBins+1)/factor
            b = hist.Integral(i,NBins+1)

            #signi = simpleAsimov(s ,b)

            #if signi > prevSigni or b < 1.0: 
            if b < 0.35 : 
                #print('better significance when merged ', i ,' to ', NBins+1, ' oldsign =  ',prevSigni , ' newsign = ',  signi)
                bestBin = i
                #   prevSigni = signi
            else : 
                #print('better significance when merged ', i-1 ,' to ', NBins+1, ' oldsign =  ',prevSigni , ' newsign = ',  signi)
                bestBin = i-1
                break 
        # this is to inforce the best bin to be 96 which corresponding to DNN >= 0.95, this will ignor the significance calculations above

        lastbinW = NBins+1 - bestBin
        beforelastbin = bestBin -1 
        otherWs = []
        otherbins = []
        for i in range(0,5) : 
            #if beforelastbin < 9000 : continue 
            otherbinW = (2+i)*lastbinW
            bin = beforelastbin - otherbinW
            otherWs.append(otherbinW)
            otherbins.append([bin,beforelastbin])
            beforelastbin = bin -1 
        #otherbins.append([800,beforelastbin-1])
        otherbins = otherbins[::-1]
        #print(otherbins)
        SRbins = otherbins
        SRbins.append([bestBin,10001])'''
        bestBin = 9980
        nbins_full = 10000
        SRbins = [[9000,9196],[9197,9392],[9393,9588],[9589,9784],[9785,9979],[9980,10001]]
        #print(SRbins)
        if not oneBin : 
            for num, bin in enumerate(SRbins) : 
                bestBin = bin[0]
                NBins = bin[1]-1
                sigerr = ROOT.Double(0.)
                bkgerr = ROOT.Double(0.)
                Fullbkgerr = ROOT.Double(0.)
                Fullsigerr = ROOT.Double(0.)
                shist.Scale(1.0/factor)
                sigRate = shist.IntegralAndError(bestBin,NBins+1,sigerr)
                bkgRate = hist.IntegralAndError(bestBin,NBins+1,bkgerr)
                FullbkgRate = hist.IntegralAndError(9000,nbins_full+1,Fullbkgerr)
                FullsigRate = shist.IntegralAndError(9000,nbins_full+1,Fullsigerr)
                for key in systList : 
                    syst_int = systList[key][0].Integral(9000,nbins_full+1) #bestBin
                    systList[key].append(syst_int) 
                    impact = round((syst_int/(FullbkgRate+0.01)),2) #bkgRate
                    invimpact = round((1/(impact+0.001)),2)
                    systList[key].append(impact if impact >= 1.0 else invimpact )

                for key in SigsystList : 
                    syst_int = SigsystList[key][0].Integral(9000,nbins_full+1) #bestBin
                    SigsystList[key].append(syst_int)
                    impact = round((syst_int/((FullsigRate)+0.01))/factor,2) #sigRate
                    invimpact = round((1/(impact+0.001)),2)
                    SigsystList[key].append(impact if impact >= 1.0 else invimpact )
                #print(systList)
                #if sigRate == 0.0 : 
                #    print ('Signal',mgo,mlsp, 'has zero rate will not write the datacard for it ')
                #    continue 
                #sigerr  = 1.0 + sqrt(sigRate)/(sigRate+0.01)
                #bkgerr  = 1.0 + sqrt(bkgRate)/(bkgRate + 0.01)
                alphaE  = float(getSFs(sfs,which='alphaE')) + 1.0
                betaE  = float(getSFs(sfs,which='betaE')) + 1.0
                gammaE  = float(getSFs(sfs,which='gammaE')) + 1.0
                #del systList , SigsystList
                datacardsdir = os.path.join(outdir,'datacards/T1tttt_Scan_mGo'+str(int(mgo))+ '_mLSP'+str(int(mlsp)))
                if not os.path.exists(datacardsdir):
                    os.makedirs(datacardsdir)
                datacard = open(datacardsdir+'/bin_'+str(num)+'.txt', 'w'); 
                datacard.write("## Datacard for signal %s (with bins from  %s to %s)\n"%('T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp)),str(bestBin),str(NBins+1)))
                datacard.write("imax 1 number of bins \n")
                datacard.write("jmax * number of processes minus 1 \n")
                datacard.write("kmax * number of nuisance parameters \n")
                datacard.write(130*'-')
                datacard.write('\n')
                datacard.write("shapes *    ch1  FAKE \n");
                datacard.write(130*'-')
                datacard.write('\n')
                datacard.write("{:<12}{:<12}".format("bin","ch1")+'\n')
                datacard.write("{:<12}{:<12}".format("observation",round(hist.Integral(bestBin,NBins+1),2))+'\n')    
                datacard.write(130*'-')
                datacard.write('\n')
                datacard.write("{:<62}{:<30}{:<30}".format("bin","ch1","ch1")+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("process",'T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp)),"background")+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("process",'0',"1")+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("rate",round(sigRate,2),round(bkgRate,2))+'\n')
                datacard.write(130*'-')
                datacard.write('\n')
                datacard.write("{:<62}{:<30}{:<30}".format("MCstats"+'T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp))+' lnN',round(1.0+(sigerr/((sigRate)+0.01)),2),'-')+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("MCstatsbackground lnN",'-',round(1.0+(bkgerr/(bkgRate+0.01)),2))+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("lumi lnN",str(1.023 if int(args.year) == 2017 else 1.025),str(1.023 if int(args.year) == 2017 else 1.025))+'\n')
                datacard.write("{:<62}{:<30}{:<30}".format("DNNshape lnN",str(1.1),str(1.1))+'\n')
                for syst in systList : 
                    if "Down" in syst : continue
                    if syst in intersection_syst : 
                        if sigRate <= 0 :  sigVar = 1.0
                        else : sigVar = 1/2*(SigsystList[syst][2*(num+1)]+SigsystList[syst.replace("Up","Down")][2*(num+1)])
                        bkgVar = 1/2*(systList[syst][2*(num+1)]+systList[syst.replace("Up","Down")][2*(num+1)])
                        #print(syst,sigVar,bkgVar)
                        datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",str(round(sigVar,2)),str(round(bkgVar,2)))+'\n')
                    else :
                        bkgVar = 1/2*(systList[syst][2*(num+1)]+systList[syst.replace("Up","Down")][2*(num+1)])
                        datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",'-',str(round(bkgVar,2)))+'\n')
                for syst in SigsystList : 
                    if "Down" in syst : continue
                    if syst in intersection_syst : continue # already taken into account in the previous loop
                    else :
                        sigVar = 1/2*(SigsystList[syst][2*(num+1)]+SigsystList[syst.replace("Up","Down")][2*(num+1)])
                        datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",str(round(sigVar,2)),'-')+'\n')

        if oneBin :
                    
            bestBin = SRbins[-1][0]
            sigerr = ROOT.Double(0.)
            bkgerr = ROOT.Double(0.)
            Fullbkgerr = ROOT.Double(0.)
            Fullsigerr = ROOT.Double(0.)    
            shist.Scale(1.0/factor)
            sigRate = shist.IntegralAndError(bestBin,NBins+1,sigerr)
            bkgRate = hist.IntegralAndError(bestBin,NBins+1,bkgerr)
            FullbkgRate = hist.IntegralAndError(9000,nbins_full+1,Fullbkgerr)
            FullsigRate = shist.IntegralAndError(9000,nbins_full+1,Fullsigerr)
            
            for key in systList : 
                syst_int = systList[key][0].Integral(9000,nbins_full+1)
                systList[key].append(syst_int)
                impact = round((syst_int/(FullbkgRate+0.01)),2)
                invimpact = round((1/(impact+0.001)),2)
                systList[key].append(impact if impact >= 1.0 else invimpact )

            for key in SigsystList : 
                syst_int = SigsystList[key][0].Integral(9000,nbins_full+1)
                SigsystList[key].append(syst_int)
                impact = round((syst_int/((FullsigRate)+0.01))/factor,2)
                invimpact = round((1/(impact+0.001)),2)
                SigsystList[key].append(impact if impact >= 1.0 else invimpact )

            if sigRate == 0.0 : 
                print ('Signal',mgo,mlsp, 'has zero rate will not write the datacard for it ')
                continue 
            #sigerr  = 1.0 + sqrt(sigRate)/(sigRate+0.01)
            #bkgerr  = 1.0 + sqrt(bkgRate)/(bkgRate + 0.01)
            alphaE  = float(getSFs(sfs,which='alphaE')) + 1.0
            betaE  = float(getSFs(sfs,which='betaE')) + 1.0
            gammaE  = float(getSFs(sfs,which='gammaE')) + 1.0
            #del systList , SigsystList
            datacard = open(datacardsdir+ '/T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp))+'.txt', 'w'); 
            datacard.write("## Datacard for signal %s (with bins from  %s to %s)\n"%('T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp)),str(bestBin),str(NBins+1)))
            datacard.write("imax 1 number of bins \n")
            datacard.write("jmax * number of processes minus 1 \n")
            datacard.write("kmax * number of nuisance parameters \n")
            datacard.write(130*'-')
            datacard.write('\n')
            datacard.write("shapes *    ch1  FAKE \n");
            datacard.write(130*'-')
            datacard.write('\n')
            datacard.write("{:<12}{:<12}".format("bin","ch1")+'\n')
            datacard.write("{:<12}{:<12}".format("observation",round(hist.Integral(bestBin,NBins+1),2))+'\n')    
            datacard.write(130*'-')
            datacard.write('\n')
            datacard.write("{:<62}{:<30}{:<30}".format("bin","ch1","ch1")+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("process",'T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp)),"background")+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("process",'0',"1")+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("rate",round(sigRate,2),round(bkgRate,2))+'\n')
            datacard.write(130*'-')
            datacard.write('\n')
            datacard.write("{:<62}{:<30}{:<30}".format("MCstats"+'T1tttt_Scan_mGo' +str(int(mgo))+ '_mLSP'+str(int(mlsp))+' lnN',round(1.0+(sigerr/((sigRate)+0.01)),2),'-')+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("MCstatsbackground lnN",'-',round(1.0+(bkgerr/(bkgRate+0.01)),2))+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("lumi lnN",str(1.023 if int(args.year) == 2017 else 1.025),str(1.023 if int(args.year) == 2017 else 1.025))+'\n')
            datacard.write("{:<62}{:<30}{:<30}".format("DNNshape lnN",str(1.1),str(1.1))+'\n')
            for syst in systList : 
                if "Down" in syst : continue
                if syst in intersection_syst : 
                    sigVar = 1/2*(SigsystList[syst][2]+SigsystList[syst.replace("Up","Down")][2])
                    bkgVar = 1/2*(systList[syst][2]+systList[syst.replace("Up","Down")][2])
                    datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",str(round(sigVar,2)),str(round(bkgVar,2)))+'\n')
                else :
                    bkgVar = 1/2*(systList[syst][2]+systList[syst.replace("Up","Down")][2])
                    datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",'-',str(round(bkgVar,2)))+'\n')
            for syst in SigsystList : 
                if "Down" in syst : continue
                if syst in intersection_syst : continue # already taken into account in the previous loop
                else :
                    sigVar = 1/2*(SigsystList[syst][2]+SigsystList[syst.replace("Up","Down")][2])
                    datacard.write("{:<62}{:<30}{:<30}".format(syst.replace("_Up","")+args.year+" lnN",str(round(sigVar,2)),'-')+'\n')
            
            #datacard.write("{:<62}{:<30}{:<30}".format("bkguncert lnN",'-',1.1)+'\n')
            #datacard.write("{:<62}{:<30}{:<30}".format("alphauncert lnN",'-',alphaE)+'\n')
            #datacard.write("{:<62}{:<30}{:<30}".format("betauncert lnN",'-',betaE)+'\n')
            #datacard.write("{:<62}{:<30}{:<30}".format("gammauncert lnN",'-',gammaE)+'\n')

    if Limit and oneBin: 
        print ('limit will be calculated on HTC from ', datacardsdir)
        print ('in this step you will need a valid cmssw dir with higgscombne tools compiled, please check https://cms-analysis.github.io/HiggsAnalysis-CombinedLimit/ ')
        print ('you have scpecified CMSSW in this place, ', cmsswdir )
        
        import socket
        host = socket.gethostname()
        JDir = outdir
        if "lxplus" in host : 
            path = "/afs/cern.ch/work/a/amohamed/anaconda3/bin"
            anaconda = "/afs/cern.ch/work/a/amohamed/anaconda3/bin/activate"
            pyth = "/afs/cern.ch/work/a/amohamed/anaconda3/envs/hepML/bin/python"
        elif "desy.de" in host : 
            path = "/nfs/dust/cms/user/amohamed/anaconda3/bin"
            anaconda = "/nfs/dust/cms/user/amohamed/anaconda3/bin/activate"
            pyth = "/nfs/dust/cms/user/amohamed/anaconda3/envs/hepML/bin/python"

        logdir = outdir+'/Logs' 
        if not os.path.exists(logdir) :  os.makedirs(logdir)
        limitOutputdir = os.path.join(outdir,'datacards/limitOutput')
        if not os.path.exists(limitOutputdir) : os.makedirs(limitOutputdir)            
        file_list = glob.glob(datacardsdir+"/*.txt")
        for i,card in enumerate(file_list) : 
            cardname = card.split('/')[-1]
            cmd = '../combinedCards/'+cardname
            confDir = os.path.join(JDir,"job_"+str(i))
            if not os.path.exists(confDir) : 
                os.makedirs(confDir)
            print(cmd)
            exec = open(confDir+"/exec.sh","w+")
            exec.write("#"+"!"+"/bin/bash"+"\n")
            exec.write("source /etc/profile"+"\n")
            exec.write("source /cvmfs/cms.cern.ch/cmsset_default.sh"+"\n")
            exec.write("echo 'running job' >> "+os.path.abspath(confDir)+"/processing"+"\n")
            exec.write("workdir="+cmsswdir+"\n")
            exec.write("melalibdir=${CMSSW_BASE}/lib/slc6_amd64_gcc630/"+"\n")
            exec.write("exedir=`echo "+os.path.join(os.getcwd(),limitOutputdir)+"`"+"\n")
            exec.write("export LD_LIBRARY_PATH=${melalibdir}:$LD_LIBRARY_PATH"+"\n")
            exec.write("cd ${workdir}"+"\n")
            exec.write("eval `scramv1 runtime -sh`"+"\n")
            exec.write("cd ${exedir}"+"\n")
            exec.write("combine -M Asymptotic "+cmd+" -n "+cardname.replace('.txt','')+" "+"\n")
            exec.write("rm -rf "+os.path.abspath(confDir))
            exec.close()
        
        subFilename = os.path.join(JDir,"submitAlllimits.conf")
        subFile = open(subFilename,"w+")
        subFile.write("executable = $(DIR)/exec.sh"+"\n")
        subFile.write("universe =  vanilla")
        subFile.write("\n")
        subFile.write("should_transfer_files = YES")
        subFile.write("\n")
        subFile.write("log = "+"{}/job_$(Cluster)_$(Process).log".format(os.path.abspath(logdir)))
        subFile.write("\n")
        subFile.write("output = "+"{}/job_$(Cluster)_$(Process).out".format(os.path.abspath(logdir)))
        subFile.write("\n")
        subFile.write("error = "+"{}/job_$(Cluster)_$(Process).err".format(os.path.abspath(logdir)))
        subFile.write("\n")
        subFile.write("when_to_transfer_output   = ON_EXIT")
        #subFile.write("\n")
        #subFile.write('Requirements  = (OpSysAndVer == "SL6")')
        subFile.write("\n")
        subFile.write("queue DIR matching dirs "+JDir+"/job_*/")
        subFile.close()
        os.system("condor_submit "+subFilename)

