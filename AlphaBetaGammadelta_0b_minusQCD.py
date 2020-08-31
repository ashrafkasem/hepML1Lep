#!/usr/bin/env python
import sys
sys.argv.append( '-b-' )
import ROOT
from ROOT import RooFit as rf
from ROOT import std
ROOT.gROOT.SetBatch(True)
sys.argv.remove( '-b-' )

ROOT.gStyle.SetOptStat(0)
ROOT.gStyle.SetOptTitle(0)
from uncertainties import unumpy
import numpy as np 
import os 
import datetime

from math import hypot, sqrt, ceil,pi

import argparse
import htcondor
import scipy
import CMS_lumi

CMS_lumi.writeExtraText = 1
CMS_lumi.extraText = "Preliminary"
iPos = 0
if( iPos==0 ): CMS_lumi.relPosX = 0.1

## Canvas sizes
cwidth = 900
cheigth = 1200


b_list = ['DiLepTT','DY','QCD','SemiLepTT','SingleT','TTV','VV','WJ','Data']
others_bkg = ['DiLepTT','SemiLepTT','DY','SingleT','TTV','VV'] #keep the QCD out 'QCD',
CRs = ['SR','CR2','CR3','QCDCR']

def getVarFromHist(hist, varname, minx = None, maxx = None):
    if not minx : 
      minX = hist.GetXaxis().GetXmin();
    else :
      minX = minx
    if not maxx : 
      maxX = hist.GetXaxis().GetXmax();
    else : 
      maxX = maxx
    nBins = hist.GetNbinsX();

    # Create observable
    var = ROOT.RooRealVar(varname, varname, minX, maxX)
    var.setBins(nBins);

    return var

def doLegend():

    #leg = TLegend(0.63,0.525,0.87,0.875)
    leg = ROOT.TLegend(0.6,0.525,0.9,0.875)
    leg.SetBorderSize(1)
    leg.SetTextFont(62)
    leg.SetTextSize(0.03321678)
    leg.SetLineColor(1)
    leg.SetLineStyle(1)
    leg.SetLineWidth(1)
    leg.SetFillColor(0)
    leg.SetFillStyle(1001)
    return leg
#def LpTemplateFit(LpTemplates, prefix="", printDir='/afs/hephy.at/user/'+username[0]+'/'+username+'/www/pngCMG2/templateFit_Phys14V3/QCDestimation'):
def LpTemplateFit(LpTemplates, prefix="", printDir='./QCDestimation/templateFit', storePlot = False,USEdata=False, delta = 1.0,lumi= "20.1"):
    if not os.path.exists(printDir):
        os.makedirs(printDir) 
    #  cWJets = samples['W']
    #  cTTJets = samples['TT']
    #  cQCD = samples['QCD']
    histoEWK = LpTemplates['EWK']
    histoQCD = LpTemplates['QCD']
    histoDATA = LpTemplates['DATA']

    #Clone histograms from LpTemplates: EWK selected and QCD anti-selected
    template_EWK   = histoEWK.Clone()
    template_QCD   = histoQCD.Clone()#histoQCDsel.Clone()#histoDATAantiSel
    
    template_QCD.Scale(delta)

    print ("Nominal yields EWK:",template_EWK.Integral())#,'EWKets_PosPdg',template_EWKets_PosPdg.Integral(),'EWKets_NegPdg',template_EWKets_NegPdg.Integral()
    print ("Nominal yields QCD:",template_QCD.Integral())

    hData = histoDATA.Clone()
    #Observable
    x = getVarFromHist(histoDATA, "L_{p}")#,minx = 0, maxx = 1.5)
    #import the contents of 'data' ROOT histogram into a RooDataHist object 
    data=ROOT.RooDataHist("data","data",ROOT.RooArgList(x),hData)

    dh_EWK=ROOT.RooDataHist("mcEWK","mcEWK",ROOT.RooArgList(x),template_EWK)
    dh_QCD=ROOT.RooDataHist("mcQCD","mcQCD",ROOT.RooArgList(x),template_QCD)

    nevents = data.sumEntries()

    yield_EWK=ROOT.RooRealVar("EWK_yield","yieldEWK",nevents,0,2*nevents)#0.1,0,10**9)
    yield_QCD=ROOT.RooRealVar("QCD_yield","yieldQCD",nevents,0,2*nevents)#0.1,0,10**9)

    #MakePDFfromMChistograms
    model_EWK=ROOT.RooHistPdf("model_EWK","model_EWK",ROOT.RooArgSet(x),dh_EWK)
    model_QCD=ROOT.RooHistPdf("model_QCD","model_QCD",ROOT.RooArgSet(x),dh_QCD)
    model=ROOT.RooAddPdf("model","model",ROOT.RooArgList(model_EWK,model_QCD),ROOT.RooArgList(yield_EWK,yield_QCD))

    #  model=ROOT.RooAddPdf("model","model",ROOT.RooArgList(model_EWK),ROOT.RooArgList(yield_EWK))
    #CombinesmyMCsintoonePDFmodel

    #Plottheimportedhistogram(s)
    dframe=x.frame(rf.Title("Data"))
    data.plotOn(dframe)

    frame_EWK=x.frame(rf.Title("EWK"))
    model_EWK.plotOn(frame_EWK)
    frame_QCD=x.frame(rf.Title("QCD"))
    model_QCD.plotOn(frame_QCD)

    nllComponents = ROOT.RooArgList("nllComponents")
    nll=model.createNLL(data,rf.NumCPU(1))
    nllComponents.add(nll)

    sumNLL = ROOT.RooAddition("sumNLL","sumNLL", nllComponents)

    ROOT.RooMinuit(sumNLL).migrad()
    ROOT.RooMinuit(sumNLL).hesse()
    ROOT.RooMinuit(sumNLL).minos()#optional
    #model.fitTo(data,ROOT.RooFit.PrintLevel(-1))#,ROOT.RooFit.Extended(),ROOT.RooFit.SumW2Error(ROOT.kTRUE))
    #model.chi2FitTo(data,ROOT.RooLinkedList())#RooFit.SumW2Error(kTRUE))
    fitFrame=x.frame(rf.Title("Fit"))
    #model.paramOn(fitFrame,rf.Layout(0.68,0.9,0.75))
    data.plotOn(fitFrame,rf.LineColor(ROOT.kBlack),ROOT.RooFit.Name('Data'))#,ROOT.RooFit.Range(-0.1,1.5))
    model.plotOn(fitFrame,rf.LineColor(ROOT.kRed),ROOT.RooFit.Name('FullFit'))#,ROOT.RooFit.Range(-0.1,1.5))
    model.plotOn(fitFrame,rf.Components("model_QCD"),rf.LineColor(ROOT.kGreen),rf.LineStyle(ROOT.kDashed),ROOT.RooFit.Name('QCDfit'))
    model.plotOn(fitFrame,rf.Components("model_EWK"),rf.LineColor(ROOT.kBlue),rf.LineStyle(ROOT.kDashed),ROOT.RooFit.Name('EWKfit'))

    if storePlot: 
        canv=ROOT.TCanvas("canv","FitModel",cwidth,cheigth)


        ROOT.SetOwnership(canv, False)
        canv.SetTopMargin(canv.GetTopMargin()*1.5)
        topsize = 0.12*600./cheigth #if doRatio else 0.06*600./height
        canv.SetTopMargin(topsize)
        canv.cd()
        stackPad = ROOT.TPad("mainpad", "mainpad", 0, 0.40, 1, 1)
        ROOT.SetOwnership(stackPad, False)
        stackPad.SetBottomMargin(0.025)
        stackPad.SetTicks(1, 1)
        stackPad.Draw()

        ratioPad_1 = ROOT.TPad("ratiodpad_1", "ratiopad_1",0,0.20,1,0.40)
        ROOT.SetOwnership(ratioPad_1, False)
        ratioPad_1.SetTopMargin(0.001)
        ratioPad_1.SetBottomMargin(0.35)
        ratioPad_1.SetTicks(1,1)
        ratioPad_1.Draw()

        ratioPad_2 = ROOT.TPad("ratiodpad_2", "ratiopad_2",0,0,1,0.20)
        ROOT.SetOwnership(ratioPad_2, False)
        ratioPad_2.SetTopMargin(0.001)
        ratioPad_2.SetBottomMargin(0.35)
        ratioPad_2.SetTicks(1,1)
        ratioPad_2.Draw()

        stackPad.cd()

        #ROOT.gROOT.SetStyle("Plain")
        #ROOT.gROOT.SetStyle("Plain")#Removesgraybackgroundfromplots
        #ROOT.gPad.SetLeftMargin(0.15)
        fitFrame.GetYaxis().SetTitleOffset(1.4)
        fitFrame.GetXaxis().SetTitle('L_{P}')
        #c1.SetLeftMargin(c1.GetLeftMargin()*1.5)
        #c1.SetRightMargin(c1.GetRightMargin()*0.8)
        if cheigth == cwidth:
            fitFrame.GetYaxis().SetTitleOffset(2.2)
        fitFrame.SetMaximum(5000)
        fitFrame.GetXaxis().SetLabelOffset(999) ## send them away
        fitFrame.GetXaxis().SetTitleOffset(999) ## in outer space
        fitFrame.Draw()
        addHists = True
        if addHists:
            doTransp = False
            if doTransp:
                alpha = 0.7
                histoEWK.SetFillColorAlpha(histoEWK.GetFillColor(),alpha)
        
                histoQCD.SetFillColorAlpha(histoQCD.GetFillColor(),alpha)
            else:
                histoEWK.SetFillStyle(3001)
        
                histoQCD.SetFillStyle(3002)

            stack = ROOT.THStack('hs','hstack')
            stack.Add(histoQCD)
            stack.Add(histoEWK)
    
            stack.Draw("histsame,nostack")
            
            ROOT.SetOwnership( stack, 0 )
        fitFrame.Draw("same")
        leg = doLegend()
        ext = ""
        leg.AddEntry(fitFrame.findObject('Data'),'Data '+str(round(histoDATA.Integral(),2)),'lp')
        leg.AddEntry(fitFrame.findObject('FullFit'),'Full fit','l')
        leg.AddEntry(fitFrame.findObject('QCDfit'),'QCD fit '+str(round(yield_QCD.getVal(),2)),'l')
        leg.AddEntry(fitFrame.findObject('EWKfit'),'EWK fit '+str(round(yield_EWK.getVal(),2)),'l')

        leg.AddEntry("","#chi^{2}/ndof = "+str(round(fitFrame.chiSquare("FullFit", "Data",2),2)),'')
        
        if addHists:
            leg.AddEntry(0,"inputs:","")
            leg.AddEntry(histoQCD,'QCD '+str(round(histoQCD.Integral(),2)),'F')
            leg.AddEntry(histoEWK,'EWK '+str(round(histoEWK.Integral(),2)),'F')
    
            

        leg.Draw()  

        ROOT.SetOwnership( leg, 0 )
        #fitFrame.SetMaximum(fitFrame.GetMaximum()*2)
        #c1.SetLogy()
        CMS_lumi.lumi_13TeV = "%s fb^{-1}" % lumi
        CMS_lumi.CMS_lumi(canv, 4, iPos,0.09)

        ROOT.gPad.Update()
        ratioPad_1.cd()

        sumbkgs = ROOT.TH1F(histoEWK.Clone())

        sumbkgs.Add(histoQCD)

        pull = ROOT.TH1F(histoDATA.Clone())

        pull.Divide(sumbkgs)
        pull.SetMarkerStyle(20)
        pull.GetYaxis().SetTitle('Data/Pred.')
        #pull.GetXaxis().SetTitle(key)
        pull.GetYaxis().SetRangeUser(0.05, 1.95)
        pull.GetYaxis().SetDecimals(True)
        pull.SetLabelSize(0.14, "XY")
        pull.GetXaxis().SetTitleSize(.14)
        pull.GetYaxis().SetTitleSize(.14)
        pull.GetYaxis().SetLabelSize(0.11)
        pull.GetXaxis().SetLabelSize(0.11)
        pull.GetYaxis().SetTitleOffset(0.25)
        pull.GetYaxis().SetNdivisions(505)

        pull.Draw("EP")
        # Draw Line at ration == 1 
        line = ROOT.TLine(pull.GetXaxis().GetXmin(),1,pull.GetXaxis().GetXmax(),1)
        line.SetLineWidth(2);
        line.SetLineColor(58);
        line.Draw()
        
        ROOT.gPad.Update()
        canv.cd()
        ratioPad_2.cd()

        #print( ROOT.TH1F(fitFrame.findObject('FullFit')))
        
        template_EWK.Scale(yield_EWK.getVal()/template_EWK.Integral())

        template_QCD.Scale(yield_QCD.getVal()/template_QCD.Integral())
        fithist =  ROOT.TH1F(template_EWK.Clone())

        fithist.Add(template_QCD)

    #fitFrame.findObject('FullFit')
        pull_2 = ROOT.TH1F(histoDATA.Clone())

        pull_2.Divide(fithist)
        pull_2.SetMarkerStyle(20)
        pull_2.GetYaxis().SetTitle('Data/Fit.')
        #pull_2.GetXaxis().SetTitle(key)
        pull_2.GetYaxis().SetRangeUser(0.05, 1.95)
        pull_2.GetYaxis().SetDecimals(True)
        pull_2.SetLabelSize(0.14, "XY")
        pull_2.GetXaxis().SetTitleSize(.14)
        pull_2.GetYaxis().SetTitleSize(.14)
        pull_2.GetYaxis().SetLabelSize(0.11)
        pull_2.GetXaxis().SetLabelSize(0.11)
        pull_2.GetYaxis().SetTitleOffset(0.25)
        pull_2.GetYaxis().SetNdivisions(505)

        pull_2.Draw("EP")
        # Draw Line at ration == 1 
        line_2 = ROOT.TLine(pull_2.GetXaxis().GetXmin(),1,pull_2.GetXaxis().GetXmax(),1)
        line_2.SetLineWidth(2);
        line_2.SetLineColor(58);
        line_2.Draw()

        #canv.Print(printDir+'/'+prefix+'_TemplateFit.png')
        canv.Print(printDir+'/'+prefix+'_'+ext+'.pdf')
        canv.Print(printDir+'/'+prefix+'_'+ext+'.root')

        dataHist = fitFrame.getHist("Data")
        pullhist = dataHist.makePullHist(fitFrame.findObject('FullFit'), True)
        reslhist = dataHist.makeResidHist(fitFrame.findObject('FullFit'))

        canv_2=ROOT.TCanvas("canv_2","PullandRes",cwidth,cheigth)

        ROOT.SetOwnership(canv_2, False)
        canv_2.SetTopMargin(canv_2.GetTopMargin()*1.5)
        #topsize = 0.12*600./cheigth #if doRatio else 0.06*600./height
        #canv_2.SetTopMargin(topsize)
        canv_2.cd()
        PullPad = ROOT.TPad("PullPad", "PullPad", 0, 0.40, 1, 0.80)
        ROOT.SetOwnership(PullPad, False)
        PullPad.SetBottomMargin(0.025)
        PullPad.SetTicks(1, 1)
        PullPad.Draw()

        ResPad = ROOT.TPad("ResPad", "ResPad",0,0,1,0.40)
        ROOT.SetOwnership(ResPad, False)
        ResPad.SetTopMargin(0.001)
        ResPad.SetBottomMargin(0.35)
        ResPad.SetTicks(1,1)
        ResPad.Draw()
        PullPad.cd()
        pullhist.GetXaxis().SetTitle('L_{P}')
        pullhist.GetXaxis().SetTitleOffset(9999)
        pullhist.GetXaxis().SetLabelOffset(9999)
        pullhist.GetYaxis().SetTitle('pulls')
        pullhist.GetYaxis().SetTitleSize(0.06)
        pullhist.GetYaxis().SetLabelSize(0.06)
        pullhist.GetYaxis().SetTitleOffset(0.7)
        line_3 = ROOT.TLine(pullhist.GetXaxis().GetXmin(),0,pullhist.GetXaxis().GetXmax(),0)
        line_3.SetLineWidth(2);
        line_3.SetLineColor(58);
        pullhist.Draw()
        line_3.Draw()
        ResPad.cd()
        reslhist.GetXaxis().SetTitle('L_{P}')
        reslhist.GetXaxis().SetTitleSize(0.07)
        reslhist.GetXaxis().SetLabelSize(0.07)
        reslhist.GetYaxis().SetTitle('residuals')
        reslhist.GetYaxis().SetTitleSize(0.06)
        reslhist.GetYaxis().SetLabelSize(0.06)
        reslhist.GetYaxis().SetTitleOffset(0.7)
        line_4 = ROOT.TLine(reslhist.GetXaxis().GetXmin(),0,reslhist.GetXaxis().GetXmax(),0)
        line_4.SetLineWidth(2);
        line_4.SetLineColor(58);
        reslhist.Draw()
        line_4.Draw()
        canv_2.cd()
        CMS_lumi.CMS_lumi(canv_2, 4, iPos,0.09)
        ROOT.gPad.Update()
        canv_2.Print(printDir+'/'+prefix+'_PullandRes_'+ext+'.pdf')

        del canv

    del nllComponents

    res = {'EWK':{'template':template_EWK, 'yield':yield_EWK.getVal(), 'yield_high':yield_EWK.getVal()+yield_EWK.getErrorHi(), 'yield_low':yield_EWK.getVal()+yield_EWK.getErrorLo(), 
                        'yieldVar':(yield_EWK.getErrorHi()-yield_EWK.getErrorLo())**2},
            
           'QCD':{'template':template_QCD,'yield':yield_QCD.getVal(), 'yield_high':yield_QCD.getVal()+yield_QCD.getErrorHi(), 'yield_low':yield_QCD.getVal()+yield_QCD.getErrorLo(),
                        'yieldVar':(yield_QCD.getErrorHi()-yield_QCD.getErrorLo())**2},
            }

    del model, data, sumNLL
    return res

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
    lumi = '20.1'
    if "_16_" in infile.split("/")[-1]:
        lumi = "35.9"
    elif "_17_" in infile.split("/")[-1]:
        lumi = "41.9"
    elif "_18_" in infile.split("/")[-1]:
        lumi = "59.71"

    CMS_lumi.lumi_13TeV = "%s fb^{-1}" % lumi

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


    WJ_1_ = 0   ; WJ_1 = 0   ; WJ_2 = 0   ; WJ_3 = 0   ; WJ_QCDCR = 0
    TTJ_1_ = 0  ; TTJ_1 = 0  ; TTJ_2 = 0  ; TTJ_3 = 0  ; TTJ_QCDCR = 0
    data_1_ = 0 ; data_3 = 0 ; data_2 = 0 ; data_3 = 0 ; data_QCDCR = 0 
    QCD_1_ = 0  ; QCD_1 =0   ; QCD_2 = 0  ; QCD_3 = 0  ; QCD_QCDCR = 0 

    err   = ROOT.Double(0.) ; err_1 = ROOT.Double(0.) ; err_3 = ROOT.Double(0.)  ;err_4 = ROOT.Double(0.)  ;err_5 = ROOT.Double(0.)  ;err_6 = ROOT.Double(0.) 
    err_7 = ROOT.Double(0.) ; err_8 = ROOT.Double(0.) ; err_9 = ROOT.Double(0.)  ;err_10 = ROOT.Double(0.) ;err_11 = ROOT.Double(0.) ;err_12 = ROOT.Double(0.) 
    err_13 = ROOT.Double(0.) ; err_14 = ROOT.Double(0.) ; err_15 = ROOT.Double(0.)  ;err_16 = ROOT.Double(0.) ;err_17 = ROOT.Double(0.) ;err_18 = ROOT.Double(0.)
    err_19 = ROOT.Double(0.) ; err_20 = ROOT.Double(0.) ; err_21 = ROOT.Double(0.) ; err_22 = ROOT.Double(0.); err_23 = ROOT.Double(0.) ; err_24 = ROOT.Double(0.) 
    err_25 = ROOT.Double(0.) ; err_26 = ROOT.Double(0.) ; err_2 = ROOT.Double(0.)  ;err_0   = ROOT.Double(0.) ;err_00   = ROOT.Double(0.);err_000   = ROOT.Double(0.)
    bkg_list = b_list
    df_dict = {}
    dftex_dict = {}
    Lp_dict = {}
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
            if 'Lp' in hname and 'CR1' in hname : 
                if ('_CR1' in hname): Lp_dict[bkg] = hist
            elif not 'sig_0b' in hname: continue
            if 'Lp' in hname : continue 
            if bkg in others_bkg : 
                if ('_SR' in hname ): otherSR.Add(hist)
                if ('_CR1' in hname ): otherCR1.Add(hist)
                if ('_CR2' in hname ): otherCR2.Add(hist)
                if ('_CR3' in hname ): otherCR3.Add(hist)
                if ('_QCDCR' in hname ): otherQCDCR.Add(hist) 
            else : 
                if ('_SR' in hname and not 'Data' in bkg): 
                    txtSR.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                    if 'WJ' in bkg   : WJ_1  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_0), 2)
                    if 'QCD' in bkg  : QCD_1 = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_15), 2)
                elif '_CR1' in hname : 
                    txtCR1.write("{:<20}{:<20}{:<20}".format(hist.GetTitle(),round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err), 2), round(err, 2))+"\n")
                    if 'WJ' in bkg   : WJ_1_  = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_23), 2)
                    if 'QCD' in bkg  : QCD_1_ = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_25), 2)
                    if 'Data' in bkg : Data_1_ = round(hist.IntegralAndError(0, hist.GetNbinsX()+1, err_26), 2)

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
    TTJ_1_ = round(otherCR1.IntegralAndError(0, otherCR1.GetNbinsX()+1, err_20),2)
    txtCR2.write((60 *('='))+'\n')
    txtCR2.write("{:<20}{:<20}{:<20}".format(otherCR2.GetTitle(),round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_8), 2), round(err_8, 2))+"\n")
    TTJ_2 = round(otherCR2.IntegralAndError(0, otherCR2.GetNbinsX()+1, err_8),2)
    txtCR3.write((60 *('='))+'\n')
    txtCR3.write("{:<20}{:<20}{:<20}".format(otherCR3.GetTitle(),round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_9), 2), round(err_9, 2))+"\n")
    TTJ_3 = round(otherCR3.IntegralAndError(0, otherCR3.GetNbinsX()+1, err_9),2)

    txtQCDCR.write((60 *('='))+'\n')
    txtQCDCR.write("{:<20}{:<20}{:<20}".format(otherQCDCR.GetTitle(),round(otherQCDCR.IntegralAndError(0, otherQCDCR.GetNbinsX()+1, err_24), 2), round(err_24, 2))+"\n")
    TTJ_QCDCR = round(otherQCDCR.IntegralAndError(0, otherQCDCR.GetNbinsX()+1, err_24),2)
    print(otherQCDCR.Integral(),WJ_QCDCR)
    #del otherSR,otherCR1,otherCR2,otherCR3

    df_dict['Category'] = ["T5qqqqWW Cat.", "TTJets Cat.","WJets Cat."]
    df_dict['TTJets'] = [str(TTJ_1)+" +/- "+str(round(err_6,2)),str(TTJ_2)+" +/- "+str(round(err_8,2)),str(TTJ_3)+" +/- "+str(round(err_9,2))]
    df_dict['Others'] = [str(WJ_1)+" +/- "+str(round(err_0,2)),str(WJ_2)+" +/- "+str(round(err_1,2)),str(WJ_3)+" +/- "+str(round(err_4,2))]
    df_dict['QCD'] = [str(QCD_1)+" +/- "+str(round(err_15,2)),str(QCD_2)+" +/- "+str(round(err_16,2)),str(QCD_3)+" +/- "+str(round(err_17,2))]
    df_dict['Data - QCD_{MC}'] =   ["-- +/- --", str(Data_2-QCD_2)+" +/- "+str(round(sqrt(err_3**2 + err_16**2),2)), str(Data_3-QCD_3)+" +/- "+str(round(sqrt(err_5**2 +err_17**2),2))]

    dftex_dict['Category &'] = ["T5qqqqWW Cat. &", "TTJets Cat.&","WJets Cat.&"]
    dftex_dict['TTJets &'] = [str(TTJ_1)+" $\pm$ "+str(round(err_6,2))+" &",str(TTJ_2)+" $\pm$ "+str(round(err_8,2))+" &",str(TTJ_3)+" $\pm$ "+str(round(err_9,2))+" &"]
    dftex_dict['Others &'] = [str(WJ_1)+" $\pm$ "+str(round(err_0,2))+" &",str(WJ_2)+" $\pm$ "+str(round(err_1,2))+" &",str(WJ_3)+" $\pm$ "+str(round(err_4,2))+" &"]
    dftex_dict['QCD &'] = [str(QCD_1)+" $\pm$ "+str(round(err_15,2))+" &",str(QCD_2)+" $\pm$ "+str(round(err_16,2))+" &",str(QCD_3)+" $\pm$ "+str(round(err_17,2))+" &"]
    dftex_dict['Data $- QCD_{MC}$\\\\ \\hline'] =   ["-- $\pm$ --"+" \\\\ \\hline", str(Data_2-QCD_2)+" $\pm$ "+str(round(sqrt(err_3**2+err_16**2),2))+"  \\\\ \\hline", str(Data_3-QCD_3)+" $\pm$ "+str(round(sqrt(err_5**2 +err_17**2),2))+"  \\\\"]
    
    a = unumpy.umatrix([[TTJ_2,WJ_2],[TTJ_3,WJ_3]],[[err_8,err_1],[err_9,err_4]])
    
    b = unumpy.umatrix([[Data_2-QCD_2],[Data_3-QCD_3]],[[sqrt(err_3**2 + err_16**2)],[sqrt(err_5**2 +err_17**2)]])

    Y_ = a.I * b 
    Y = np.squeeze(np.asarray(unumpy.nominal_values(Y_)))
    Yerr = np.squeeze(np.asarray(unumpy.std_devs(Y_)))
    #print (np.allclose(np.dot(a, Y), b))
    alpha, beta , delta = round(Y[0],2) , round(Y[1],2), 1.0
    alphaerr, betaerr , deltaerr  = round(Yerr[0],2) , round(Yerr[1],2),0.0
    
    # make WJ histogram
    EWK = ROOT.TH1F('EWK','EWK',30, -0.5, 2.5)
    for key in Lp_dict: 
        if key in others_bkg:
            temp = ROOT.TH1F(Lp_dict[key].Clone())
            temp.Scale(alpha)
            EWK.Add(temp)
        elif key == 'WJ' : 
            temp = ROOT.TH1F(Lp_dict[key].Clone())
            temp.Scale(beta)
            EWK.Add(temp)

    LpTemplates = {}
    Lp_dict['EWK'] = EWK
    
    LpTemplates['EWK'] =  ROOT.TH1F(Lp_dict['EWK'].Clone())
    LpTemplates['QCD'] =  ROOT.TH1F(Lp_dict['QCD'].Clone())
    LpTemplates['DATA']=  ROOT.TH1F(Lp_dict['Data'].Clone())
    
    abgd_dict = {
                0 : {"alpha" :alpha , 
                     "beta"  : beta , 
                     "delta" : delta }
    }

    txtyields_ = open(mass_dir+'/yields_evolution.txt',"w+")
    if 1==1:#not "mucha" in infile: 
        for i in range(0,10):
            Fit = LpTemplateFit(LpTemplates, prefix="_"+str(i), printDir=mass_dir+'/QCDtemplateFit', storePlot = True,USEdata=False,delta = delta,lumi = lumi)
            QCD_norm = Fit['QCD']['yield']/(LpTemplates['QCD'].Integral())
            EWK_norm = Fit['EWK']['yield']/(LpTemplates['EWK'].Integral())
            print('QCD norm before the fit: ', LpTemplates['QCD'].Integral(), 'after the fit:' ,Fit['QCD']['yield'],'ratio', QCD_norm)
            print('EWK norm before the fit: ', LpTemplates['EWK'].Integral(), 'after the fit:' ,Fit['EWK']['yield'],'ratio', EWK_norm)
            txtyields_.write('QCD norm before the fit: '+ str(LpTemplates['QCD'].Integral())+ ' after the fit: ' +str(Fit['QCD']['yield'])+' ratio '+ str(QCD_norm)+"\n")
            txtyields_.write('EWK norm before the fit: '+ str(LpTemplates['EWK'].Integral())+ ' after the fit: ' +str(Fit['EWK']['yield'])+' ratio '+ str(EWK_norm)+"\n")

            a = unumpy.umatrix([[TTJ_2,WJ_2],[TTJ_3,WJ_3]],[[err_8,err_1],[err_9,err_4]])
        
            b = unumpy.umatrix([[Data_2-QCD_norm*QCD_2],[Data_3-QCD_norm*QCD_3]],[[sqrt(err_3**2 + err_16**2)],[sqrt(err_5**2 +err_17**2)]])

            Y_ = a.I * b 
            Y = np.squeeze(np.asarray(unumpy.nominal_values(Y_)))
            Yerr = np.squeeze(np.asarray(unumpy.std_devs(Y_)))

            alpha, beta , delta = round(Y[0],2) , round(Y[1],2), round(QCD_norm,2)
            alphaerr, betaerr , deltaerr  = round(Yerr[0],2) , round(Yerr[1],2),round(0.03,2)

            abgd_dict[i+1] = {}
            abgd_dict[i+1]['alpha'] = alpha;     abgd_dict[i+1]['beta'] = beta;    abgd_dict[i+1]["delta"] = delta
    

    import pandas as pd 
    df = pd.DataFrame(df_dict)
    dftex = pd.DataFrame(dftex_dict)
    
    txtalphabetagamma_ = mass_dir+'/alphabetagamma.txt'
    with open(txtalphabetagamma_, 'w') as f: df.to_string(f, col_space=10,index=False)
    txtalphabetagamma = open(txtalphabetagamma_, "a")
    txtalphabetagamma.write("\n")
    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ','alpha',' ','beta' ,' ', 'delta',' ')+"\n")
    if 1==1:#not "mucha" in infile: 
        txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ',abgd_dict[1]['alpha'], '+/-'+str(alphaerr),abgd_dict[1]['beta'] ,'+/-'+str(betaerr),abgd_dict[1]['delta'] ,'+/-'+str(deltaerr))+"\n")
        textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format(' ',abgd_dict[1]['alpha'], '+/-'+str(alphaerr),abgd_dict[1]['beta'] ,'+/-'+str(betaerr),abgd_dict[1]['delta'] ,'+/-'+str(deltaerr))+"\n")
    #else : 
    #    txtalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' ',abgd_dict[0]['alpha'], '+/-'+str(alphaerr),abgd_dict[0]['beta'] ,'+/-'+str(betaerr),"" ,"")+"\n")
    #    textalphabetagamma.write("{:<20}{:<12}{:<15}{:<12}{:<15}{:<12}{:<15}".format(' ',abgd_dict[0]['alpha'], '+/-'+str(alphaerr),abgd_dict[0]['beta'] ,'+/-'+str(betaerr),"" ,"")+"\n")

    latexalphabetagamma_ = mass_dir+'/alphabetagammaTable.tex'
    latexalphabetagamma = open(latexalphabetagamma_, "a+")
    latexalphabetagamma.write("\\documentclass[a4paper,11pt]{article} \n")
    latexalphabetagamma.write("\\usepackage{longtable} \n")
    latexalphabetagamma.write("\\begin{document} \n")
    latexalphabetagamma.write("\\tiny{\\begin{longtable}{|c || c | c | c|c| c|} \n")
    latexalphabetagamma.write("\\hline \n")
    dftex.to_string(latexalphabetagamma, col_space=10,index=False)
    latexalphabetagamma.write('\\hline\\hline \n')
    
    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format('   & ','$\\alpha$',' & ','$\\beta$' ,'&' ,'' , '& ')+"\\\\ \n")
    latexalphabetagamma.write('\\hline\\hline \n')
    if 1==1:#not "mucha" in infile:
        for i in range(0,10):
            latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' i = '+str(i)+' & ',str(abgd_dict[i]['alpha']), ' $\pm$ '+str(alphaerr)+' & ',str(abgd_dict[i]['beta']) ,' $\pm$ '+str(betaerr)+' & ',str(abgd_dict[i]['delta']) ,''+' & \\\\' )+"\n")
            latexalphabetagamma.write('\\hline\\hline\\hline \n')
    #else : 
    #    latexalphabetagamma.write("{:<20}{:<12}{:<10}{:<12}{:<10}{:<12}{:<10}".format(' & ',str(abgd_dict[0]['alpha']), ' $\pm$ '+str(alphaerr)+' & ',str(abgd_dict[0]['beta']) ,' $\pm$ '+str(betaerr)+' & ',"" ,' & \\\\' )+"\n")
    #    latexalphabetagamma.write('\\hline\\hline\\hline \n')
    latexalphabetagamma.write("\\caption{ Normalizaton parameters for backgrounds, calculated from control categories.}\n")
    latexalphabetagamma.write("\\label{normalize}\n")
    latexalphabetagamma.write("\\end{longtable}}\n")
    latexalphabetagamma.write("\\end{document} \n")