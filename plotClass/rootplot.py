import ROOT
import os

def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result

class rootplot(object) : 
    def __init__(self,input,outdir,All_files = {},outtext="sample") : 
        self.input      = input  
        self.outdir     = outdir 
        self.rootList   = find_all_matching('.root',self.input)
        self.group      = {}
        self.All_files  = All_files
        self.outtext=outtext
        if not os.path.exists(self.outdir) : os.makedirs(self.outdir)

        textsample = open(self.outdir+"/"+self.outtext+".txt", "w+")

        for g in self.All_files :
            textsample.write(g+" : "+"\n") 
            textsample.write("{:<20}{:<20}".format(" ",self.All_files[g]['scale'])+"\n")
            for f in self.All_files[g]['files']:
                    for rf in self.rootList:
                        if f in rf : 
                            textsample.write("{:<20}{:<20}".format(" ",rf)+"\n")
                            self.group.setdefault(g,[]).append(rf)
        
    
    def makeChain (self, filesList = [],Tname = "sf/t") : 
        ''' This function takes a group of files and add then to one Tchain
            to get one chain for each background kind'''
        chain = ROOT.TChain(Tname)
        for f in filesList : 
            chain.Add(f)
        return chain

    def makecuts(self,Tree = None,cutstring = "",extraCuts = "",name=None ) :
        ''' This Function is to apply selection on specific tree'''
        from ROOT import TCut
        #print(name)
        CUTtext = open(self.outdir+"/cuts.txt", "w+")
        CUTtext.write(cutstring+extraCuts+'\n')
        cutstring = TCut(cutstring+extraCuts)
        #print(cutstring)
        cut = cutstring.GetTitle()
        # when doing the cutflow, make a temp root file to avoid using a lot of memory
        if name : 
            fname = self.outdir+"/"+name+".root"
            #check if the file is already existing to avoid reproducing it
            if os.path.exists(fname) : 
                #open it and try to check if the tree is stored inside
                tt_out = self.makeChain( filesList = [fname],Tname = "t")
                if tt_out.LoadTree(0) >= 0 : 
                    return tt_out
                else: 
                    print(name,"the chain is not saved properly, will reproduce it")
                    tempFile = ROOT.TFile(self.outdir+"/"+name+".root","recreate")
                    #print(Tree.LoadTree(0))
                    tt_out = Tree.CopyTree(cut)
                    # need to use Write() here otherwise it will not be overwritten
                    tt_out.Write()
                    tempFile.Close()
                    # Stupied to write/Close the file to be able to load the chain correctly
                    tt_out = self.makeChain( filesList = [fname],Tname = "t")
                    return tt_out
            # if file is not existing, poduce it
            else : tempFile = ROOT.TFile(self.outdir+"/"+name+".root","recreate")
        # the normal workflow
        if Tree ==None : 
            print ('no Tree founded please cheack')
            pass  
        else : 
            tt_out = Tree.CopyTree(cut)
            return tt_out
        pass
    def chain_and_cut (self,filesList,Tname,cutstring ,extraCuts,name =None): 
        '''This function is to make TChain and apply selections at the same time'''
        chain = self.makeChain(filesList ,Tname)
        chain = self.makecuts(chain,cutstring,extraCuts,name)
        return chain