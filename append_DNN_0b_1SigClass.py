#!/usr/bin/env python
import argparse
import os

def load_model(pathToModel):
    '''Load a previously saved model (in h5 format)'''
    #print (" Loading the model from ",pathToModel)
    # load json and create model
    json_file = open(pathToModel+'.json', 'r')
    loaded_model_json = json_file.read()
    json_file.close()
    model = model_from_json(loaded_model_json)
    # load weights into new model
    model.load_weights(pathToModel+'.h5')
    return model


def Predict_Keras(infile,outdir,var_list,class_list,model = None) : 

    model = load_model(model)
    model.compile(loss='sparse_categorical_crossentropy',metrics=['accuracy'],optimizer='adam')
    
    it = uproot.open(infile)["sf/t"]
    file_in = ROOT.TFile(infile, "READ")
    tree_in = file_in.Get("sf/t")
    file_out = ROOT.TFile(os.path.join(outdir,infile.split("/")[-1]), "RECREATE")
    file_out.mkdir('sf')
    file_out.cd('sf')
    #tree_out = file_out.Get("sf/t")
    TT1l_val = array.array('f', [0.])
    WJet_val = array.array('f', [0.])
    Sig_val = array.array('f', [0.])


    if ('T1tttt' in infile) : 
        return
    else :
        p_df = it.pandas.df(var_list+['Event','Run','Lumi','nLep','Selected','nVeto'])
        p_df = p_df.loc[(p_df['nLep'] == 1) & (p_df['nJets30Clean'] >= 3)& (p_df['nVeto'] == 0)& (p_df['HT'] > 500)& (p_df['LT'] > 250)]
        p_df = p_df.reset_index(drop=True)
        prediction = pd.DataFrame(model.predict_proba(p_df[var_list].values),columns=['TTJ', 'WJ','Sig'])

    tree_out = tree_in.CopyTree("(nLep == 1) && (nJets30Clean >= 3)&& (nVeto == 0)&& (HT > 500)&& (LT > 250)")

    TT1l  = tree_out.Branch('TTJ_0b', TT1l_val, 'TTJ_0b/F')
    WJet  = tree_out.Branch('WJ_0b', WJet_val, 'WJ_0b/F')
    Sig  = tree_out.Branch('sig_0b', Sig_val, 'sig_0b/F')
    
    prediction['Event'] = p_df['Event']
    prediction['Run'] =  p_df['Run']
    prediction['Lumi'] =  p_df['Lumi']
    #print(prediction)

    t_start = time.time()

    #print prediction.size
    
    for i_ev in range(tree_out.GetEntries()):
        #print i_ev
        if i_ev % 10000 == 0:
            print ('Event', i_ev,"/",tree_out.GetEntries())
            t_current = time.time()
            print ('Time', t_current - t_start)
            t_start = t_current
        #tree_out.GetEntry(i_ev)
        tree_out.GetEntry(i_ev)
        #if not ((tree_out.nLep == 1) & (tree_out.nJets30Clean >= 3)& (tree_out.Selected == 1)& (tree_out.nVeto == 0)& (tree_out.HT > 500)& (tree_out.LT > 250)) : continue
        """ the selections cuts in pandas seeme to keep all the events order similer to Tcut in root,
        to save time i keep the ordering with event numbers which i checked to be 100% similer to the order in pandas
        you can check that by un commenting the next couple of lines and you will see"""
        #df_idx = prediction.loc[(prediction['Event'] == tree_out.Event )&(prediction['Run'] == tree_out.Run)&(prediction['Lumi'] == tree_out.Lumi)].index
        #print(df_idx, i_ev)

        TT1l_val[0] = prediction['TTJ'][i_ev]
        WJet_val[0] = prediction['WJ'][i_ev]
        Sig_val[0] = prediction['Sig'][i_ev]

        TT1l.Fill()
        WJet.Fill()
        Sig.Fill()

    #file_out.Delete("t;1")
    #file_out.Delete("sf/t;1")
    tree_out.Write("t", ROOT.TObject.kOverwrite)
    file_out.Close()

def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Evaluate DNN with scikit-learn')
    parser.add_argument('--indir', help='List of datasets to process',default=None, metavar='indir')
    parser.add_argument('--outdir', help='output directory',default=None, metavar='outdir')
    parser.add_argument('--infile', help='infile to process',default=None, metavar='infile')
    parser.add_argument('--exec', help="wight directory", default='./batch/append_exec_0b1Sig.sh', metavar='exec')
    parser.add_argument('--batchMode','-b', help='Batch mode.',default=False, action='store_true')
    parser.add_argument('--model', help='model to be used',default=None, metavar='model')

    args = parser.parse_args()

    var_list = ['MET', 'MT', 'Jet2_pt','Jet1_pt', 'Lep_pt', 'LT', 'HT','nTop_Total_Combined', 'nJets30Clean', 'dPhi',"Lep_relIso","Lep_miniIso","iso_pt","iso_MT2","nWLoose","nWMedium","nWTight"]
    
    wdir = os.getcwd()
    
    if not args.batchMode and args.infile: 
        import ROOT
        import time
        import array
        import operator
        import uproot
        os.environ["KERAS_BACKEND"] = "tensorflow"
        import keras
        from keras.models import model_from_json
        import pandas as pd
        
        outdir = args.outdir
        if not os.path.exists(outdir):
            os.makedirs(outdir) 

        if ( "T1tttt" in args.infile or 'T5qqqq' in args.infile) :
            Predict_Keras(args.infile,outdir,var_list,['TTJ_0b','WJ_0b','sig_0b'],model = args.model)
        else : 
            Predict_Keras(args.infile,outdir,var_list,['TTJ_0b','WJ_0b','sig_0b'],model = args.model)
    else : 
        
        outdir = os.path.join(args.outdir,os.path.basename(os.path.normpath(args.indir)))
        if not os.path.exists(outdir):
            os.makedirs(outdir) 

        logdir = outdir+'/Logs' 
        if not os.path.exists(logdir):
            os.makedirs(logdir) 

        Filenamelist = find_all_matching(".root",args.indir) 
        #print (Filenamelist)
        import socket
        host = socket.gethostname()

        if "lxplus" in host : 
            path = "/afs/cern.ch/work/a/amohamed/anaconda3/bin"
            anaconda = "/afs/cern.ch/work/a/amohamed/anaconda3/bin/activate"
            pyth = "/afs/cern.ch/work/a/amohamed/anaconda3/envs/hepML/bin/python"
        elif "desy.de" in host : 
            path = "/nfs/dust/cms/user/amohamed/anaconda3/bin"
            anaconda = "/nfs/dust/cms/user/amohamed/anaconda3/bin/activate"
            pyth = "/nfs/dust/cms/user/amohamed/anaconda3/envs/hepML/bin/python"
        
        for i,fc in enumerate(Filenamelist) : 
            #if not "Run2018" in fc : continue
            confDir = os.path.join(outdir,"job_"+str(i))
            if not os.path.exists(confDir) : 
                os.makedirs(confDir)
            exec = open(confDir+"/exec.sh","w+")
            exec.write("#"+"!"+"/bin/bash"+"\n")
            exec.write("eval "+'"'+"export PATH='"+path+":$PATH'"+'"'+"\n")
            exec.write("source "+anaconda+" hepML"+"\n")
            exec.write("cd "+wdir+"\n")
            exec.write("echo 'running job' >> "+confDir+"/processing"+"\n")
            exec.write("echo "+wdir+"\n")
            exec.write(pyth+" append_DNN_0b_1SigClass.py --infile "+fc+"  --model "+args.model+" --outdir "+outdir+" --indir "+args.indir)
            exec.write("\n")
            # let the script deletes itself after finishing the job
            exec.write("rm -rf "+confDir)

            exec.close()
        
        subFilename = os.path.join(outdir,"submitAllgrid.conf")
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
        subFile.write("\n")
        subFile.write('Requirements  = ( OpSysAndVer == "CentOS7")')
        subFile.write("\n")
        subFile.write("queue DIR matching dirs "+outdir+"/job_*/")
        if "lxplus" in host : 
            subFile.write("\n")
            subFile.write('+JobFlavour = "longlunch"')
        subFile.close()
        os.system("condor_submit "+subFilename)
        
