#!/usr/bin/env python
import sys,os
#sys.path.append("/nfs/dust/cms/user/amohamed/anaconda3/envs/hepML/lib/python3.6/site-packages/")
import argparse
import glob

import time 
def find_all_matching(substring, path):
    result = []
    for root, dirs, files in os.walk(path):
        for thisfile in files:
            if thisfile.startswith("."): continue 
            if substring in thisfile:
                result.append(os.path.join(root, thisfile ))
    return result
    
def getSFs(sffile,mass='1900_100',which='alpha') : 
    token = open(sffile,'r')
    linestoken=token.readlines()
    if which == 'alpha' :  tokens_column_number = 1
    if which == 'alphaE' :  tokens_column_number = 2
    if which == 'beta' :  tokens_column_number = 3
    if which == 'betaE' :  tokens_column_number = 4
    if which == 'gamma' :  tokens_column_number = 5
    if which == 'gammaE' :  tokens_column_number = 6
    resulttoken=[]
    masstoken=[]
    for x in linestoken:
        resulttoken.append(x.split()[tokens_column_number])
        masstoken.append(x.split()[0])
    token.close()
    idx = masstoken.index(mass)
    return resulttoken[idx]
    
def getSFs_(sffile,which='alpha') : 
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

def getSFs_0b(sffile,mass='1900_100',which='alpha') : 
    token = open(sffile,'r')
    linestoken=token.readlines()
    if which == 'alpha' :  tokens_column_number = 1
    if which == 'alphaE' :  tokens_column_number = 2
    if which == 'beta' :  tokens_column_number = 3
    if which == 'betaE' :  tokens_column_number = 4
    resulttoken=[]
    masstoken=[]
    for x in linestoken:
        resulttoken.append(x.split()[tokens_column_number])
        masstoken.append(x.split()[0])
    token.close()
    idx = masstoken.index(mass)
    return resulttoken[idx]
    
def getSFs_0b_(sffile,which='alpha') : 
    token = open(sffile,'r')
    linestoken=token.readlines()
    if which == 'alpha' :  tokens_column_number = 0
    if which == 'alphaE' :  tokens_column_number = 1
    if which == 'beta' :  tokens_column_number = 2
    if which == 'betaE' :  tokens_column_number = 3
    if which == 'delta' :  tokens_column_number = 4
    if which == 'deltaE' :  tokens_column_number = 5
    resulttoken=[]
    #masstoken=[]
    for x in linestoken:
        resulttoken.append(x.split()[tokens_column_number])
        #masstoken.append(x.split()[0])
    token.close()
    #idx = masstoken.index(mass)
    return resulttoken[1]

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Runs a NAF batch system for 1l plotter', formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('--indir', help='input root files', metavar='indir')
    parser.add_argument('--lumi','-L' ,help="", default="35.9", metavar='lumi')
    parser.add_argument('--outdir' ,help="outputdir", default=None, metavar='outdir')
    parser.add_argument('--scale' ,help="if YES scale MC to Data, if not use alpha/beta/gamma", default=False, action='store_true')
    parser.add_argument('--param','--parameter' ,help="parameters directory", default=None, metavar='param')
    parser.add_argument('--abg' ,help="text file that has the alpha/beta/gamma values to be used", default=None, metavar='abg')
    parser.add_argument('--only' ,help="choose a cut text to run with it only", default=None, metavar='only')
    parser.add_argument('--blind' ,help="to blind the inclusive distributions",default=False,action='store_true')
    parser.add_argument('--blindall' ,help="to blind the all distributions",default=False,action='store_true')
    parser.add_argument("--showSF", default=False, help="show the SF or not",action='store_true')    
    parser.add_argument("--showCount", default=False, help="show the counts in legend",action='store_true')  
    parser.add_argument("--HEM", default=False, help="apply the HEM Fix, only for 2018",action='store_true')    

    args = parser.parse_args()
    HEMtext = ""
    if args.lumi == "35.9" : 
        year = "2016"
    elif args.lumi == "41.9" : 
        year = "2017"
    elif args.lumi == "59.74" : 
        year = "2018"
        if args.HEM : HEMtext = " --HEM "
    elif args.lumi == "137.54" : 
        year = "20161718"
    else : print("lumi must be in [35.9,41.9,59.74]") ; sys.exit() 
    
    cmd = " --indir "+args.indir+" --lumi "+args.lumi+ " --YmaX 0.0  --YmiN 0.1 --rmax 1.95 --rmin 0.05 --doRatio --year "+year + HEMtext

    if args.showSF : 
        cmd += " --showSF "

    if args.showCount : 
        cmd += " --showCount "

    if not "_0b" in args.param : 
        cmd+= ' --mb '

    cmd+=" --cuts "+args.param+"/inclusive.txt --varList "+args.param+"/baseplots.py"
    if args.scale : 
        cmd+= " --scale_bkgd_toData "
    else : 
        cmd+=" --do_alphabetagamma --scale_bkgd_toData "
        if args.abg == None: 
            print("you've decided to use the alpha/beta/gamm to scale background to data, but you haven't used --abg to provide the text file that contain them")
            sys.exit()
    if args.blindall : 
        cmd+= " --blind "
    list_masses = [x[0] for x in os.walk(args.param)]
    list_masses_ = list_masses[1:]
    cmd_ = cmd
    cmd_array = []
    for mdir in list_masses_ : 
        #to retrive the common options in the cmd
        cmd = cmd_
        if args.only != None : 
            myFiles = glob.glob(mdir+"/*"+args.only+"*.txt")
        else : myFiles = glob.glob(mdir+"/*.txt")
        if "_" in str(mdir.split("/")[-1]):
            mgo = str(mdir.split("/")[-1]).split("_")[0]
            mlsp = str(mdir.split("/")[-1]).split("_")[1]
            cmd+=" --Smass "+mgo+"_"+mlsp+" --outdir "+os.path.join(args.outdir,mdir.split("/")[-1])
            if not args.scale : 
                if not "_0b" in args.param : 
                    alpha = (getSFs(args.abg,mass=mdir.split("/")[-1],which="alpha"))
                    beta = (getSFs(args.abg,mass=mdir.split("/")[-1],which="beta"))
                    gamma = (getSFs(args.abg,mass=mdir.split("/")[-1],which="gamma"))
                    cmd+= " --alpha "+alpha+" --beta "+beta+" --gamma "+gamma
                else : 
                    alpha = (getSFs_0b(args.abg,mass=mdir.split("/")[-1],which="alpha"))
                    beta = (getSFs_0b(args.abg,mass=mdir.split("/")[-1],which="beta"))
                    cmd+= " --alpha "+alpha+" --beta "+beta
        else :
            cmd+=" --outdir "+os.path.join(args.outdir,mdir.split("/")[-1]) + " --Smass 1500_1000  2300_100  1800_1300  2300_800 "
            if not args.scale : 
                if not "_0b" in args.param :     
                    alpha = (getSFs_(args.abg,which="alpha"))
                    beta = (getSFs_(args.abg,which="beta"))
                    gamma = (getSFs_(args.abg,which="gamma"))
                    delta = (getSFs_(args.abg,which="delta"))
                    cmd+= " --alpha "+alpha+" --beta "+beta+" --gamma "+gamma+" --delta "+delta
                else : 
                    alpha = (getSFs_0b_(args.abg,which="alpha"))
                    beta = (getSFs_0b_(args.abg,which="beta"))
                    delta = (getSFs_0b_(args.abg,which="delta"))
                    cmd+= " --alpha "+alpha+" --beta "+beta+" --delta "+delta

        if args.blind and os.path.exists(mdir+"/mplots_blind.py"): 
            incl_cmd = cmd+" --mvarList "+mdir+"/mplots_blind.py "
        elif not args.blind and os.path.exists(mdir+"/mplots.py"): 
            incl_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
        else : incl_cmd = cmd
        if args.only == None : 
            if os.path.exists(args.param+"/0bCS.txt") : 
                zerob_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                zerob_cmd = zerob_cmd.replace("inclusive.txt","0bCS.txt")
                cmd_array.append(zerob_cmd)
            if os.path.exists(args.param+"/2LCS.txt") : 
                diLepCS_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                diLepCS_cmd = diLepCS_cmd.replace("inclusive.txt","2LCS.txt")
                cmd_array.append(diLepCS_cmd)
            if os.path.exists(args.param+"/QCD_ech_sel.txt") : 
                QCD_e_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCD_e_cmd = QCD_e_cmd.replace("inclusive.txt","QCD_ech_sel.txt")
                cmd_array.append(QCD_e_cmd)
            if os.path.exists(args.param+"/QCD_mch_sel.txt") : 
                QCD_m_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCD_m_cmd = QCD_m_cmd.replace("inclusive.txt","QCD_mch_sel.txt")
                cmd_array.append(QCD_m_cmd)
            if os.path.exists(args.param+"/QCD_ech_antisel.txt") : 
                QCDA_e_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCDA_e_cmd = QCDA_e_cmd.replace("inclusive.txt","QCD_ech_antisel.txt")
                cmd_array.append(QCDA_e_cmd)
            if os.path.exists(args.param+"/QCD_mch_antisel.txt") : 
                QCDA_m_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCDA_m_cmd = QCDA_m_cmd.replace("inclusive.txt","QCD_mch_antisel.txt")
                cmd_array.append(QCDA_m_cmd)
            if os.path.exists(args.param+"/inclusive_ech_antisel.txt") : 
                incA_e_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                incA_e_cmd = incA_e_cmd.replace("inclusive.txt","inclusive_ech_antisel.txt")
                cmd_array.append(incA_e_cmd)
            if os.path.exists(args.param+"/inclusive_mch_antisel.txt") : 
                incA_m_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                incA_m_cmd = incA_m_cmd.replace("inclusive.txt","inclusive_mch_antisel.txt")
                cmd_array.append(incA_m_cmd)
            if os.path.exists(args.param+"/QCD_antisel.txt") : 
                QCD_AntiAll_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCD_AntiAll_cmd = QCD_AntiAll_cmd.replace("inclusive.txt","QCD_antisel.txt")
                cmd_array.append(QCD_AntiAll_cmd)
            if os.path.exists(args.param+"/QCD_sel.txt") : 
                QCD_all_cmd = cmd+" --mvarList "+mdir+"/mplots.py "
                QCD_all_cmd = QCD_all_cmd.replace("inclusive.txt","QCD_sel.txt")
                cmd_array.append(QCD_all_cmd)
            cmd_array.append(incl_cmd)
        
        for mcut in myFiles : 
            othercmd = cmd
            #print(mcut)
            if ("_postHEM" in str(mcut) or "_preHEM" in str(mcut)) and not year == "2018" : continue
            if ("QCD_mchsel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","QCD_mch_sel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("QCD_echsel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","QCD_ech_sel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("QCD_mchantisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","QCD_mch_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("QCD_echantisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","QCD_ech_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("Sig_SR_ech_antisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","inclusive_ech_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("Sig_VR_ech_antisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","inclusive_ech_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("Sig_SR_mch_antisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","inclusive_mch_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("Sig_VR_mch_antisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","inclusive_mch_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif ("antisel" in str(mcut)) : 
                othercmd = othercmd.replace("inclusive.txt","inclusive_antisel.txt")
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif "inclusive_njseq6" in str(mcut) and os.path.exists(mdir+"/mplots.py") :
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "
            elif "inclusive" in str(mcut) and os.path.exists(mdir+"/mplots_blind.py") :
                othercmd = othercmd+" --mvarList "+mdir+"/mplots_blind.py "
            elif not "inclusive" in str(mcut) and os.path.exists(mdir+"/mplots.py") : 
                othercmd = othercmd+" --mvarList "+mdir+"/mplots.py "


            othercmd+= " --mcuts "+mcut
            if ('Sig.txt' in mcut or 'Sig_ech.txt' in mcut or 'Sig_mch.txt' in mcut or "Sig_lastbin" in mcut or "Sig_ge" in mcut or 'Sig_nj7.txt' in mcut or 'Sig_nj6.txt' in mcut or 'Sig_SR.txt' in mcut or 'Sig_SR_ech.txt' in mcut or 'Sig_SR_mch.txt' in mcut or 'Sig_bin' in mcut ) and (not "Anti" in mcut or "QCD" in mcut) : othercmd +=' --blind '
            cmd_array.append(othercmd)
            
    cmd_array = [x.replace("//","/") for x in cmd_array]
    
    logdir = os.path.join(args.outdir,'Logs')
    wdir = os.getcwd()
    if not os.path.exists(logdir):
        os.makedirs(logdir) 
    outtext = open(args.outdir+"/commands.txt", "w+")
    for cmd in cmd_array : 
        outtext.write(cmd+'\n')
    
    JDir = args.outdir
    
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

    for i,comd in enumerate(cmd_array) : 
        confDir = os.path.join(JDir,"job_"+str(i))
        if not os.path.exists(confDir) : 
            os.makedirs(confDir)
        print(comd)
        exec = open(confDir+"/exec.sh","w+")
        exec.write("#"+"!"+"/bin/bash"+"\n")
        exec.write("eval "+'"'+"export PATH='"+path+":$PATH'"+'"'+"\n")
        exec.write("source "+anaconda+" hepML"+"\n")
        exec.write("cd "+wdir+"\n")
        exec.write("echo 'running job' >> "+confDir+"/processing"+"\n")
        exec.write("echo "+wdir+"\n")
        exec.write(pyth+" RoPlotter.py "+comd)
        exec.write("\n")
        # let the script deletes itself after finishing the job
        exec.write("rm -rf "+confDir)
        exec.close()

    subFilename = os.path.join(JDir,"submitAllplots.conf")
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
    subFile.write('Requirements  = ( OpSysAndVer == "CentOS7" || OpSysAndVer == "SL6")')
    subFile.write("\n")
    subFile.write("queue DIR matching dirs "+JDir+"/job_*/")
    subFile.close()
    os.system("condor_submit -name s13 "+subFilename)
