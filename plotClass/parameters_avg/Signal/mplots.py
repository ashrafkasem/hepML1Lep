varList.append(["TTDi","TTDi","DNN classifier t#bar{t} ll",[20,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig","sig","DNN classifier T1t^{4}",[20,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["WJ","WJ","DNN classifier W+jets",[20,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["TTS","TTS","DNN classifier t#bar{t} l",[20,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])

varList.append(["sig_100bins","sig","DNN classifier T1t^{4}",[100,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_10bins","sig","DNN classifier T1t^{4}",[10,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_20bins","sig","DNN classifier T1t^{4}",[20,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_30bins","sig","DNN classifier T1t^{4}",[30,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_40bins","sig","DNN classifier T1t^{4}",[40,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_50bins","sig","DNN classifier T1t^{4}",[50,0,1],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["sig_VR","sig","DNN classifier T1t^{4}",[3,0,1],"LogY",["MoreY",1000],'IncludeOverflows',['varbin',[0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.5,0.6,0.7,0.8,0.9,0.98,0.99,0.995,1.0],True]])
varList.append(["sig_VR2","sig","DNN classifier T1t^{4}",[3,0,1],"LogY",["MoreY",1000],'IncludeOverflows',['varbin',[0.,0.1,0.2,0.3,0.4,0.5,0.6,0.7,0.8,0.9,0.98,0.99,0.995,1.0],True]])
varList.append(["sig_VR3","sig","DNN classifier T1t^{4}",[3,0,1],"LogY",["MoreY",1000],'IncludeOverflows',['varbin',[0.98,0.99,0.998,1.0],False]])

varList.append(["CatTT1Lep","(TTS >TTDi ) && (TTS >sig) && (TTS >WJ)","t#bar{t} l + jets Event Category",[2,0.0,2.0],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["CatTT2Lep","(TTDi > TTS) && (TTDi >sig) && (TTDi >WJ )","t#bar{t} ll + jets Event Category",[2,0.0,2.0],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["CatWJ","(WJ >TTDi ) && (WJ >sig)&& (WJ  > TTS )","W+jets Event Category",[2,0.0,2.0],"LogY",["MoreY",1000],'IncludeOverflows'])
varList.append(["CatSig","(sig >TTDi ) && (sig >TTS) && (sig >WJ)","T1t^{4} Event Category",[2,0.0,2.0],"LogY",["MoreY",1000],'IncludeOverflows'])
