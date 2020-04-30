varList.append(["HT", "HT", "H_{T} [GeV]", [42, 0, 3000], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["LT", "LT", "L_{T} [GeV]", [35, 0, 1200] , "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["MET", "MET", "MET [GeV]", [35, 0, 1200] , "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["Lep_pt", "Lep_pt", "Lep p_{T} [GeV]", [32, 0, 800], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["nJets30Clean", "nJets30Clean", "jet multiplicity", [15, 0, 15] , "LogY",["MoreY",2000]])
varList.append(["nBJet","nBJet","b-jet multiplicity (Med)",[8, 0, 8], "LogY",["MoreY",2000]])
varList.append(["Jet2_pt","Jet2_pt","2. Leading jet p_{T}",[35, 0, 1000], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["Jet1_pt","Jet1_pt","1. Leading jet p_{T}",[35, 0, 1000], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["nTop_Total_Combined","nTop_Total_Combined","Top multiplicity", [5, 0, 5] , "LogY",["MoreY",2000]])
varList.append(["dPhi", "fabs(dPhi)", "#Delta #varphi (lep,W)" ,[30, 0, 3.142],"LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["dPhi_1bins", "fabs(dPhi)", "#Delta #varphi (lep,W)" ,[1, 0, 3.142],"LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["dPhi_2bins", "fabs(dPhi)", "#Delta #varphi (lep,W)" ,[2, 0, 3.142],"LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["dPhi_4bins", "fabs(dPhi)", "#Delta #varphi (lep,W)" ,[4, 0, 3.142],"LogY",["MoreY",2000],"IncludeOverflows"])

varList.append(["nVtx", "nVtx", "Nvert", [100, 0, 100], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["iso_MT2", "iso_MT2", "MT_{2} [GeV]", [20, 0, 1000], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["Lep_miniIso", "Lep_miniIso", "MiniISO(l) [GeV]", [50, 0, 0.5], "LogY",["MoreY",2000],"IncludeOverflows"])

varList.append(["Lp_log", "Lp", "L_{p}", [30, -0.5, 2.5], "LogY",["MoreY",2000],"IncludeOverflows"])
varList.append(["Lp", "Lp", "L_{p}", [30, -0.5, 2.5],["MoreY",1.6],"IncludeOverflows"])

varList.append(["HT_AN", "HT", "H_{T} [GeV]", [42, 0, 3000], "LogY",["MoreY",2000],"IncludeOverflows",['varbin', [500,750,1000,1250,2500],True]])
varList.append(["LT_AN", "LT", "L_{T} [GeV]", [35, 0, 1200] , "LogY",["MoreY",2000],"IncludeOverflows",['varbin', [250,350,450,600,1000],True]])
varList.append(["dPhi_AN", "fabs(dPhi)", "#Delta #varphi (lep,W)" ,[30, 0, 3.142],"LogY",["MoreY",2000],"IncludeOverflows",['varbin', [0,0.1,0.2,0.3,0.4,0.5,0.75,1,1.5,2.0,3.142],True]])
