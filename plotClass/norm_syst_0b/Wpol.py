
Wpol = {
    'DiLepTT' :
        {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF*nISRttweight',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF*nISRttweight',
        },

    'SemiLepTT' : 
        {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF*nISRttweight',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF*nISRttweight',
        },

    'SingleT' : 

        {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },

    'VV' : 
         {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },

    'TTV' : 
         {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },

    'QCD' : 
         {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },
        
    'WJ' : 
        {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },

    'DY' : 
         {
        'scale_up' : '1000.0/sumOfWeights2*WpolWup*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        'scale_dn' : '1000.0/sumOfWeights2*WpolWdown*genWeight*Xsec*1*btagSF*puRatio*lepSF',
        },
    'Signal_1' : 
        {
        'scale_up' : '1000.0*genWeight*susyXsec/susyNgen*btagSF*lepSF',
        'scale_dn' : '1000.0*genWeight*susyXsec/susyNgen*btagSF*lepSF',
        },
}
