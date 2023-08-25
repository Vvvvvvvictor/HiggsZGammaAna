## example how to read the electron format v2
from correctionlib import _core

evaluator = _core.CorrectionSet.from_file('./../POG/EGM/2016postVFP_UL/electron.json')

#Reconstruction (pT< 20 GeV) Run-2 scale factor
valsf= evaluator["UL-Electron-ID-SF"].evaluate("2016postVFP","sf","RecoBelow20",1.1, 15.0)
print("sf is:"+str(valsf))

#Reconstruction (pT> 20 GeV) Run-2 scale factor
valsf= evaluator["UL-Electron-ID-SF"].evaluate("2016postVFP","sf","RecoAbove20",1.1, 25.0)
print("sf is:"+str(valsf))

#Medium ID Run-2 systematic uncertainty
valsyst= evaluator["UL-Electron-ID-SF"].evaluate("2016postVFP","syst","Medium",1.1, 34.0)
print("syst is:"+str(valsyst))
