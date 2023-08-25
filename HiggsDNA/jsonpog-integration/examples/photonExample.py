## example how to read the photon format v2
from correctionlib import _core

evaluator = _core.CorrectionSet.from_file('./../POG/EGM/photon.json')

#Medium ID Run-2 scale factor
valsf= evaluator["UL-Photon-ID-SF"].evaluate("2017","sf","Medium",1.1, 34.0)
print("sf is:"+str(valsf))

#Loose ID Run-2 syst
valsyst= evaluator["UL-Photon-ID-SF"].evaluate("2016postVFP","syst","Loose",1.1, 34.0)
print("syst is:"+str(valsyst))

#Loose ID Run-2 CSEV scale factor
valsf= evaluator["UL-Photon-CSEV-SF"].evaluate("2016postVFP","sf","Loose","EBInc")
print("sf is:"+str(valsf))

#Loose ID Run-2 CSEV scale factor
valsf= evaluator["UL-Photon-PixVeto-SF"].evaluate("2016postVFP","sf","Loose","EBInc")
print("sf is:"+str(valsf))

#Loose ID Run-2 CSEV systematic uncertainty
valsf= evaluator["UL-Photon-PixVeto-SF"].evaluate("2016postVFP","syst","Loose","EBInc")
print("sf is:"+str(valsyst))

