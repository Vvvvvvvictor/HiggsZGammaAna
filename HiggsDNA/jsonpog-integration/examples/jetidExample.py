#! /usr/bin/env python
# Example of how to read the JetID JSON file
# For more information, see the README in
# https://gitlab.cern.ch/cms-nanoAOD/jsonpog-integration/-/tree/master/POG/JME
from correctionlib import _core

# Load CorrectionSet
fname = "../POG/JME/2022_Summer22/jetid.json.gz"
if fname.endswith(".json.gz"):
    import gzip
    with gzip.open(fname,'rt') as file:
        data = file.read().strip()
        evaluator = _core.CorrectionSet.from_string(data)
else:
    evaluator = _core.CorrectionSet.from_file(fname)


##### Jet Identification (JetID)
eta = -0.23
chHEF = 0.4
neHEF = 0.25
chEmEF = 0.1
neEmEF = 0.15
muEF = 0.1
chMultiplicity = 1.0
neMultiplicity = 2.0
multiplicity = chMultiplicity + neMultiplicity
print(f"for a jet with eta: {eta}, chHEF: {chHEF}, neHEF: {neHEF}, chEmEF: {chEmEF}, neEmEF: {neEmEF}, muEF: {muEF}, chMultiplicity: {chMultiplicity}, neMultiplicity: {neMultiplicity}, multiplicity: {multiplicity}")

criteria_name = "AK4PUPPI_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

print("="*50)
print()

eta = 4.5
chHEF = 0.4
neHEF = 0.25
chEmEF = 0.1
neEmEF = 0.15
muEF = 0.1
chMultiplicity = 1.0
neMultiplicity = 2.0
multiplicity = chMultiplicity + neMultiplicity
print(f"for a jet with eta: {eta}, chHEF: {chHEF}, neHEF: {neHEF}, chEmEF: {chEmEF}, neEmEF: {neEmEF}, muEF: {muEF}, chMultiplicity: {chMultiplicity}, neMultiplicity: {neMultiplicity}, multiplicity: {multiplicity}")

criteria_name = "AK4PUPPI_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

criteria_name = "AK4CHS_TightLeptonVeto"
result = evaluator[criteria_name].evaluate(eta, chHEF, neHEF, chEmEF, neEmEF, muEF, chMultiplicity, neMultiplicity, multiplicity)
print(f"jetid: {criteria_name} is {result} ({'pass' if result else 'fail'})")

print("="*50)
