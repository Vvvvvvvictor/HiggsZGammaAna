print()
print("################################")
print("Starting testJSONfiles.py")
print("I will check all .json files in this folder")
print()

from correctionlib.schemav1 import CorrectionSet
from pydantic import ValidationError

def testJSONfile( fileName ):
    print()
    print(f"{fileName} -> test started")
    try:
        obj = CorrectionSet.parse_file(fileName)
    except ValidationError:
        # Possibly print error text
        print(f"{fileName} -> test failed")
        return False
    print(f"{fileName} -> test passed")
    return True

import glob, os
onlyfiles = glob.glob("POG/*/*.json")
failedTests = []
for fileName in onlyfiles:
    if os.stat(fileName).st_size == 0:
        print("Skipping empty file %s"%fileName)
        continue
    if not testJSONfile(fileName): failedTests.append(fileName)

if len(failedTests)>0: 
    print("\ntestJSONfiles.py is FAILED for the following files: %s \n"%str(failedTests))
    CorrectionSet.parse_file(failedTests[0])

print()
print("testJSONfiles.py is DONE without errors.")
print("################################")
print()
