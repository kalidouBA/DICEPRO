from sys import *
"""
TODO: 
    verify that all cell type labels are properly interpreted
    Manually check addition is done correctly
"""
"""
Purpose of this script:
Depending on the tool, reference, and mixture, the columns of the
outputs may need to me reordered or combined to match the answerkey
"""

def ConvertNames(namesList):
    newNames = []
    for oldName in namesList:
       newNames.append(Convert1Name(oldName))
    return newNames

def Convert1Name(Name):
  name = Name.lower()
  if "endothelial" in name\
          or "endo" in name:
     return "Other"
  if "neutrophil" in name:
     return "Neutrophils"
  if "macrophage" in name:
     return "Other"
  if "fibroblast" in name\
     or "cafs" in name:
     return "Other"
  if "cd4" in name\
     or "cd8" in name\
     or "T Cell" in name:
     return "Lymphocyte"
  if "mono" in name\
     or "cd14" in name:
     return "Monocytes"
  if "nk" in name\
          or "natural" in name:
     return "Lymphocyte"
  if "mast" in name:
     return "Other"
  if name == "b" or \
          "cd19b" in name\
          or "b cell" in name\
          or "b_cell" in name\
          or "b.cell" in name\
          or "bcell" in name\
          or "b lineage" in name\
          or "b-cell" in name:
     return "Lymphocyte"
  if name in ["p-value", "correlation", "rmse","absolute score (sig.score)"]:
     return "trash"
  return "Other"

def combineColumns(predMat, colsToCombine):
    ConvertedMat = []
    for lineParts in predMat:
      ConvertedLine = [lineParts[0]]
      for group in colsToCombine:
        groupSum = 0
        for col in group:
          try:
            groupSum += float(lineParts[col+1])
          except:
            continue
        ConvertedLine.append(str(groupSum))
      ConvertedMat.append(ConvertedLine)
    return ConvertedMat

def ReformatPredictions(predictionsFile, AnsCTs):
    """
    Identifies columns in predictions file that correspond to the columns of the answer file
    in the case of multiples pred CTS that map to the same AnsCT, takes the sum
    """
    predMat = []
    for line in open(predictionsFile,"r"):
      splitLine = line.strip().split("\t")
      if len(splitLine) == 1:
         splitLine = line.strip().split(",")
      predMat.append(splitLine)


    if len(predMat[0]) == len(predMat[1]):
       #i.e. something like "gene" fills top left corner
       convertedNames = ConvertNames(predMat[0][1:])
    else:
       convertedNames = ConvertNames(predMat[0])

    #initialize list of columns to combine
    colsToCombine = []
    for i in range(len(AnsCTs)):
        colsToCombine.append([])

    for i in range(len(AnsCTs)):
      AnsCT = AnsCTs[i]
      for j in range(len(convertedNames)):
        predCT = convertedNames[j]
        if predCT == AnsCT:
          colsToCombine[i].append(j)

    Reformatted = combineColumns(predMat[1:], colsToCombine)
    return Reformatted       

    
def main():
   AnsFile = open(argv[2],"r")
   for line in AnsFile:
     AnsCTs = line.strip().split("\t")[1:]
     break
   NewPreds = ReformatPredictions(argv[1],AnsCTs + ["Other"])
   print "\t" + "\t".join(AnsCTs + ["Other"])
   for line in NewPreds:
      print "\t".join(line)
main()
