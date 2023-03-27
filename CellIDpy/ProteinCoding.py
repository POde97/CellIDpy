def MgProteinCoding():
  import pickle
  with open("CellIDpy/CellIDpy/MgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 

def HgProteinCoding():
  import pickle
  with open("CellIDpy/CellIDpy/HgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 
