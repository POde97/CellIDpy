def MgProteinCoding():
  import pickle
  with open("MgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 

def HgProteinCoding():
  import pickle
  with open("HgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 
