def MgProteinCoding():
  import pickle
  with open("data/MgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 

def HgProteinCoding():
  import pickle
  with open("data/HgProteinCoding", "rb") as fp:   # Unpickling
    l = pickle.load(fp)
  
  return l 
