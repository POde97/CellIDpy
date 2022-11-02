import scipy.stats as stats
import mapply
import scanpy as sc
import pandas as pd
import numpy as np
import scanpy as sc
import scipy
#mapply.init(
 #   n_workers=-1,
  #  chunk_size=1,
  #  max_chunks_per_worker=16,
  #  progressbar=True
#)
class Hypergeom(object):
  def __init__(self,ad,reference,prediction = True,NUniverse = None ):

     
    
    if isinstance(reference, pd.DataFrame) == True:
      c = reference
    else:
      c = reference.obs

    

    #Distance adata
    self.Lsad = len(ad.obs["signature"].iloc[0])
    
    Distance1 = pd.DataFrame(self.GeneCellDistance(ad.obsm["X_MCA"].T,ad.varm["GenesCoordinates"]),columns=ad.var_names,index=ad.obs.index)
    #Indicator matrix1
    Distance1 = Distance1.mask(Distance1.rank(axis=1, method='max', ascending=True) > self.Lsad, 0)
    Distance1 = Distance1.mask(Distance1!=0,1)


    #Create dummy-IndicatorMatrix matrix for Reference
    
    pets_iter = (set(p) for p in c["signature"])
    pets = sorted(set.union(*pets_iter))
    dummies = pd.DataFrame(np.zeros((len(c), len(pets)), dtype=np.int), columns=pets)
    for i, pet in enumerate(c["signature"]):
      dummies.loc[i, pet] = 1
    dummies.index = list(c.index)
    


    #Allineate dummy and matrix 1
    
    #Order dummies genes in the same way of Distance1
    Allign = pd.DataFrame([0 for z in range(len(Distance1.T))],index = list(Distance1.columns),columns = ["todrop"])
    dummies = pd.concat([Allign,dummies.T],axis=1).fillna(0)
    dummies = dummies.drop("todrop",axis=1).T
    
    dummies = dummies.loc[:,dummies.columns.isin(list(Distance1.columns))]


    #Allign = pd.DataFrame([0 for z in range(len(dummies.T))],index = list(dummies.columns),columns = ["todrop"])
    #Distance1 = pd.concat([Distance1.T,Allign],axis=1).fillna(0)
    #Distance1 = Distance1.drop("todrop",axis=1).T

    Intersection = pd.DataFrame(np.matmul(Distance1.to_numpy(),dummies.to_numpy().T),index = list(Distance1.index),columns = list(dummies.index))
    

    lenS = pd.DataFrame(c["signature"].apply(lambda x: len(self.intersection(x,list(ad.var_names))))).T
    #Add len of reference gene signature as first row of Intersection
    Intersection = pd.concat([lenS,Intersection])
    
    if NUniverse == None:
      self.NUniverse = np.shape(Distance1)[1]
      
    else:
      self.NUniverse = NUniverse 
      
    


    self.Intersection = Intersection.mapply(lambda x: self.HyperG(x))
    self.Intersection.index = Distance1.index
    self.Intersection = -np.log10(self.Intersection.apply(lambda x: self.BHcorrection(x)))


     
    pred = pd.DataFrame(self.Intersection.apply(lambda x: x.nlargest(1).index.tolist()[0], axis=1),columns=["gs_prediction"])
    pval = pd.DataFrame(self.Intersection.apply(lambda x: x.max(), axis=1),columns=["pred_pval"])

      #Save prediction in adata 
    if prediction ==True:
      ad.obs["gs_prediction"] = pred["gs_prediction"].copy()
      ad.obs["pred_Pval"] = pval["pred_pval"].copy()
      #Non significant prediction 
      ad.obs.loc[ad.obs["pred_Pval"] <= 2.0, 'gs_prediction'] = "unassigned" 
      

    




  def GeneCellDistance(self,cellC,FeaturesCoordinates):
    Dist = scipy.spatial.distance_matrix(cellC.T, FeaturesCoordinates, p=2)
    return Dist
  
  def intersection(self,lst1, lst2):
    return list(set(lst1) & set(lst2))

  

  def HyperG(self,x):
    #lst3 = [value for value in lst1 if value in lst2]
    #g = 200 ## Number of submitted genes
    k = self.Lsad ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
    m = x[0] ## Number of "marked" elements, i.e. genes associated to this biological process
    N = int(self.NUniverse) ## Total number of genes with some annotation in GOTERM_BP_FAT.
    n = N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
    x = x[1:]## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process


    Pval = stats.hypergeom(M=N, 
                  n=m, 
                  N=k).sf(x-1)
# 4.989682834451419e-12
    return Pval


  def BHcorrection(self,x):
    d = x.sort_values(ascending = False)*len(x)
    x = d/pd.DataFrame([i for i in range(1,len(d)+1)],index = d.index)[0]
    return x
