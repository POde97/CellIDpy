import scanpy as sc
from sklearn.utils.extmath import randomized_svd
import pandas as pd
import numpy as np
from sklearn.utils.extmath import randomized_svd
from scipy.spatial import distance
import scipy
import heapq

class CellID():

  def __init__(self, ad,ngenesxcell=200):

		
    self.n = ngenesxcell
    self.ad = ad 
    if scipy.sparse.issparse(self.ad.X)== True:
      self.X = pd.DataFrame.sparse.from_spmatrix(self.ad.X.T,index =self.ad.var_names, columns=self.ad.obs.index)
    else:  
      self.X = pd.DataFrame(self.ad.X.T,index =self.ad.var_names, columns=self.ad.obs.index )
    #self.X = pd.DataFrame(ad.X.T,index =ad.var_names, columns=ad.obs.index )
    
    self.X = self.X.loc[self.X.var(axis=1) !=0]
    bool_series = self.X.index.duplicated(keep = False)

    self.X = self.X.iloc[~bool_series]
    self.X = self.X.reset_index().drop("index",axis=1)
    
    self.ad = self.ad[:,list(self.X.index)]
    
    
    print("Computing Fuzzy Matrix")
    self.Dc,self.Z,self.FM = self.MCAStep1()
    print("Computing SVD")
    U,s,self.V = randomized_svd(self.FM.to_numpy(), 
                              n_components=50,
                              random_state=None)
    print("Computing coordinates")
    self.CoordinatesCell,self.CoordinatesGenes = self.MCAStep2()
    print("Computing Cell-Genes Distances")
    self.Dist = self.GeneCellDistance()
    print("Build signature with TOP "+str(ngenesxcell)+" closest genes")
    self.signaturedf = self.XcellGeneS()

    print("Storing MCA in adata object")
    self.ad.obsm["X_MCA"] = self.cellC.T.copy()
    self.ad.varm['GenesCoordinates'] = self.FeaturesCoordinates.copy()
    self.ad.obs["signature"] = self.signaturedf["signature"].copy()
    
    del self.signaturedf,self.cellC

  def MCAStep1(self): 
    rmin = self.X.min(axis=1)
    rmax = self.X.max(axis=1)
    range = rmax - rmin
    self.X  = self.X.subtract(rmin,axis=0)
    self.X = self.X.divide(range,axis=0)
    self.FM = pd.concat([self.X,1-self.X],axis= 0)
    self.FM = self.FM.reset_index().drop("index",axis=1)

    total = self.FM.sum().sum()
    colsum = self.FM.sum(axis=0)
    rowsum = self.FM.sum(axis=1)

    self.FM = self.FM.divide(np.sqrt(rowsum),axis=0)

    self.Z = self.FM.divide(np.sqrt(colsum),axis=1)
    self.Dc = 1/(np.sqrt(rowsum/total))
    return self.Dc,self.Z,self.FM
  
  def MCAStep2(self):
    FeatureC = np.matmul(self.Z.to_numpy(),self.V.T)
    Zcol = np.shape(self.FM)[1]
    del self.Z
    FeaturesCoordinates = pd.DataFrame(FeatureC).multiply(self.Dc,axis=0)
    self.cellC = np.sqrt(Zcol)*self.V
    self.FeaturesCoordinates=FeaturesCoordinates.to_numpy()[:int(len(FeaturesCoordinates)/2),:]
    return self.cellC,self.FeaturesCoordinates

  def GeneCellDistance(self):
    Dist = scipy.spatial.distance_matrix(self.cellC.T, self.FeaturesCoordinates, p=2)
    return Dist

  def XcellGeneS(self):
    
    df = pd.DataFrame(self.Dist,columns =self.ad.var_names, index=self.ad.obs.index )
    self.signaturedf = pd.DataFrame(df.apply(lambda x: x.nsmallest(int(self.n)).index.tolist(), axis=1),columns=["signature"])
    
    return self.signaturedf
