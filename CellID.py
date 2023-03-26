import scanpy as sc
from sklearn.utils.extmath import randomized_svd
import pandas as pd
import numpy as np
from sklearn.utils.extmath import randomized_svd
from scipy.spatial import distance
import scipy
import heapq
from scipy.sparse.linalg import svds
import warnings
from sklearn.neighbors import NearestCentroid
import numpy as np
    
warnings.simplefilter(action='ignore', category=pd.errors.PerformanceWarning)

#Arguments
#ad: anndata object 
#ngenesxcell : number of genes per cell for per cell gene signature
#dim: dimension of Multiple correspondence analysis embedding
#verbose: progress bar and status
#Group_key: calculate per group gene signature using metadata (see CellID (g) in the CellID paper) 

class CellID():

  def __init__(self, ad,ngenesxcell=200,dim=50,verbose = True,group_key=None):
    
    self.n = ngenesxcell
    self.ad = ad.copy() 
    if group_key !=None:
      self.clusterS = list(ad.obs[group_key])

    if scipy.sparse.issparse(self.ad.X)== True:
      self.X = pd.DataFrame.sparse.from_spmatrix(self.ad.X.T,index =self.ad.var_names, columns=self.ad.obs.index).copy()
    else:  
      self.X = pd.DataFrame(self.ad.X.T,index =self.ad.var_names, columns=self.ad.obs.index ).copy()
    
    self.X = self.X.fillna(0.)

    self.X = self.X.astype(np.float64)
    self.X.replace([np.inf, -np.inf], 0., inplace=True)

    self.X = self.X.loc[self.X.var(axis=1) !=0].copy()
    bool_series = self.X.index.duplicated(keep = False).copy()

    self.X = self.X.iloc[~bool_series].copy()
    self.ad = self.ad[:,list(self.X.index)].copy()
    
    if verbose == True:
      print("Computing Fuzzy Matrix")
    self.Dc,self.Z,self.FM = self.MCAStep1()
    if verbose == True:
      print("Computing SVD")
    #sklearn Svd with ordered dimensions
    U,s,self.V = randomized_svd(self.FM.to_numpy(), 
                              n_components=dim+1,
                              random_state=None)
    self.V = self.V[1:,:].copy()
    
    #scipy svd not ordered (seems to be in reversed order ) 
    #U,s,self.V1 = svds(self.FM.to_numpy(),k=51)#,return_singular_vectors="vh",maxiter =1000,tol = 0.00001
    #self.V = self.V[:-1,:].copy()
    
    #irlbpy svd as R cellid 
    #U,s,self.V,it1,it2 = irlb(self.FM.to_numpy(),51,tol=0.00001,maxit=1000)
    #print(np.shape(self.V))
    #self.V = self.V[:,1:]
    #self.V = self.V.T[1:,:]
    if verbose == True:
      print("Computing coordinates")
    self.CoordinatesCell,self.CoordinatesGenes = self.MCAStep2()
    if verbose == True:
      print("Computing Cell-Genes Distances")
    self.Dist = self.GeneCellDistance()
    #print(np.shape(self.Dist))
    if verbose == True:
      print("Build signature with TOP "+str(ngenesxcell)+" closest genes")
    self.signaturedf = self.XcellGeneS()
    if verbose == True:
      print("Storing MCA in adata object")
    self.ad.raw = self.ad
    self.ad.obsm["X_MCA"] = self.cellC.T.copy()#np.fliplr(
    self.ad.varm['GenesCoordinates'] = self.FeaturesCoordinates.copy()#np.fliplr(
    self.ad.obs["signature"] = self.signaturedf["signature"].copy()
    
    

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

    self.FM = self.FM.divide(np.sqrt(colsum),axis=1)
    self.FM = self.FM.divide(np.sqrt(rowsum),axis=0)
    self.Z = self.FM.copy()
    
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
    Dist = scipy.spatial.distance_matrix(self.cellC.T.astype(np.double), self.FeaturesCoordinates.astype(np.double), p=2)
    return Dist
  




  def XcellGeneS(self):
    
    df = pd.DataFrame(self.Dist,columns =self.X.index, index=self.X.columns )
    self.signaturedf = pd.DataFrame(df.apply(lambda x: x.nsmallest(int(self.n)).index.tolist(), axis=1),columns=["signature"])
    
    return self.signaturedf   

  def XclusterGeneS(self):
      
    
    Dist = scipy.spatial.distance_matrix(self.cellC.T.astype(np.double), self.FeaturesCoordinates.astype(np.double), p=2)
    df = pd.DataFrame(self.Dist,columns =self.X.index, index=self.X.columns )

    clf = NearestCentroid()
    #centroid coordinates : n clustre x ngenes 
    centroid_coordinates = clf.fit(df.values, self.clusterS).centroids_
    
    centr_df = pd.DataFrame(centroid_coordinates,columns = list(df.columns),index = list(np.unique(self.clusterS)))
    signature_c_df = pd.DataFrame(centr_df.apply(lambda x: x.nsmallest(int(self.n)).index.tolist(), axis=1),columns=["signature"])

    return signature_c_df





