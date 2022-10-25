import scipy.stats as stats
import mapply
import scanpy as sc
mapply.init(
    n_workers=-1,
    chunk_size=1,
    max_chunks_per_worker=8,
    progressbar=True
)
class Hypergeom(object):
  def __init__(self,ad,reference,prediction = True,NUniverse = None ):

    self.signaturedf = pd.DataFrame(ad.obs["signature"])
    self.intersectiondf=pd.DataFrame()
     
    
    if prediction == True:
      c = reference
      for i in range(len(c)): 
  #intersectiondf[c.iloc[i].name]=signaturedf.apply(lambda x: len(intersection(x["signature"],c.iloc[i]["official gene symbol"])) ,axis=1)
        self.intersectiondf[c.iloc[i].name]=self.signaturedf.apply(lambda x: self.HyperG(x["signature"],c.iloc[i]["official gene symbol"]) ,axis=1)

      pred = pd.DataFrame(self.intersectiondf.apply(lambda x: x.nlargest(1).index.tolist()[0], axis=1),columns=["gs_prediction"])
      pval = pd.DataFrame(self.intersectiondf.apply(lambda x: x.max(), axis=1),columns=["pred_pval"])

    #Save prediction in adata 
    
      ad.obs["gs_prediction"] = pred["gs_prediction"].copy()
      ad.obs["pred_Pval"] = pval["pred_pval"].copy()

    else:



      #Concatenate 2 matri Same number of genes needed

      ad.obsm["X_MCA"].T
      ad.varm["GenesCoordinates"]
      


      Distance1 = pd.DataFrame(self.GeneCellDistance(ad.obsm["X_MCA"].T,ad.varm["GenesCoordinates"]),columns=ad.var_names,index=ad.obs.index)
      #Indicator matrix1
      Distance1 = Distance1.mask(Distance1.rank(axis=1, method='max', ascending=True) > 200, 0)
      Distance1 = Distance1.mask(Distance1!=0,1)

      Distance2 = pd.DataFrame(self.GeneCellDistance(reference.obsm["X_MCA"].T,reference.varm["GenesCoordinates"]),columns=reference.var_names,index=reference.obs.index)
      #Indicator matrix1

      Distance2 = Distance2.mask(Distance2.rank(axis=1, method='max', ascending=True) > 200, 0)
      Distance2 = Distance2.mask(Distance2!=0,1)
      
      if NUniverse == None:
        
        self.NUniverse = max(np.shape(Distance2)[1],np.shape(Distance1)[1])
      else:
        
        self.NUniverse = NUniverse 
        
      if np.shape(Distance2)[1]>np.shape(Distance1)[1]:
        Distance2= Distance2[list(Distance1.columns)]
      elif np.shape(Distance1)[1]>np.shape(Distance2)[1]:
        Distance1= Distance1[list(Distance2.columns)]


      
      
      Intersection = pd.DataFrame(np.matmul(Distance1.to_numpy(),Distance2.to_numpy().T),index = list(Distance1.index),columns = list(Distance2.index))
  
  #intersectiondf[c.iloc[i].name]=signaturedf.apply(lambda x: len(intersection(x["signature"],c.iloc[i]["official gene symbol"])) ,axis=1)
      #self.intersectiondf[c.iloc[i].name]=self.signaturedf.mapply(lambda x: self.HyperG(x) ,axis=1)
      self.Intersection = Intersection.mapply(lambda x: self.HyperG(x))





  def GeneCellDistance(self,cellC,FeaturesCoordinates):
    Dist = scipy.spatial.distance_matrix(cellC.T, FeaturesCoordinates, p=2)
    return Dist
  
  def intersection(lst1, lst2):
    return list(set(lst1) & set(lst2))

  

  def HyperG(self,x):
    #lst3 = [value for value in lst1 if value in lst2]
    g = 200 ## Number of submitted genes
    k = 200 ## Size of the selection, i.e. submitted genes with at least one annotation in GO biological processes
    m = 200 ## Number of "marked" elements, i.e. genes associated to this biological process
    N = int(self.NUniverse) ## Total number of genes with some annotation in GOTERM_BP_FAT.
    n = N - m ## Number of "non-marked" elements, i.e. genes not associated to this biological process
    x = x## Number of "marked" elements in the selection, i.e. genes of the group of interest that are associated to this biological process


    Pval = stats.hypergeom(M=N, 
                  n=m, 
                  N=k).sf(x-1)
# 4.989682834451419e-12
    return -np.log10(Pval)
