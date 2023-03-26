import itertools
import statsmodels
import statsmodels.api as sm
import scipy.stats as stats
import scanpy as sc
import pandas as pd
import numpy as np
import scipy
import mapply
import torch

#Arguments:
#ad : anndata object of single cell dataset   
#reference: pandas dataframe (reference gene list) or anndata object for label transfer between dataset 
#NUniverse: custom NUniverse of Hypergeometric test. If None NUniverse == #genes in per cell gene signature of 
#           ad and reference 
#min_size: exclude reference gene list with a number of genes low than min_size 
#verbose: progress bar and operation
#gpu_cuda: if True use gpu for speed calculation.
#adj_pval: if true correct Pvalue for multiple testing with Benjamini and Hochberg correction
#log_trans: if true -log10(Pvalue) 
#mply_p: max_chunks_per_worker of mapply parallelization https://pypi.org/project/mapply/


class Hypergeom(object):
  def __init__(self,ad,reference,prediction = True,NUniverse = None,min_size = 10,verbose=True,gpu_cuda =None ,adj_pval =True,log_trans = True,mply_p=16):

     
    
    if isinstance(reference, pd.DataFrame) == True:
      c = reference
    else:
      c = reference.obs

    

    #Distance adata
    self.Lsad = len(ad.obs["signature"].iloc[0])
    
    #More efficient for Cell-ID but not for Cell2CellMatch (Old)
    #Distance1 = pd.DataFrame(self.GeneCellDistance(ad.obsm["X_MCA"].T,ad.varm["GenesCoordinates"]),columns=ad.var_names,index=ad.obs.index)
    #Indicator matrix1
    #Distance1 = Distance1.mask(Distance1.rank(axis=1, method='max', ascending=True) > self.Lsad, 0)
    #Distance1 = Distance1.mask(Distance1!=0,1)

    #Latest version : we create dummy for both so we can speed up Hyper and performe a single hyper Test
    pets_iter = (set(p) for p in ad.obs["signature"])
    pets = sorted(set.union(*pets_iter))
    Distance1 = pd.DataFrame(np.zeros((len(ad.obs), len(pets))), columns=pets)#, dtype=np.int
    for i, pet in enumerate(ad.obs["signature"]):
      Distance1.loc[i, pet] = 1
    Distance1.index = list(ad.obs.index)

    #Create dummy-IndicatorMatrix matrix for Reference
    if gpu_cuda == None:
      pets_iter = (set(p) for p in c["signature"])
      pets = sorted(set.union(*pets_iter))
      dummies = pd.DataFrame(np.zeros((len(c), len(pets)), dtype=np.int), columns=pets)
      for i, pet in enumerate(c["signature"]):
        dummies.loc[i, pet] = 1
      dummies.index = list(c.index)
      #Order dummies genes in the same way of Distance1
      Allign = pd.DataFrame([0 for z in range(len(Distance1.T))],index = list(Distance1.columns),columns = ["todrop"])
      dummies = pd.concat([Allign,dummies.T],axis=1).fillna(0)
      dummies = dummies.drop("todrop",axis=1).T
    
      dummies = dummies.loc[:,dummies.columns.isin(list(Distance1.columns))]
  
    if gpu_cuda == True:
      idx_inter = list(Distance1.index)
      col_inter = list(Distance1.index)
      #gpu cuda multiplication only for local topoly Network
      cuda = torch.device('cuda')
      Distance1 = torch.tensor(Distance1.to_numpy()).to(cuda)
      Distance1 = torch.mm(Distance1,torch.transpose(Distance1,0,1))
      Intersection = pd.DataFrame(Distance1.cpu().numpy(),index = idx_inter,columns=col_inter)
      del Distance1
    else: 
      Intersection = pd.DataFrame(np.matmul(Distance1.to_numpy(),dummies.to_numpy().T),index = list(Distance1.index),columns =  list(dummies.index))
    

    lenS = pd.DataFrame(c["signature"].apply(lambda x: len(self.intersection(x,list(ad.var_names))))).T
    #Add len of reference gene signature as first row of Intersection
    Intersection = pd.concat([lenS,Intersection])
    
    if NUniverse == None:
      self.NUniverse = np.shape(ad.X)[1]
    
    else:
      self.NUniverse = NUniverse 
      
    if verbose ==True:
      mapply.init(
        n_workers=-1,
        chunk_size=1,
        max_chunks_per_worker=mply_p,
        progressbar=True
      )
    else:
      mapply.init(
         n_workers=-1,
         chunk_size=1,
         max_chunks_per_worker=mply_p,
         progressbar=False
       )

    
    self.Intersection = Intersection.mapply(lambda x: self.HyperG(x))
    
    if gpu_cuda == None:
      self.Intersection.index = Distance1.index
    else:
      self.Intersection.index = idx_inter
    
    if adj_pval==True:
      self.Intersection = self.Intersection.apply(lambda x: self.BHcorrection(x))
    if log_trans ==True:
      self.Intersection = -np.log10(self.Intersection) 
    #self.Intersection = self.Intersection.apply(lambda x: self.BHcorrection(x),axis=1)

    
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

    return Pval

  def BHcorrection(self,x):
    return statsmodels.stats.multitest.multipletests(list(x), alpha=0.01, method='fdr_bh', is_sorted=False, returnsorted=False)[1]



  def correction(self,x):
    d = x.sort_values()*len(x)#
    x = d/pd.DataFrame([i for i in range(1,len(d)+1)],index = d.index)[0]
    x[x>1.] = 1.
    return x



