import scanpy as sc
import sys
import matplotlib.pyplot as plt
import numpy as np
def find_badcells(adata,MT_cutoff=0.18, gpc_cutoff=2000, UMIpc_cutoff=0,
    do_plot='True'):
    
    sys.stderr.write('Removing bad cells... ')
    
    # remove cells with high MT genes transcription 
    MT_symb = ['MT-ND1','MT-ND2','MT-ND3','MT-ND4','MT-ND4L','MT-ND5','MT-ND6',
        'MT-CYB','MT-CO1','MT-CO2','MT-CO3','MT-ATP6','MT-ATP8','MT-RNR2']
    MT = []
    for mt in MT_symb:
        try:
            MT.append(adata.var.loc[mt].name)
        except KeyError:
            pass
          
    find = lambda searchList, elem: [np.where(searchList==e)[0][0] for e in elem]
    
    
    pMT = adata[:,MT].X.todense().sum(axis=1)/adata.X.todense().sum(axis=1)
     
    #gbm.matrix[:,find(gbm.ensembl,MT)].sum(axis=1)/gbm.matrix.sum(axis=1)
    pMT = np.squeeze(np.asarray(pMT.T))
    plotpMT = np.sort(pMT)[::-1]
  
    # save plot for MT distribution
    if do_plot.lower()=='True'.lower():
        plt.figure()
        plt.plot(plotpMT)

        try:
            plt.axvline(x=np.where(plotpMT>MT_cutoff)[0][-1],color='k')
        except IndexError:
            pass
        plt.ylim(0,1)
        plt.xlabel('cell index')
        plt.ylabel('MT proportion')
        plt.savefig('MT_distribution.png')
    
    # genes_per_cell
    genes_per_cells = np.squeeze(np.asarray(
        adata.X.todense().astype(bool).sum(axis=1)))
    pltgpc = np.sort(genes_per_cells)[::-1]
    if do_plot.lower()=='True'.lower():
      plt.figure()
      plt.plot(pltgpc,c='blue')
    
      try:
          plt.axvline(x=np.where(pltgpc>gpc_cutoff)[0][-1],color='k')
      except IndexError:
          pass
      plt.title('genes per cell')
      plt.xlabel('cell')
      plt.ylabel('number of different genes')
      plt.savefig('genes_per_cell.png',bbox_inches='tight',dpi=100)  
  
    # UMI per cell  
    UMI_per_cell = np.squeeze(np.asarray(adata.X.todense().sum(axis=1)))
    pltUpc = np.sort(UMI_per_cell)[::-1]
    if do_plot.lower()=='True'.lower():
      plt.figure()
      plt.subplot(2,1,1)
      plt.hist(UMI_per_cell, bins=100)
      plt.ylabel('number of cells')
      plt.xlabel('UMI count')
      plt.subplot(2,1,2)
      plt.plot(pltUpc)
      try:
          plt.axvline(x=np.where(pltUpc>UMIpc_cutoff)[0][-1],color='k')
      except IndexError:
          pass
      plt.ylabel('UMI count')
      plt.xlabel('cells')
      plt.savefig('UMI_per_cell.png',bbox_inches='tight',dpi=100)
            
    sys.stderr.write('done\n')
    return (pMT<MT_cutoff) & (genes_per_cells>gpc_cutoff) & (UMI_per_cell>UMIpc_cutoff)
