import matplotlib.pyplot as plt
import numpy as np
#MCA plot
#Simoultaneous plot of genes and cell in MCA space

#Arguments:
#ciD: CellID object
#cell_type_attr: x cell metadata (example cell-type) 
#genep: list of genes to plot (example: genep = ["CTRL", "INS", "MYZAP", "CDH11"])




#from CellIDpy import mcaplot
import matplotlib.pyplot as plt
import seaborn as sns

def MCAplot(ciD,cell_type_attr,genep=[]):
  fig,ax = plt.subplots(figsize=(8, 8))
  #genep = ["CTRL", "INS", "MYZAP", "CDH11"]
  idxg = []
  for gene in genep:
    idxg.append(list(ciD.ad.var_names).index(gene))

  cmap = ["blue","grey","green","yellow","orange","red","purple","brown",
          "plum","slategray","tan","floralwhite","sandybrown",
          "rosybrown","gold","olivedrab"]
  #if ciD.ad.obs[cell_type_attr].nunique() < len(cmap):
   # cmapz = cmap[:ciD.ad.obs[cell_type_attr].nunique()]
  #else:
  colors = sns.hls_palette(ciD.ad.obs[cell_type_attr].nunique())


  cmapz = colors.as_hex()

  colorr = list(ciD.ad.obs[cell_type_attr].map(dict(zip(list(ciD.ad.obs["cell.type"].unique()),cmapz))))
  ax.scatter(ciD.ad.obsm['X_MCA'][:,0],ciD.ad.obsm['X_MCA'][:,1],3,c=colorr)

  ax.plot(ciD.ad.varm["GenesCoordinates"][idxg,0],ciD.ad.varm["GenesCoordinates"][idxg,1],"x")

  for i,coord in enumerate(zip(ciD.ad.varm["GenesCoordinates"][idxg,0],ciD.ad.varm["GenesCoordinates"][idxg,1])):
    ax.annotate(genep[i], coord,size = 16)

  for v in zip(list(ciD.ad.obs["cell.type"].unique()),cmapz):
    ax.scatter([],[], c=v[1], label=v[0])

  ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))


  plt.ylabel("mca_2")
  plt.xlabel("mca_1")


#mcaplot.MCAplot(cID,"cell.type",genep = ["PAX6", "SOX18"])  
  
  
  
  



