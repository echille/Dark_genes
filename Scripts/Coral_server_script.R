#Set up workspace

getwd() #Display the current working directory

#If necessary, change the path below to the directory where the data files are stored. "." means current directory. On Windows use a forward slash / instead of the usual \.

workingDir = ".";
setwd(WGCNA);
library(WGCNA) #Load the WGCNA package
options(stringsAsFactors = FALSE) #The following setting is important, do not omit.
enableWGCNAThreads() #Allow multi-threading within WGCNA.

#Load the data saved in the first part

adjTOM <- load(file="datExpr.RData")
adjTOM

#Run analysis

softPower=6 #Set softPower to 6
adjacency=adjacency(datExpr, power=softPower,type="signed") #Calculate adjacency
TOM= TOMsimilarity(adjacency,TOMType = "signed") #Translate adjacency into topological overlap matrix
dissTOM= 1-TOM #Calculate dissimilarity in TOM
save(adjacency, TOM, dissTOM, file = "adjTOM.RData") #Save to bluewaves dir
save(dissTOM, file = ".dissTOM.RData") #Save to bluewaves dir

### Clustering using TOM

geneTree= flashClust(as.dist(dissTOM), method="average") Form distance matrix
pdf(file="dissTOMClustering.pdf", width=20, height=20)
plot(geneTree, xlab="", sub="", main= "Gene Clustering on TOM-based dissimilarity", labels= FALSE,hang=0.04)
dev.off()

### Module identification using dynamicTreeCut

minModuleSize = 30
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
deepSplit = 2, pamRespectsDendro = FALSE,
minClusterSize = minModuleSize)
table(dynamicMods) #list modules and respective sizes
save(dynamicMods, geneTree, file = "dyMod_geneTree.RData") #Save to bluewaves dir and scp to 3-WGCNA/Input dir

