library("ggplot2")
library("bootnet")

nBoots <- 1

# read in data
Data <- read.csv("C:\\Users\\HearneL\\AnacondaProjects\\FBR_Trauma_networks\\TraumaEffect_data_wAge_prepro.csv", stringsAsFactors = FALSE)

# remove unwanted variables
Data$X <- NULL
Data$ID <- NULL

# generate network
Network <- estimateNetwork(Data,default = "EBICglasso")

# plot network
Groups <- c(rep("CIDI",3),rep("Kessler",1),
            rep("TSCL",6),rep("PTCI",3))

tiff(file = "C:\\Users\\HearneL\\AnacondaProjects\\FBR_Trauma_networks\\Net1.tiff",
     width = 1700, 
     height = 900, 
     units = "px", 
     pointsize = 12,
     compression = c("none"),
     bg = "transparent", 
     res = 300)

plot(Network,
     layout = "spring",
     labels = TRUE,
     nodeNames = colnames(Data),
     groups = Groups,
     legend.cex = 0.3, 
     vsize = 5,
     esize = 15, 
     palette = 'pastel',
     posCol = "#FF9933",
     negCol = "#84c4ff", 
     borders = FALSE, 
     vTrans = 200,
     details = FALSE)
dev.off()

library("qgraph")
centralityPlot(Network)

# bootstrap test for CIs of edge weights
boot1 <- bootnet(Network, nBoots=nBoots, type = "nonparametric")
plot(boot1,labels = TRUE, order = "sample")
CS <- corStability(boot1)

# bootstrap test for CIs of centrality measures
boot2 <- bootnet(Network, nBoots=nBoots, type = "case")
plot(boot2)
CS <- corStability(boot2)