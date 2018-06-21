library("ggplot2")
library("bootnet")

nBoots <- 100

# read in data
Data <- read.csv("C:\\Users\\HearneL\\AnacondaProjects\\FBR_Trauma_networks\\TraumaEvents_data.csv", stringsAsFactors = FALSE)

# remove unwanted variables
Data$X <- NULL
Data$ID <- NULL

#only want to include variables with at least 8 instances in this case.
Data1 = Data[,colSums(Data) > 20]

# generate network
Network <- estimateNetwork(Data1,default = "IsingFit")

# plot network
Groups <- c(rep("CIDI",8),rep("Cultural",8))
tiff(file = "C:\\Users\\HearneL\\AnacondaProjects\\FBR_Trauma_networks\\Net2.tiff",
     width = 1700, 
     height = 900, 
     units = "px", 
     pointsize = 14,
     compression = c("none"),
     bg = "transparent", 
     res = 300)

plot(Network,
     layout = "spring",
     labels = TRUE,
     nodeNames = colnames(Data1),
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

# bootstrap test for CIs of edge weights
boot1 <- bootnet(Network, nBoots=nBoots)
plot(boot1,labels = TRUE, order = "sample")
