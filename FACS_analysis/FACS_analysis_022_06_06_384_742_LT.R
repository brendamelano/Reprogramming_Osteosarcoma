library(flowCore)
library(latticeExtra)
library(ggcyto)
library(tidyverse)
library(knitr)
library(dplyr)
library(flowWorkspace)
library(svglite)
library(ggplot2)


# Read in flow data from the FACS sort for OS384 and OS742 for the pilot LT experiment
fs <- read.flowSet(path = "~/Desktop/FACS_analysis/FACS_2022_06_06_384_742_LT/", pattern = ".fcs", alter.names = T)



##############    OS384     #############



# Running analysis for only OS384 first
fs <- fs[1:2]


# Printing the fcs file specimens
# Should be only the 384 samples
pData(fs)[1:2,]


# Removing specimen prefix
pData(fs)$well <- gsub("Specimen_001_","",sampleNames(fs)) 


# removing the .fcs suffix
pData(fs)$well <- gsub(".fcs","", pData(fs)$well) 


# printing the colnames or colors of the different colors assayed
colnames(fs)


# Changing the pacific blue color to BFP
colnames(fs)[colnames(fs) == "Pacific.Blue.A"] <- "BFP"


# converting the flowset object to a gating set object
gs <- GatingSet(fs)


# defining the gate for singlets
g.singlets <- polygonGate(filterId = "Singlets",
                          "FSC.A"=c(1.4e4,16e4,14e4,.5e4),
                          "FSC.H"=c(0.9e4,10e4,15e4,1.5e4)) 


# Checking where the gate plots
ggcyto(gs[[1]], aes(x=FSC.A, y=FSC.H), subset="root") + 
  geom_hex(bins = 200)+
  geom_gate(g.singlets)+
  ggcyto_par_set(limits = "instrument") +
  theme_bw()#+
  #geom_point(color="blue", shape=18 ,size=8)


# Add gate to GatingSet
gs_pop_add(gs, g.singlets) 


# Recompute GatingSet
recompute(gs) 


# Plotting the singlets
ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets") +
  geom_hex(bins = 300) +
  ggcyto_par_set(limits = list(x = c(0,1.2e5), y = c(-1, 3e4))) +
  theme_bw()


# plotting all the singlets
ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="root") + 
  geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8) +
  facet_wrap(~well,ncol = 10)


# gating on the live cells
g.live <- polygonGate(filterId = "Live","FSC.A"=c(5e3,1e5,1e5,8e3),"SSC.A"=c(0,2e3,5e4,1e4)) # define gate


# plotting the live cells
ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets") +
  geom_hex(bins = 300) + 
  geom_gate(g.live) +
  ggcyto_par_set(limits = list(x = c(0,1.2e5), y = c(-1, 3e4)))


# adding the live gate to the gating set
gs_pop_add(gs,g.live,parent="Singlets") 


# recompute GatingSet
recompute(gs) 


# plotting all the singlets
ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="Singlets")+geom_hex(bins = 100)+geom_gate("Live")+
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument")+
  facet_wrap(~well,ncol = 10)


# set gate
g.gfp <- rectangleGate(filterId="BFP positive","BFP"=c(2000, Inf)) 


# adding the gfp gate
gs_pop_add(gs, g.gfp,parent="Live") # add gate to GatingSet


# recomputing the gating set object after filtering on GFP
recompute(gs) # recalculate Gatingset


# plotting the control and BFP positive cells
ggcyto(gs, aes(x = BFP), subset = "Live") +
  geom_density(fill = "forestgreen", alpha = 0.5) +  # Set the alpha value to control transparency
  geom_gate("BFP positive") +
  ggcyto_par_set(limits = "instrument") +
  scale_x_flowjo_biexp() +
  facet_wrap(~ well, ncol = 4) + 
  theme_bw() +
  geom_stats(
    size = 3,
    color = "red",
    adjust = 0.6,
    label.format = format,
    label.args = list(digits = 2)
  )


ggcyto(gs, aes(x = BFP), subset = "Live") +
  geom_histogram(fill = "forestgreen", alpha = 0.5, position = "identity") +
  geom_gate("BFP positive") +
  ggcyto_par_set(limits = "instrument") +
  scale_x_flowjo_biexp() +
  theme_bw() +
  #facet_null() +
  geom_stats(
    size = 3,
    color = "red",
    adjust = 0.6,
    label.format = format,
    label.args = list(digits = 2)
  )

  
ggcyto(gs, aes(x = BFP), subset = "Live") + 
  geom_density(aes(y = ..count..)) + 
  geom_gate("BFP positive") + facet_null()

fs_ctrl <- gs_pop_get_data(gs[1], "Live")

fs_LT <- gs_pop_get_data(gs[2], "Live")

fs_ctrl_BFP <- gs_pop_get_data(gs[1], "BFP positive")
fs_ctrl_BFP > 2000


OS384_ctrl <- as.data.frame(OS384_ctrl)

OS384_ctrl <- subset(gs[1], subset = ('BFP positive'>0))

OS384_LT <- gs[2]


# saving a plot with the overlapped histograms
OS384_FACS_overlap <- ggplot() +
  geom_density(data = fs_ctrl, aes(x = BFP, fill = "Control"), alpha = 0.5) +
  geom_density(data = fs_LT, aes(x = BFP, fill = "LT"), alpha = 0.5) +
  scale_fill_manual(values = c("red", "blue"),
                    labels = c("Control", "LT")) +
  scale_x_flowjo_biexp() +
  geom_vline(xintercept = 2400, linetype = "dashed", color = "red") +
  ggtitle('OS384 BFP positive cells') +
  theme_bw() +
  theme(legend.title = element_blank()) 



print(OS384_FACS_overlap)

# Calculate percentage above threshold line
ctrl_above_threshold <- sum(fs_ctrl$BFP > 2400) / nrow(fs_ctrl) * 100
lt_above_threshold <- sum(fs_LT$BFP > 2400) / nrow(fs_LT) * 100

# Add percentage above threshold line
OS384_FACS_overlap <- OS384_FACS_overlap +
  annotate("text", x = 2400, y = 0.05, label = paste0("Control: ", percent(ctrl_above_threshold)),
           vjust = 0, hjust = -0.5, color = "black") +
  annotate("text", x = 2400, y = 0.03, label = paste0("LT: ", percent(lt_above_threshold)),
           vjust = 0, hjust = -0.5, color = "black")


# Save the plot as an SVG file
ggsave("~/Desktop/OS384_FACS_overlap.svg", plot = OS384_FACS_overlap, device = "svg")


percentage_ctrl <- sum(fs_ctrl[[BFP]] > threshold) / nrow(fs_ctrl) * 100
percentage_LT <- sum(fs_LT$BFP > threshold) / nrow(fs_LT) * 100

ggplot() +
  geom_density(data = fs_ctrl, aes(x = BFP), fill = "forestgreen", alpha = 0.5) +
  geom_density(data = fs_LT, aes(x = BFP), fill = "blue", alpha = 0.5) +
  geom_vline(xintercept = threshold, linetype = "dashed", color = "red") +
  geom_text(x = threshold, y = 0, label = paste0("Positive: ", round(percentage_ctrl, 2), "%"), 
            vjust = -1, hjust = -0.5, color = "forestgreen") +
  geom_text(x = threshold, y = 0, label = paste0("Positive: ", round(percentage_LT, 2), "%"), 
            vjust = -1, hjust = 1.5, color = "blue") +
  scale_x_flowjo_biexp() +
  theme_bw() +
  ggtitle('OS384 BFP positive cells')



ps <- gs_pop_get_count_with_meta(gs)

ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
ps %>% select(sampleName,well,Population,Count,ParentCount,percent_of_parent) %>% head() %>% kable

