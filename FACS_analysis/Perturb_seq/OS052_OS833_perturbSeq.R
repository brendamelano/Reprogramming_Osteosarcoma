library(flowCore)
library(latticeExtra)
library(ggcyto)
library(tidyverse)
library(knitr)
library(dplyr)
library(flowWorkspace)
library(svglite)
library(ggplot2)


### OS8333

# read in flow data from the FACS sort for OS384 and OS742 for the pilot LT experiment
fs <- read.flowSet(path = "~/Desktop/23_06_12_052_833_pertur_seq_screen/OS833/", pattern = ".fcs", alter.names = T)


# removing the .fcs suffix
pData(fs)$well <- gsub(".fcs","", sampleNames(fs)) 


# Changing the pacific blue color to BFP
colnames(fs)[colnames(fs) == "Pacific.Blue.A"] <- "BFP"


# converting the flowset object to a gating set object
gs <- GatingSet(fs)


# defining the gate for singlets
# this is manual gating, try non-manual as well (high and low)
g.singlets <- polygonGate(filterId = "Singlets","FSC.A"= c(1.4e4,16e4,14e4,.5e4),
                          "FSC.H"=c(0.9e4,10e4,15e4,1.5e4)) 


# checking where the gate plots
# get rid of grid lines
ggcyto(gs[[1]], aes(x=FSC.A, y=FSC.H), subset="root") + 
  geom_hex(bins = 200)+
  geom_gate(g.singlets)+
  theme_bw() +
  ggcyto_par_set(limits = "instrument") 
   #+
#geom_point(color="blue", shape=18 ,size=8)


# add gate to GatingSet
gs_pop_add(gs, g.singlets) 


# recompute GatingSet
recompute(gs) 


# perhaps delete
# plotting the singlets
ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets") +
  geom_hex(bins = 300) +
  ggcyto_par_set(limits = list(x = c(0,1.2e5), y = c(-1, 3e3))) +
  theme_bw()


# plotting all the singlets
ggcyto(gs,aes(x=FSC.A,y=FSC.H),subset="root") + 
  geom_hex(bins = 100)+geom_gate("Singlets")+
  geom_stats(adjust = 0.8) +
  facet_wrap(~well,ncol = 10) +
  theme_bw()


# Creating a fate for the live cells
g.live <- polygonGate(filterId = "Live","FSC.A"=c(5e3,1e5,1e5,8e3),"SSC.A"=c(0,10,2e3,2e3))


# plotting the live cells
ggcyto(gs[[1]],aes(x=FSC.A,y=SSC.A),subset="Singlets") +
  geom_hex(bins = 300) + 
  geom_gate(g.live) +
  ggcyto_par_set(limits = list(x = c(0,1.2e5), y = c(-1, 3e3))) +
  theme_bw()


# adding the live gate to the gating set
gs_pop_add(gs,g.live,parent="Singlets") 


# recompute GatingSet
recompute(gs) 


# Plotting all the singlets
ggcyto(gs,aes(x=FSC.A,y=SSC.A),subset="Singlets") +
  geom_hex(bins = 100)+geom_gate("Live") +
  geom_stats(adjust = 0.8)+ggcyto_par_set(limits = "instrument") +
  facet_wrap(~well,ncol = 10) +
  ggcyto_par_set(limits = list(x = c(0,1.2e5), y = c(-1, 3e3)))


# set gate
g.gfp <- rectangleGate(filterId="BFP positive","BFP"=c(100, Inf)) 

# Changing the pacific blue color to BFP
colnames(gs)[colnames(gs) == "Pacific.Blue.A"] <- "BFP"

# adding the gfp gate
gs_pop_add(gs, g.gfp,parent="Live") 
# add gate to GatingSet

colnames(gs)

# checking the gate
ggcyto(gs[[1]], aes(x=BFP), subset="Live") +
  geom_density(fill="forestgreen") +
  geom_gate(g.gfp)+ggcyto_par_set(limits = "instrument") +
  scale_x_flowjo_biexp() +
  geom_stats( size = 6,  color = "red", fill = "green", adjust = 0.7) +
  theme_bw() +
  ggtitle('OS384 control')


# recomputing the gating set object after filtering on GFP
recompute(gs) # recalculate Gatingset


# plotting 
ggcyto(gs,aes(x=BFP),subset="Live") +
  geom_density(fill="forestgreen") +
  geom_gate(g.gfp) +
  ggcyto_par_set(limits = "instrument") +
  scale_x_flowjo_biexp() +
  facet_wrap(~well,ncol = 4) + 
  theme_bw() +
  geom_stats(size = 3,  color = "red", adjust = 0.6, label.format = function(x) sprintf("%.2f%%", x * 100))



# Calculate statistics outside the plot
stats <- ggcyto(gs, aes(x = BFP), subset = "Live") +
  geom_density(fill = "forestgreen") +
  geom_gate(g.gfp) +
  ggcyto_par_set(limits = "instrument") +
  scale_x_flowjo_biexp() +
  facet_wrap(~ well, ncol = 4) + 
  theme_bw()


# Create a data frame for custom labels
label_df <- data.frame(
  well = rep(levels(gs$well), each = 1),  # Replace `gs$well` with appropriate column name
  label = round(runif(length(levels(gs$well))), 2)  # Replace `gs$well` with appropriate column name
)

# Plot with custom labels
stats +
  geom_text(data = label_df, aes(label = paste0(label, "%")), color = "red", size = 4, vjust = -0.5)


ggcyto(gs, aes(x = BFP), subset = "Live") +
  geom_density(fill = "forestgreen", alpha = 0.5) +  # Set the alpha value to control transparency
  geom_gate(g.gfp) +
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
  geom_density(aes(y = ..count..)) + 
  geom_gate("BFP positive") + facet_null()

fs_ctrl <- gs_pop_get_data(gs[1], "Live")
fs_LT <- gs_pop_get_data(gs[2], "Live")


OS384_ctrl <- as.data.frame(OS384_ctrl)


OS384_ctrl <- subset(gs[1], subset = ('BFP positive'>0))


OS384_LT <- gs[2]


OS384_FACS_overlap <- ggplot() +
  geom_density(data = fs_ctrl, aes(x = BFP), fill = "forestgreen", alpha = 0.5) +
  geom_density(data = fs_LT, aes(x = BFP), fill = "blue", alpha = 0.5) +
  scale_x_flowjo_biexp() +
  theme_bw() +
  geom_vline(xintercept = 2400, linetype = "dashed", color = "red") +
  ggtitle('OS384 BFP positive cells')


# Save the plot as an SVG file
ggsave("~/Desktop/OS384_FACS_overlap.svg", plot = OS384_FACS_overlap, device = "svg")


percentage_ctrl <- sum(fs_ctrl$BFP > threshold) / nrow(fs_ctrl) * 100
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


# Create a new column to indicate the condition (control or transduced)
gs$condition <- ifelse(gs$sample == "ctrl", "Control", "Transduced")

# Combine the data from control and transduced samples into a single data frame
combined_data <- rbind(
  subset(gs, sample == "ctrl"),
  subset(gs, sample == "LT")
)


ps <- gs_pop_get_count_with_meta(gs)

ps <- ps %>% mutate(percent_of_parent=Count/ParentCount)
ps %>% select(sampleName,well,Population,Count,ParentCount,percent_of_parent) %>% head() %>% kable


### OS052

