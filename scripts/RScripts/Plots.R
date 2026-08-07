require(ggplot2)
require(tidyr)
require(dplyr)

# 1.0 Varying coverage thresholds ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data")
data <- read.csv("variedCovThresh.csv")
data$Node <- as.factor(data$Node)
data$CoverageThreshPercent <- as.factor(data$CoverageThreshPercent)

ggplot(data, aes(x = Node, y = nRegions, fill = CoverageThreshPercent)) + 
  geom_col(position = "dodge") +
  theme_minimal() +
  theme(legend.position = "inside",             
        legend.position.inside = c(0.97, 0.95),
        legend.justification = c("right", "top"), 
        legend.background = element_rect(color = "gray")) + 
  labs(y = "Number of Regions") +
  scale_fill_brewer(palette = "Set2", 
                    name = "Coverage Threshold (%)")

# 2.0 Genic fraction ----
## 2.1 Genic vs. intergenic ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data/regionGenomicAnntotations")
genFrac <- read.csv("genicFraction.csv")
ggplot(genFrac, aes(x = Node, y = FractionGenic)) + 
  geom_col(fill = "#74B39B") + 
  theme_minimal() + 
  labs(y = "Fraction overlapping genic annotations")

# Reformat dataframe
genFrac_long <- genFrac %>%
  pivot_longer(
    cols = `GenicRegions`:`IntergenicRegions`, 
    names_to = "regionType",
    values_to = "nRegions"
  )

ggplot(genFrac_long, aes(x = Node, y = nRegions, fill = regionType)) + 
  geom_col(position = "dodge") + 
  theme_minimal() + 
  labs(y = "Number of regions") + 
  scale_fill_brewer(palette = "Set2", 
                    name = "Region Type")

## 2.2 Genic annotations ----
setwd("/Users/haileynapier/Work/VertGenLab/Projects/vertCons/data/regionGenomicAnntotations")
genAnn <- read.table("annTableCat.txt", sep = " ", header = F)
names(genAnn) <- c("nRegions","regionType", "node")

ggplot(genAnn, aes(x = node, y = nRegions, fill = regionType)) + 
  geom_col(position = "stack") + 
  theme_minimal() + 
  scale_fill_brewer(palette = "Dark2") +
  labs(y = "Number of Regions", 
       x = "Node") +
  scale_x_continuous(breaks = 0:17) + 
  theme(panel.grid.minor = element_blank()) 

## 2.3 How many annotations per region (on average)?
annsOverRegions <- genAnn %>%
  group_by(node) %>%
  summarize(total_sum = sum(nRegions))
names(annsOverRegions) <- c("Node", "SumAnns")
totalRegions <- genFrac %>%
  select(Node, TotalRegions)
totalAnnsRegions <- inner_join(totalRegions,annsOverRegions)
totalAnnsRegions$AvAnnsPerRegion <- totalAnnsRegions$SumAnns/totalAnnsRegions$TotalRegions
ggplot(totalAnnsRegions, aes(x = Node, y = AvAnnsPerRegion)) + 
  geom_col(fill = "#74B39B") + 
  theme_minimal() + 
  labs(y = "Average Annotations Per Region")

## 2.4 Percent of regions associated with each annotation
genAnn$Node <- genAnn$node
genAnn <- full_join(genAnn, annsOverRegions)
genAnn$percentAnns <- genAnn$nRegions/genAnn$SumAnns
ggplot(genAnn, aes(x = Node, y = percentAnns, fill = regionType)) + 
  geom_col() +
  theme_minimal() +
  labs(y = "Percent of Total Overlapping Annotations", 
       fill = "Annotation Description") + 
  scale_fill_brewer(palette = "Dark2")
