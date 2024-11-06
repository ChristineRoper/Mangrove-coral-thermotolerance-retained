# Setup
#########
library(car) # For Shapiro test
library(ARTool) # for ART ANOVA

setwd("Methylated DNA")

input <- "MethylatedDNA_raw.csv"
DNA_data <- read.csv(input, header = TRUE)

# make factor
DNA_data$Sample_group <- as.factor(DNA_data$Sample_group)
DNA_data$Group <- as.factor(DNA_data$Group)
DNA_data$Temperature <- as.factor(DNA_data$Temperature)

# Check for normality
##########################

# DNA_data$Methyl_percent_log <- log(DNA_data$Methyl_percent)
shapiro.test(DNA_data$Methyl_percent) # did not pass
# qqnorm(DNA_data$Methyl_percent, pch = 1, frame = FALSE)
# qqline(DNA_data$Methyl_percent, col = "steelblue", lwd = 2)

leveneTest(Methyl_percent ~ Sample_group, DNA_data) # passed

# show means and standard error
stdErr <- function(x) sd(x) / sqrt(length(x))
DNA_data %>%
  group_by(Group, Temperature) %>%
  summarise(mean_methyl_percent = mean(Methyl_percent), SE = stdErr(Methyl_percent)) %>%
  View()


# Non-parametric alternative to two-way Anova
##############################################

Art_anova <- art(Methyl_percent ~ Group + Temperature + Group:Temperature, data = DNA_data)
anova(Art_anova)
# post-hoc comparisons
art.con(Art_anova, "Group:Temperature", adjust = "none")


### Plots ###

library(ggplot2)
source('../plotTheme.r')

# Bar plot
DNA_data %>%
  group_by(Group, Temperature, Sample_group) %>%
  summarise(mean_methyl_percent = mean(Methyl_percent), SE = stdErr(Methyl_percent)) %>%
  ggplot(aes(x = Sample_group, y = mean_methyl_percent, fill=Temperature )) +
  geom_bar(stat = "identity") +
  geom_errorbar(aes(ymin = mean_methyl_percent - SE, ymax = mean_methyl_percent + SE), width = 0.1,  colour = "gray40") +
  scale_fill_manual(values=c("#88C673", "#FCDC6A","#FF9B06")) +
  labs(x = "Group", y = "DNA Methylation", fill = "Temperature (ºC)") + # Label axes
  #theme(axis.text.x = element_text(angle = 20, vjust = 0.7, hjust=0.65)) +
  scale_x_discrete(labels = rep(c('mangrove\nto reef','reef\nto reef','wild\nmangrove'), 3))  + # Change the x-axis labels
  plotTheme
ggsave("Methylated DNA bargraph.png", width = 7, height = 5)


# Boxplot - did not use
ggplot(DNA_data, aes(x = Sample_group, y = Methyl_percent)) +
  geom_point() + # Add points
  geom_boxplot() + # Add boxplot for better visualization of distribution
  labs(x = "Sample Group", y = "Methyl Percent") + # Label axes
  ggtitle("Methyl Percent by Sample Group") # Add title

