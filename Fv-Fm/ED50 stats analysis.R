
library(readr)
library(dplyr)
library(emmeans)
library(dunn.test)
library(FSA)
library(ARTool)

# input: ED50.csv
# output: provides stats results, but not output to files

setwd("Fv-Fm")

#### read csv ####
ED50_data <- read_csv("ED50.csv") %>% 
  # Split treatment and make new columns for origin and destination
  separate(Treatment, c("Origin", "Destination"), remove=FALSE) %>%
  # make all factors
  mutate(across(c(Treatment, TreatmentYear, PamID, Year, RnaID), as.factor))

ED50_data %>% group_by(Year, Treatment) %>% summarise(number=length(ED50)) # QC check for sample sizes/group names


#### ED50 comparisons ####

#filter based on groups of interest
ED50_2022 = filter(ED50_data, Year==2022)
ED50_2023 = filter(ED50_data, Year==2023)


#### Compare the two groups from 2022 ####

# check assumptions for each treatment
# aggregate ED50's based on their treatment, Fun = a function that returns the shapiro test p-value
aggregate(ED50 ~ Treatment, ED50_2022, FUN = function(x) shapiro.test(x)$p.value) #passed
bartlett.test(ED50 ~ Treatment, ED50_2022) #passed

#perform unpaired t-test
t.test(ED50 ~ Treatment, ED50_2022, var.equal = TRUE) # var.equal = TRUE because data passed Bartlett test


#### Compare all 4 groups in 2023 #####

#check assumptions
# aggregate ED50's based on their treatment, Fun = a function that returns the shapiro test p-value
aggregate(ED50 ~ Treatment, ED50_2023, FUN = function(x) shapiro.test(x)$p.value) #passed
bartlett.test(ED50 ~ Treatment, ED50_2023) #did not pass

#perform Kruskal-Wallis tests (non-parametric alternative to one-way anova since data did not pass Barlett test)
kruskal.test(ED50 ~ Treatment, data = ED50_2023)
dunnTest(ED50 ~ Treatment, data = ED50_2023, method = "bonferroni") # recommended post-hoc comparison with bonferroni correction


#### Compare 2022 vs 2023 wild corals (e.g. mangrove 2022 vs wild mangrove 2023) ####

ED50_wildcontrols = filter(
  ED50_data,
  TreatmentYear %in% c("mangrove-reef-2022", "reef-reef-2022", "wild-mangrove-2023", "wild-reef-2023" )
)

#check assumptions (as above)
aggregate(ED50 ~ Treatment*Year, ED50_wildcontrols, FUN = function(x) shapiro.test(x)$p.value) #passed
bartlett.test(ED50 ~ TreatmentYear, ED50_wildcontrols) #passed

# Check how many in each group
ED50_wildcontrols %>% 
  group_by(Treatment, Year) %>% 
  summarise(number=length(ED50),mean=mean(ED50))

#T-tests
t.test(ED50 ~ Year, data = filter(ED50_wildcontrols, Treatment %in% c("mangrove-reef", "wild-mangrove")), var.equal = TRUE) # var.equal = TRUE because data passed Bartlett test
t.test(ED50 ~ Year, data = filter(ED50_wildcontrols, Treatment %in% c("reef-reef", "wild-reef")), var.equal = TRUE)

ggplot(ED50_wildcontrols, aes(y = ED50, x = TreatmentYear, col = TreatmentYear)) + geom_boxplot()


#### Pairwise comparison of the same 2022 vs 2023 individuals (e.g. mangrove 2022 vs mangrove-reef 2023) ######
library(rstatix)

#filter out any corals that were not replicated in both 2022 and 2023
ED50_2022vs2023 = filter(
    ED50_data,
    Treatment %in% c("mangrove-reef", "reef-reef"),
   !PamID %in% c(24, 26, 27, "T", 4, 7, "O", "P", "W", "X") # keep all samples not in this list
  ) %>%
  mutate(Treatment = droplevels(Treatment), Year = droplevels(Year), RnaID = droplevels(RnaID)) %>% # update factors to only have levels that exist in the data
  # Sort by RnaID. It doesn't change from year to year like PamID does so we can use it to make sure the samples are in the same order each year
  arrange(RnaID)

#check assumptions as above
aggregate(ED50 ~ Treatment*Year, ED50_2022vs2023, FUN = function(x) shapiro.test(x)$p.value) #did not pass for reef-reef 2022
bartlett.test(ED50 ~ TreatmentYear, ED50_2022vs2023) #passed

# Check how many in each group
ED50_2022vs2023 %>% group_by(Treatment, Year) %>% summarise(number=length(ED50))
ED50_2022vs2023 %>% group_by(RnaID) %>% summarise(number=length(ED50)) # check data to ensure each sample appears in both years

#paired T-test
t.test(ED50 ~ Year, data = filter(ED50_2022vs2023, Origin=='mangrove'), paired = TRUE)
wilcox.test(ED50 ~ Year, data = filter(ED50_2022vs2023, Origin=='reef'), paired = TRUE) # non-parametric test since shapiro wilk test failed

# quick visualisation
#ggplot(ED50_2022vs2023, aes(y = ED50, x = TreatmentYear, col = TreatmentYear)) + geom_boxplot()


### Boxplots (note: not sytlised) ###

# 2022 boxplot
ggplot(ED50_2022, aes(x = Treatment, y = ED50)) +
  geom_point() + # Add points
  geom_boxplot() + # Add boxplot for better visualization of distribution
  labs(x = "Treatment group", y = "ED50") # Label axes

# 2023 boxplot
ggplot(ED50_2023, aes(x = Treatment, y = ED50)) +
  geom_point() + # Add points
  geom_boxplot() + # Add boxplot for better visualization of distribution
  labs(x = "Treatment group", y = "ED50") # Label axes

# 2022 vs 2023 boxplot
ggplot(ED50_2022vs2023, aes(x = TreatmentYear, y = ED50)) +
  geom_point() + # Add points
  geom_boxplot() + # Add boxplot for better visualization of distribution
  labs(x = "Treatment group", y = "ED50") # Label axes


### Show descriptive stats ###
stdErr <- function(x) sd(x) / sqrt(length(x)) # create standard error function

ED50_2022 %>%
  group_by(Treatment) %>%
  summarise(mean_ED50 = mean(ED50), SE = stdErr(ED50)) %>%
  View()

ED50_2023 %>%
  group_by(Treatment) %>%
  summarise(mean_ED50 = mean(ED50), SE = stdErr(ED50)) %>%
  View()

ED50_2022vs2023 %>%
  group_by(Treatment) %>%
  summarise(mean_ED50 = mean(ED50), SE = stdErr(ED50)) %>%
  View()

