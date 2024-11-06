# original script by Daniel J Barshis
# modified by Christian R Voolstra and Daniel J Barshis
# Adapted by Christine Roper and Jake Crosby

# input: PAM_Data_tidy.csv
# output: ED50.csv & plots

library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(tibble)
library(emmeans)
library(drc)
source("../plotTheme.r")

setwd("~/Documents/PhD/Chapter 4/FINAL DATA AND CODE/Fv-Fm")


#### ReadIn / QAQC ####
pamdata <- read_csv(
  "PAM_Data_tidy.csv"
) %>%
  # Make columns into numeric and factors
  mutate(FVFM = as.numeric(FVFM), across(!c(FVFM, Temperature), as.factor)) %>% # anything that isn't FVFM or temperature is a factor

  # create a new column where FVFM values are converted to decimals
  mutate(PAM = FVFM / 1000) %>%
  # Add a treatment column
  unite(Treatment, NativeHabitat, ExperiemntalSite, sep = "-", remove = FALSE) %>%
  # Add a TreatmentYear column
  unite(TreatmentYear, Treatment, Year, sep = "-", remove = FALSE) %>%
  # We are only dealing with Pocillopora acuta
  filter(Species == "Pocillopora acuta")

# gives descriptive stats
aggregate(PAM ~ Treatment + Temperature + Year, data = pamdata, summary)
# Checking sample sizes
aggregate(PAM ~ Treatment + Temperature + Year, data = pamdata, length)

#### DRC Curve Fitting ####

# 2022

DRCpam <- drm(PAM ~ Temperature,
  data = filter(pamdata, Timepoint == "T1", Year == 2022),
  curveid = Treatment,
  fct = LL.3(names = c("hill", "max", "ed50")) # log logistic with 3 parameters
)
summary(DRCpam)
compParm(DRCpam, "ed50") # statistical test comparing ED50 between groups
compParm(DRCpam, "ed50", "-")
plot(DRCpam) # gives response curve plot
points(pamdata$Temperature, pamdata$PAM) # adds data points to plot
ED(DRCpam, c(50))[, 1] # prints calculated ED50 value


# 2023

DRCpam <- drm(PAM ~ Temperature,
  data = filter(pamdata, Timepoint == "T1", Year == 2023),
  curveid = Treatment,
  fct = LL.3(names = c("hill", "max", "ed50"))
)
summary(DRCpam)
compParm(DRCpam, "ed50") # parameter ratios are compared (default)
compParm(DRCpam, "ed50", "-") # parameter differences are compared
plot(DRCpam)
points(pamdata$Temperature, pamdata$PAM)
ED(DRCpam, c(50))[, 1]

# Both years (but not wild)

DRCpam <- drm(PAM ~ Temperature,
  data = filter(pamdata, Timepoint == "T1", NativeHabitat != "wild"), # all habitats which are not wild
  curveid = TreatmentYear,
  fct = LL.3(names = c("hill", "max", "ed50"))
)
summary(DRCpam)
compParm(DRCpam, "ed50")
compParm(DRCpam, "ed50", "-")
plot(DRCpam)
points(pamdata$Temperature, pamdata$PAM)
ED(DRCpam, c(50))[, 1]


# All samples (compares all individuals to one another - not grouped by habitat etc)

DRCpam <- drm(PAM ~ Temperature,
  data = filter(pamdata, Timepoint == "T1"),
  curveid = PamID,
  fct = LL.3(names = c("hill", "max", "ed50"))
)
summary(DRCpam)
compParm(DRCpam, "ed50")
compParm(DRCpam, "ed50", "-")

# plot just to QC/visualise data
# plot(DRCpam)
# p <- filter(pamdata, Timepoint == "T1", PamID == 1)
# points(p$Temperature, p$PAM)

# Create ED50.csv
# Get ED50 values for each sample
all_ed50s <- ED(DRCpam, c(50)) %>%
  # Make it a dataframe
  as.data.frame() %>%
  # make a column called rowname from the row names so we can work with them
  rownames_to_column("rowname") %>%
  # split the rowname from e:XX:50 format to 3 cols where we want to keep the XX (PamID)
  separate_wider_delim(rowname, ":", names = c("e", "PamID", "ed")) %>%
  # Create new column called ED50 using the 'Estimate' column values
  mutate(ED50 = Estimate) %>%
  # Drop the columns we don't want ('Estimate' col, and the useless bits we split off e:XX:50 rownames)
  dplyr::select(!c(e, ed, Estimate)) %>%
  # Add initial FV/FM from coefficients
  cbind(
    coef(DRCpam) %>%
      as.data.frame() %>%
      rename(initial = ".") %>%
      filter(startsWith(rownames(.), "max:")) # "max" is initial
  ) %>%
  # Add in Year, Treatment, TreatmentYear based on "PamID"
  left_join(dplyr::select(pamdata, RnaID, PamID, Year, Treatment, TreatmentYear), "PamID", multiple = "first")
# Write to file
write.csv(all_ed50s, "ED50.csv")

### Check for differences in initial (initial FV/FM) ###

# For 2022, by treatment
t.test(initial ~ Treatment,  data=filter(all_ed50s, Year==2022))

### Plot ###

# renaming group names to what we want them to be
pamdata <- mutate(pamdata, Treatment = recode_factor(TreatmentYear,
  "mangrove-reef-2022" = "Mangrove 2022",
  "mangrove-reef-2023" = "Mangrove → reef 2023",
  "reef-reef-2022" = "Reef 2022",
  "reef-reef-2023" = "Reef → reef 2023",
  "wild-mangrove-2023" = "Wild mangrove 2023",
  "wild-reef-2023" = "Wild reef 2023"
))

plotx <- "2022 vs 2023 II" # change this based on which plot you want from the list below

if (plotx == "all") {
  # All 6 groups
  plot_pam <- filter(pamdata, Timepoint == "T1")
} else if (plotx == "2022") {
  # 2022
  plot_pam <- filter(pamdata, Timepoint == "T1", Year == 2022)
} else if (plotx == "2023") {
  # 2023
  plot_pam <- filter(pamdata, Timepoint == "T1", Year == 2023)
} else if (plotx == "2022 vs 2023") {
  # 2022 vs 2023
  plot_pam <- filter(pamdata, Timepoint == "T1", NativeHabitat != "wild")
} else if (plotx == "2022 vs 2023 II") {
  # 2022 vs 2023 II
  targets <- c(
    "Mangrove 2022",
    "Reef 2022",
    "Wild mangrove 2023",
    "Wild reef 2023"
  )
  plot_pam <- filter(pamdata, Treatment %in% targets)
}

# # test code - don't need
# targets <- c(
#   "Mangrove 2022",
#   "Wild mangrove 2023",
#   "Mangrove → reef 2023"
# )
# plot_pam <- filter(pamdata, Treatment %in% targets)


DRCpam <- drm(PAM ~ Temperature,
  data = plot_pam, # make plot based on target from above
  curveid = Treatment,
  fct = LL.3(names = c("hill", "max", "ed50"))
)

temperatures <- seq(25, 45, 0.5) # making a sequence from 25 to 45 in steps of 0.5
groups <- levels(droplevels(plot_pam$Treatment))
plot_data <- data.frame(
  Treatment = rep(groups, each = length(temperatures)), # repeat the group for the sequence of temperatures
  Temperature = rep(temperatures, times = length(groups)) #
)
plot_data$FVFM <- predict(DRCpam, newdata = plot_data)
plot_data$Year <- as.factor(ifelse(endsWith(plot_data$Treatment, "2022"), 2022, 2023)) # created for plot below but no longer needed

ED50s <- ED(DRCpam, c(50)) %>% # this gets the ED50 values to provide below for the dashed lines
  # Make it a dataframe
  as.data.frame() %>%
  # make a column called rowname from the row names so we can work with them
  rownames_to_column("rowname") %>%
  # split the rowname from e:XX:50 format to 3 cols where we want to keep the XX (PamID)
  separate_wider_delim(rowname, ":", names = c("e", "Treatment", "ed")) %>%
  # Rename Estimate to ED50
  mutate(ED50 = Estimate) %>%
  # Drop the cols we don't want (old Estimate col, and the useless bits we split off e:XX:50 rownames)
  dplyr::select(!c(e, ed, Estimate))
ED50s$Year <- as.factor(ifelse(endsWith(ED50s$Treatment, "2022"), 2022, 2023))

# creates data points of mean FVFM of each group for plot below
means <- aggregate(PAM ~ Temperature * Treatment, data = filter(pamdata, Treatment %in% targets), mean)

# plot data and model fit
ggplot(plot_data, aes(Temperature, FVFM, col = Treatment)) +
  labs(x = "Temperature (ºC)", y = "Fv/Fm") +
  # Creates ED50 dashed vertical lines
  geom_vline(aes(xintercept = ED50, col = Treatment), ED50s, linetype = "dashed", show.legend = FALSE) +
  # Curve
  geom_line(lwd = 1) +
  # Points
  geom_point(aes(y = PAM), data = means) +
  scale_x_continuous(n.breaks = 6) +
  coord_cartesian(
    xlim = c(31, 41)
    #  ylim = c(0, max(sample_data[[target]]))
  ) +
  # facet_grid(cols=vars(Year)) +
  scale_colour_manual(
    values = c(
      "#f5a2ae",
      "#28e0ced2",
      "#f48703ae",
      "#30be4fa2",
      "#51b9df",
      "#ea2424ca"
    ),
    breaks = c( # setting the names in the legend based on colours
      "Mangrove 2022",
      "Reef 2022",
      "Mangrove → reef 2023",
      "Reef → reef 2023",
      "Wild reef 2023",
      "Wild mangrove 2023"
    )
  ) +
  theme(legend.key = element_blank()) +
  plotTheme

ggsave(paste0("plots/ED50_plot_", plotx, "_points.png"), width = 6, height = 4) # changes plot name based on 'plotx' selection above


# Boxplot horizontal
boxdata <- all_ed50s
# changes the order of the groups
boxdata$TreatmentYear <- ordered(boxdata$TreatmentYear, levels = c(
  "reef-reef-2022",
  "mangrove-reef-2022",
  "reef-reef-2023",
  "wild-reef-2023",
  "mangrove-reef-2023",
  "wild-mangrove-2023"
))
boxdata$TreatmentYear <- recode_factor(boxdata$TreatmentYear,
  # changing the name for the groups in the plot
  "reef-reef-2022" = "wild reef (2022)",
  "mangrove-reef-2022" = "wild mangrove (2022)",
  "reef-reef-2023" = "reef to reef",
  "wild-reef-2023" = "wild reef (2023)",
  "mangrove-reef-2023" = "mangrove to reef",
  "wild-mangrove-2023" = "wild mangrove (2023)"
)

boxdata %>%
  ggplot(aes(ED50, TreatmentYear, fill = TreatmentYear)) +
  labs(y = "Treatment", x = "ED50 (ºC)") +
  geom_boxplot() +
  facet_grid(rows = vars(Year), scales = "free_y", space = "free_y", switch = "y") + # free_y allows scales to be different
  # geom_hline(aes(yintercept = ED50, col = Treatment), ED50s, linetype = "dashed") +
  scale_fill_manual(
    values = c(
      "#4EE5D7",
      "#F5A2AE",
      "#6FCA83",
      "#51b9df",
      "#F7AD53",
      "#EE5252"
    )
  ) +
  scale_y_discrete(position = "right") +
  theme(legend.position = "none") +
  plotTheme

ggsave("plots/ED50_boxplot.png", width = 8, height = 4)


# Boxplot vertical - not using
# plot data and model fit
all_ed50s %>%
  ggplot(aes(Treatment, ED50, fill = Treatment)) +
  labs(x = "Treatment", y = "Fv/Fm") +
  # ED50
  geom_boxplot() +
  facet_grid(cols = vars(Year), scales = "free_x", space = "free_x") +
  # geom_hline(aes(yintercept = ED50, col = Treatment), ED50s, linetype = "dashed") +
  # facet_grid(cols=vars(Year)) +
  scale_colour_manual(
    values = c(
      "#f5a2ae",
      "#28e0ced2",
      "#f48703ae",
      "#30be4fa2",
      "#51b9df",
      "#ea2424ca"
    ),
    breaks = c(
      "Mangrove 2022",
      "Reef 2022",
      "Mangrove → reef 2023",
      "Reef → reef 2023",
      "Wild reef 2023",
      "Wild mangrove 2023"
    )
  ) +
  theme(legend.position = "none") +
  plotTheme

ggsave("plots/ED50_boxplot_vertical.png", width = 8, height = 10)
