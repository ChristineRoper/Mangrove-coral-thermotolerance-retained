library(readr) # for read_csv
library(tibble) # for tibble
library(dplyr) # for summarise and others
library(tidyr) # for drop_na
library(ggplot2)
library(lubridate) # for floor_date

options(readr.show_col_types = FALSE)

AEST <- "Australia/Brisbane"

# 2022-2023 Mangrove HOBO_Temp_pH
logger_data <- read_csv("Temp Feb 22 - Feb 23 Low Isles/2022-2023 Mangrove HOBO_Temp_pH.csv") %>%
    mutate(
        # Parse date
        timestamp = as.POSIXct(Date_Time, format = "%Y-%m-%d %H:%M:%S", tz = AEST),
        temperature = Temp,
        pH = pH,
        # mV=mV,
        site = "Mangrove"
    ) %>%
    # skip before 2022-02-09 11:16:20 - they are suspicious
    filter(timestamp > as.POSIXct("2022-02-09 11:16:20")) %>%
    select(timestamp, temperature, pH, site)

# 2022 Mangrove 2nd logger
logger_data <- read_csv("Temp Feb 22 - Feb 23 Low Isles/2022-2023 Mangrove2 HOBO_Temp_pH.csv") %>%
    mutate(
        # Parse date
        timestamp = as.POSIXct(Date_Time, format = "%d/%m/%Y %H:%M", tz = AEST),
        temperature = Temperature,
        site = "Mangrove"
    ) %>%
    select(timestamp, temperature, site) %>%
    # skip first 24 - they are suspicious
    filter(timestamp > as.POSIXct("2022-02-06 17:00")) %>%
    bind_rows(logger_data)

# 2022-2023 Reef HOBO_Temp_pH
logger_data <- read_csv("Temp Feb 22 - Feb 23 Low Isles/2022-2023 Reef HOBO_Temp_pH.csv") %>%
    mutate(
        # Parse date
        timestamp = as.POSIXct(Date_Time, format = "%Y-%m-%d %H:%M:%S", tz = AEST),
        temperature = Temp,
        pH = pH,
        # mV=mV,
        site = "Reef"
    ) %>%
    # skip before 2022-02-09 12:23:51 - they are suspicious
    filter(timestamp > as.POSIXct("2022-02-09 12:23:51")) %>%
    select(timestamp, temperature, pH, site) %>%
    bind_rows(logger_data)

# 2022-2023 Reef 2nd logger
logger_data <- read_csv("Temp Feb 22 - Feb 23 Low Isles/2022-2023 Reef2 HOBO_Temp_pH.csv", skip = 1, name_repair = "universal") %>%
    mutate(
        # Parse date
        timestamp = as.POSIXct(Date.Time..GMT.11.00, format = "%d/%m/%Y %H:%M", tz = AEST),
        temperature = Temp...C..LGR.S.N..20848679..SEN.S.N..20848679..LBL..Temperature.,
        site = "Reef"
    ) %>%
    drop_na(temperature) %>%
    select(timestamp, temperature, site) %>%
    # skip first 24 - they are suspicious
    filter(timestamp > as.POSIXct("2022-02-06 17:00")) %>%
    bind_rows(logger_data)

# logger_data <- logger_data %>%
#     mutate(
#         # add end of each 8-day period (SST data is grouped by 8 days)
#         eight_days = floor_date(timestamp, "8 days")
#     ) %>%
#     # average each 8-day block as mean_temp
#     group_by(eight_days, site) %>%
#     summarise(mean_temp = mean(temperature, na.rm = TRUE))

# logger_means <- logger_data %>%
#     filter(eight_days > as.POSIXct("2022-02-09 00:00:00"), eight_days <= as.POSIXct("2023-02-09 00:00:00")) %>%
#     group_by(site) %>%
#     summarise(mean = mean(mean_temp), min = min(mean_temp), max = max(mean_temp), min_date = min(eight_days), max_date = max(eight_days))

# coeff <- 5
# offset = -15

cutoff_temp <- 35

temp_data <- logger_data %>%
    filter(temperature < cutoff_temp) %>%
    # Call temperature "value" and add a column variable = "Temperature (ºC)"
    mutate(value = temperature, variable = "Temperature (ºC)")

ph_data <- logger_data %>%
    filter(!is.na(pH), temperature < cutoff_temp) %>%
    mutate(value = pH, variable = "pH")

all_data =  bind_rows(temp_data, ph_data)

means <- all_data %>%
    group_by(site, variable) %>%
    summarise(mean = mean(value), min = min(value), max = max(value), min_date = min(timestamp), max_date = max(timestamp))

means %>% write.csv("logger_stats.csv")

plot <- all_data %>%
    # average each time_period block
    # mutate(time_period = floor_date(timestamp, "7 day")) %>%
    # group_by(time_period, site, variable) %>%
    # summarise(mean_value = mean(value, na.rm = TRUE)) %>%
    ggplot(aes(x = timestamp, y = value, color = site)) +
    # Means
    geom_hline(
        data = means,
        aes(yintercept = mean),
        color='#000000',
        linewidth = 0.3,
        linetype = "dashed"
    ) +
    geom_line(linewidth = 0.1) +
    facet_grid(
        cols = vars(site),
        rows = vars(variable),
        scales = "free_y"
    ) +
    ylab(element_blank()) +
    xlab("Date") +
    # xlim(as.POSIXct("2021-01-01 00:00:00"), as.POSIXct("2023-03-05 00:00:00")) +
    scale_color_manual(values = c("#f56363", "#5c8ef1")) +
    # Christine Theme
    theme(
        legend.position = "none",
        panel.background = element_blank(),
        panel.border = element_rect(fill = NA),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        strip.background = element_blank(),
        axis.text.x = element_text(colour = "black"),
        axis.text.y = element_text(colour = "black"),
        axis.ticks = element_line(colour = "black"),
        plot.margin = unit(c(0.1, 0.1, 0.1, 0.1), "line")
    )

# print(plot)

ggsave(
    "plot.png",
    plot,
    width = 8, height = 5
)
