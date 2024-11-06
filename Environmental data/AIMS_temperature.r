library(readr) # for read_csv
library(tibble) # for tibble
library(dplyr) # for summarise and others
library(tidyr) # for drop_na
library(ggplot2)
library(lubridate) # for floor_date

options(readr.show_col_types = FALSE)

AEST <- "Australia/Brisbane"


aims_data <- read_csv("~/Documents/PhD/Chapter 4/Data/R analysis/Environmental data/Low Isles LOWFL1 Water Temperature @1.0m from 2021-01-01 to 2023-09-06.csv") %>%
    mutate(
        # Parse date
        timestamp = as.POSIXct(time, format = "%Y-%m-%d %H:%M:%S"),
        temperature = qc_val,
        pH = 0,
        site = "Reef"
    ) %>%
    # We only want Feb 2021 to feb 2023
    filter(timestamp >= as.POSIXct("2021-02-01 00:00:00", tz="UTC"), timestamp < as.POSIXct("2023-02-01 00:00:00", tz="UTC")) %>%
    mutate(year = if_else(timestamp < as.POSIXct("2022-02-01 00:00:00", tz="UTC"), 2021, 2022)) %>%
    select(timestamp, temperature, year)

# summary

aims_summary = aims_data %>% group_by(year) %>% summarise(
    min = min(temperature), 
    max = max(temperature),
    mean = mean(temperature),
    min_date = min(timestamp),
    max_date = max(timestamp)
) 

# Plot

timestamp = aims_data$timestamp[400]
timestamp
paste(as.Date(timestamp), hour(timestamp))

hourly = aims_data %>%
    mutate(hour = paste(as.Date(timestamp), hour(timestamp))) %>%
    group_by(hour) %>%
    summarise(
        temperature = mean(temperature),
        timestamp = max(timestamp),
        year = max(year)
    )

plot <- hourly %>%
    ggplot(aes(x = timestamp, y = temperature, color = as.factor(year))) +
    # Means
    geom_hline(
        data = aims_summary,
        aes(yintercept = mean,  color = as.factor(year)),
        # linewidth = 0.3,
        linetype = "dashed"
    ) +
    geom_line(linewidth = 0.1) +
   
    ylab("Temperature (ºC)") +
    xlab("Date") +
    # xlim(as.POSIXct("2021-01-01 00:00:00"), as.POSIXct("2023-03-05 00:00:00")) +
    scale_color_manual(values = c("black", "black")) +
    facet_grid(cols=vars(year), scales="free_x")+
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

plot

ggsave(
    "~/Documents/PhD/Chapter 4/Data/R analysis/Environmental data/plot_grey.png",
    plot,
    width = 8, height = 5
)
