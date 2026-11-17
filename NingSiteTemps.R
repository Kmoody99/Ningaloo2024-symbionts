#Ningaloo temperature data for Bruboodjoo and Tantabiddi, NOAA CRW data
#https://coastwatch.pfeg.noaa.gov/erddap/griddap/NOAA_DHW.htmlTable?CRW_BAA%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_BAA_mask%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_BAA_7D_MAX%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_BAA_7D_MAX_mask%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_DHW%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_DHW_mask%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_HOTSPOT%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_HOTSPOT_mask%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_SEAICE%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_SST%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_SSTANOMALY%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D,CRW_SSTANOMALY_mask%5B(2004-11-05T12:00:00Z):1:(2024-11-05T12:00:00Z)%5D%5B(89.975):1:(-89.975)%5D%5B(-179.975):1:(179.975)%5D
#
library(tidyverse)
library(ggpubr)
library(ggrepel)

# Load SST data file
sst_data <- read.csv("/Users/laurenhowe-kerr/Dropbox/Lauren Howe-Kerr/MERL spawning planning 2024/MS- overview Fig1/20yr_SST_CRW.csv")

# Check the dimensions of sst_data
print(dim(sst_data))

# Check for missing data or NaN values
missing_data <- is.na(sst_data)

# Print the number of missing data and NaN values
print(sum(missing_data))

sst_data$site <- as.factor(sst_data$site)

# Ensure 'time' is a datetime column (if it's not already)
sst_data$time <- as.POSIXct(sst_data$time)

####################################
############ Past average temperature leading up to sampling time in April
####################################

# Calculate day of year starting from May 1st as Day 1
df <- sst_data %>%
  mutate(
    date = as.Date(time),                          # Extract the date
    day_of_year = yday(time),                      # Standard day of year (1-365 or 366)
    day_of_year_reordered = if_else(
      month(time) >= 5,                            # If the month is May or later
      yday(time) - yday(as.Date(paste0(year(time), "-05-01"))),  # Subtract from May 1st of the same year
      yday(time) + (if_else(leap_year(year(time)), 366, 365) - yday(as.Date(paste0(year(time), "-05-01"))))  # Subtract from May 1st of the previous year
    )
  )

df_n <- df %>%
  filter(
    time >= as.Date("2004-05-01") & time <= as.Date("2023-04-30")
  )


# Calculate the mean temperature and standard deviation per day of the year for each site
df_stats <- df_n %>%
  group_by(site, day_of_year_reordered) %>%
  summarise(
    mean_temperature = mean(CRW_SST, na.rm = TRUE),
    sd_temperature = sd(CRW_SST, na.rm = TRUE),  # Standard deviation
    .groups = "drop"
  )


# Filter the data for 2023-2024 period (May 1, 2023 - April 30, 2024) for comparison for each site
df_2023_2024_B <- df %>%
  filter(site == "Bruboodjoo") %>% filter(
    (year(time) == 2023 & month(time) >= 5) | year(time) == 2024 | 
      (year(time) == 2023 & month(time) < 5)
  ) %>%
  group_by(day_of_year_reordered) %>%
  summarise(
    temp_2023_2024 = mean(CRW_SST, na.rm = TRUE),  # Mean temperature for the 2023-2024 period
    mean_dhw = mean(CRW_DHW, na.rm = TRUE),  # New variable for DHW
    sd_dhw = sd(CRW_DHW, na.rm = TRUE),  # Standard deviation for DHW
    .groups = "drop"
  )

df_2023_2024_T <- df %>%
  filter(site == "Tantabiddi") %>% filter(
            (year(time) == 2023 & month(time) >= 5) | year(time) == 2024 | 
            (year(time) == 2023 & month(time) < 5)
  ) %>%
  group_by(day_of_year_reordered) %>%
  summarise(
    temp_2023_2024 = mean(CRW_SST, na.rm = TRUE),    # Mean temperature for the 2023-2024 period
    mean_dhw = mean(CRW_DHW, na.rm = TRUE),  # New variable for DHW
    sd_dhw = sd(CRW_DHW, na.rm = TRUE),  # Standard deviation for DHW
    .groups = "drop"
  )

# Plot the results with ±1 standard deviation and a black line for the 2023-2024 period (Supp Fig.)
ggplot(df_stats, aes(x = day_of_year_reordered, y = mean_temperature, color = site)) +
  geom_line() +  # Line for mean temperature
  geom_ribbon(aes(ymin = mean_temperature - sd_temperature, 
                  ymax = mean_temperature + sd_temperature, 
                  fill = site), 
              alpha = 0.2, color = NA) +  # Shaded ribbon for ±1 SD
  labs(
    x = "Day of Year (May 1st to April 30th)",
    y = "Mean Temperature (°C)"
  ) +
  theme_classic() +
  scale_y_continuous(
    breaks = c(20,22,24,26,28)) +  # Set breaks every 1 unit
  scale_x_continuous(breaks = c(0,61,123,184,244,303)) +
  scale_color_manual(values = c("#4D2D18", "#AC844E")) + 
  scale_fill_manual(values = c("#4D2D18", "#AC844E")) + 
  geom_line(data = df_2023_2024_T, aes(x = day_of_year_reordered, y = temp_2023_2024), color = "#AC844E", size = 0.5, linetype='dashed') +  # Black line for 2023-2024
  geom_line(data = df_2023_2024_B, aes(x = day_of_year_reordered, y = temp_2023_2024), color = "#4D2D18", size = 0.5, linetype='dashed')  # Black line for 2023-2024


####################################
############ Variability 
####################################

df <- sst_data

##VARIABILITY ANNUALLY - in terms of standard deviation and yearly range
# Filter data for the past 20 years (2004-2024)
df_filtered <- df %>%
  filter(
    time >= as.Date("2004-05-01") & time <= as.Date("2024-04-30")
  )

# Extract the year from the date column
df_filtered$year <- as.integer(format(df_filtered$time, "%Y"))

# Calculate the yearly standard deviation of temperature for each site, as well as mean and max temps, max dhw, and range
df_yearly_variability <- df_filtered %>%
  group_by(site, year) %>%
  summarise(
    yearly_mean = mean(CRW_SST, na.rm = TRUE),  
    yearly_sd = sd(CRW_SST, na.rm = TRUE),  # Calculate standard deviation for each year
    yearly_min = min(CRW_SST, na.rm = TRUE),
    yearly_max = max(CRW_SST, na.rm = TRUE),
    yearly_range = (yearly_max - yearly_min),
    yearly_dhw = mean(CRW_DHW, na.rm = TRUE),
    yearly_maxdhw=max(CRW_DHW,na.rm = TRUE),
    .groups = "drop"
  )

# Calculate the average yearly variability (mean of the standard deviations) for each site
df_avg_variability <- df_yearly_variability %>%
  group_by(site) %>%
  summarise(
    avg_yearly_sd = mean(yearly_sd, na.rm = TRUE),  # Average of the yearly standard deviations
    avg_yearly_range = mean(yearly_range, na.rm = TRUE),
    .groups = "drop"
  )

# View the result
print(df_avg_variability)

# Boxplot of yearly range for each site
ggplot(df_yearly_variability, aes(x = site, y = yearly_range)) +
  geom_point(aes(color = year)) +
  geom_boxplot(alpha = 0.3) +  # Boxplot with fill color
  labs(
    x = "Site",
    y = "Yearly Range"
  ) +
  theme_minimal() +
  scale_color_gradient(low = "darkblue", high = "lightblue") +  # Apply color gradient
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)) 

anova_result <- aov(yearly_sd ~ site, data = df_yearly_variability)
anova_result <- aov(yearly_range ~ site, data = df_yearly_variability)
summary(anova_result)

####################################
############ DHW


# DHW
ggplot(df_yearly_variability, aes(x = site, y = yearly_maxdhw, color = year, group = year)) +
  geom_point() +  # Color points according to 'year'
  geom_line() +  # Add lines connecting points for the same year
  labs(
    x = "Site",
    y = "Yearly Maximum dHW"
  ) +
  theme_classic() +
  #scale_color_gradient(low = "lightgray", high = "black") +  # Apply color gradient
  scale_color_gradient(low = "darkblue", high = "red") +  # Apply color gradient
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if necessary
    legend.position = "right"
  )

label_data <- df_yearly_variability[df_yearly_variability$year %in% c(2011, 2012, 2013, 2021), ]


dhwplot <- ggplot(df_yearly_variability, aes(x = site, y = yearly_maxdhw, color = site, group = year)) +
  geom_line() +  # Add lines connecting points for the same year
  geom_point(size=2) +  # Color points according to 'year'
  theme_classic() +
  labs(
    x = NULL,
    y = "Degree heating weeks (°C-weeks)") +
  #scale_color_gradient(low = "lightgray", high = "black") +  # Apply color gradient
  scale_color_manual(values = c("#4D2D18", "#AC844E")) + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),  # Rotate x-axis labels if necessary
    axis.ticks.x = element_blank(),
    legend.position = "none"  # Remove the legend (since it's redundant in this case)
  ) +
  geom_hline(yintercept = 8, linetype ="dashed", size = 0.75) +
  geom_text_repel(
    data = label_data,
    aes(label = year),
    size = 4,
    direction = "y",        # Only adjust along the y-axis
    nudge_x = -0.3,         # Push labels to the left
    segment.color = NA,     # Remove the pointing line
    hjust = 1,              # Align text to the right since it's to the left of the point
    show.legend = FALSE
  )
dhwplot

# Check normality for one site
shapiro.test(df_yearly_variability$yearly_maxdhw[df_yearly_variability$site == "Bruboodjoo"])

# Check normality for the other site
shapiro.test(df_yearly_variability$yearly_maxdhw[df_yearly_variability$site == "Tantabiddi"])

# Calculate the average yearly DHW for each site
df_avg_variability <- df_yearly_variability %>%
  group_by(site) %>%
  summarise(
    avg_yearly_dhw = mean(yearly_dhw, na.rm = TRUE),  # Average of the yearly standard deviations
    avg_yearly_maxdhw = mean(yearly_maxdhw, na.rm = TRUE),
    sd_yearly_dhw = sd(yearly_dhw, na.rm = TRUE),  # Average of the yearly standard deviations
    sd_yearly_maxdhw = sd(yearly_maxdhw, na.rm = TRUE),
    .groups = "drop"
  )


df_yearly_variability$year <- as.factor(df_yearly_variability$year)

# MAXDHW
# Extract the yearly_maxdhw values for each site
site1_values <- df_yearly_variability$yearly_maxdhw[df_yearly_variability$site == "Bruboodjoo"]
site2_values <- df_yearly_variability$yearly_maxdhw[df_yearly_variability$site == "Tantabiddi"]

# Ensure we have 21 years of data for both sites (check if they match)
length(site1_values)  # Should be 21
length(site2_values)  # Should be 21

# Calculate the number of positive differences
positive_differences <- sum(site1_values > site2_values)
negative_differences <- sum(site1_values < site2_values)

# Perform a binomial test (Sign test)
n <- length(site1_values)  # Number of pairs (21 years)
binom_result <- binom.test(positive_differences, n, p = 0.5, alternative = "two.sided")

# Print the result
print(binom_result)



####################################
############ Maximum Monthly means

df <- sst_data

##VARIABILITY ANNUALLY
# Filter data for the past 20 years (2004-2024)
df_filtered <- df %>%
  filter(
    time >= as.Date("2004-05-01") & time <= as.Date("2024-04-30")
  )


# Create year, month, and site columns
df <- df_filtered %>%
  mutate(
    year = year(time),   # Extract year from date
    month = month(time)  # Extract month from date
  )

# Calculate the monthly means for each site, year, and month
monthly_means <- df %>%
  group_by(site, year, month) %>%
  summarise(monthly_mean = mean(CRW_SST, na.rm = TRUE)) %>%
  ungroup()

# Find the highest monthly mean for each site and year
MMM <- monthly_means %>%
  group_by(site, year) %>%
  summarise(highest_monthly_mean = max(monthly_mean, na.rm = TRUE)) %>%
  ungroup()

# View the result
print(MMM)

#gradient version
ggplot(MMM, aes(x = site, y = highest_monthly_mean, color = year, group = year)) +
  geom_point() +  # Color points according to 'year'
  geom_line() +  # Add lines connecting points for the same year
  theme_classic() +
  #scale_color_gradient(low = "lightgray", high = "black") +  # Apply color gradient
  scale_color_gradient(low = "darkblue", high = "red") +  # Apply color gradient
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),  # Rotate x-axis labels if necessary
    legend.position = "right" 
  )

#site color version
MMMplot <- ggplot(MMM, aes(x = site, y = highest_monthly_mean, color = site, group = year)) +
  geom_line() +  # Add lines connecting points for the same year
  geom_point(size=2) +  # Color points according to 'year'
  theme_classic() +
  labs(
    x = NULL,
    y = "Max. monthly mean temp. (°C)") +
  #scale_color_gradient(low = "lightgray", high = "black") +  # Apply color gradient
  scale_color_manual(values = c("#4D2D18", "#AC844E")) + 
  theme(
    axis.text.y = element_text(size = 14),
    axis.title.y = element_text(size = 14),
    axis.text.x = element_blank(),  # Rotate x-axis labels if necessary
    axis.ticks.x = element_blank(),
    legend.position = "none" # Remove the legend (since it's redundant in this case)
  )

site1_values <- MMM$highest_monthly_mean[MMM$site == "Bruboodjoo"]
site2_values <- MMM$highest_monthly_mean[MMM$site == "Tantabiddi"]
# Check normality for one site
shapiro.test(MMM$highest_monthly_mean[MMM$site == "Bruboodjoo"])

# Check normality for the other site
shapiro.test(MMM$highest_monthly_mean[MMM$site == "Tantabiddi"])

# Perform a paired t-test
t_test_result <- t.test(site1_values, site2_values, paired = TRUE)

# Print the result
print(t_test_result)
sd(site1_values)
mean(site1_values)
sd(site2_values)

mean(site2_values)

combined_plot <- ggarrange(MMMplot, dhwplot, 
                           ncol = 2, nrow = 1)

combined_plot 
#saved as 5x3 PDF
