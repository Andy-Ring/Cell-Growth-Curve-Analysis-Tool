# Cell Growth Curve Analysis Tool

## Overview

This Shiny app provides an interactive interface for analyzing cell growth data, fitting regression models, performing statistical tests, and visualizing growth trends. Users can upload experimental data, run statistical tests, and generate high-resolution plots.

## Features

- **Data Input:** Upload growth data in CSV format.
- **Template Download:** Download a preformatted CSV template for data entry.
- **Growth Curve Analysis:** Fit logistic or exponential growth models.
- **Statistical Testing:** Perform T-tests, ANOVA, Kruskal-Wallis, and more.
- **Customizable Plots:** Adjust plot titles, axis labels, and statistical comparisons.
- **Downloadable Results:** Save high-resolution plots and analysis summaries.

## Web Acess
This app can be acessed from https://01951577-cbc7-9ef2-1a78-48f59e4541d7.share.connect.posit.cloud/

## Running the App Locally
### Installation
To run this Shiny app, ensure you have R and the required packages installed. You can install the necessary dependencies using:

```r
install.packages(c("shiny", "ggplot2", "dplyr", "car", "stats", "ggpubr", "multcomp", "bslib", "thematic"))
```

### Running the App
1. Clone this repository or download the script files.
2. Open R or RStudio and navigate to the app directory.
3. Run the following command:
   ```r
   shiny::runApp()
   ```

## User Guide

### 1. Data Input
- Download a **CSV template** for data entry.
- Upload a CSV file containing **time, growth, and group** information.
- Ensure the uploaded data follows the required format.

### 2. Growth Curve Analysis
- The app tries to fit a **logistic growth model**, if that fails, it will default to an **exponential growth model**.
- Displays **growth rate and doubling time** calculations.
- Provides regression equations for further analysis.

### 3. Statistical Testing
- Click **Recommend Best Test** to determine the most suitable statistical test. This function runs tests for normality and equal variance which the recommendation is based on.
- Choose from **T-test, ANOVA, Welch’s ANOVA, Kruskal-Wallis, or Pairwise Wilcoxon**.
- View test results in the **Statistical Test Results** section.

### 4. Visualization Options

- Customize the **plot title, x-axis, and y-axis labels**.
- Generate and download a **growth curve plot**.
- Display statistical comparisons directly on the plots.

## Downloads

- **Growth Curve Plot** (PNG)
- **Statistical Test Results** (CSV)
- **Regression Equations** (CSV)
- **Analysis Summary** (TXT)

## Credits

- **Author:** Andy Ring
- **Version:** 1.0.4
- **Date:** February 14th, 2025

## License

This project is licensed under the MIT License. Feel free to use and modify it for your research and analysis needs.

