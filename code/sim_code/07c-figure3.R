# Script for computing R2 stats (e.g., mean % changes, tests between models)
# Author: Nate 
# Date: July 18, 2025

# load packages
library(dplyr)
library(tidyr)
library(readxl)
library(writexl)
library(data.table)

# DEFINE FUNCTIONS ----

# input: result df subsetted by nsnps and h2_total
devdf.pheno.format = function(df) {
  kvals = unique(df$k)
  models = unique(df$Model)
  full_res_list = list()
  
  for (kval in kvals) {
    for (model in models) {
      df_sub = df %>% filter(Model == model, k == kval)
      
      # Calculate mean and sd R2 (handle if no rows)
      mean_r2 = if (nrow(df_sub) > 0) mean(df_sub$R2, na.rm = TRUE) else NA
      sd_r2   = if (nrow(df_sub) > 0) sd(df_sub$R2, na.rm = TRUE) else NA
      
      # Default results (no t-test)
      res_list = list(
        test_model = paste(model, kval, sep = "_"),
        base_model = NA, 
        model_type = model,
        h2_total = NA,
        k = kval, 
        prop_add = NA, 
        prop_dom = NA,
        mean_r2 = mean_r2,
        sd_r2 = sd_r2
      )
      
      record = TRUE
      
      if (model == "ADD") {
        record = FALSE
        # record pheno characteristics
        res_list$h2_total = unique(df_sub$h2_total)
        res_list$prop_add = unique(df_sub$prop_add)
        res_list$prop_dom = unique(df_sub$prop_dom)
      }
      
      if (record) {
        df_ref = df %>% filter(Model == "ADD", k == kval)
        if (nrow(df_ref) > 0 && nrow(df_sub) > 0) {
          res_list$base_model = paste(unique(df_ref$Model), unique(df_ref$k),  sep = "_")
          res_list$h2_total = unique(df_sub$h2_total)
          res_list$prop_add = unique(df_sub$prop_add)
          res_list$prop_dom = unique(df_sub$prop_dom)
        }
      }
      
      outname = paste(model, kval, sep = "_")
      full_res_list[[outname]] = res_list
    }
  }
  return(full_res_list)
}

# runs for all heritability values for devsims
devdf.within.pheno.results = function(df, h2_vals = c(0.5, 0.2, 0.1)) {
  result_dfs = lapply(h2_vals, function(h2) {
    # Filter the dataset
    df_sub = df %>% filter(h2_total == h2)
    
    # Run the function
    res_list = devdf.pheno.format(df_sub)
    
    # Convert result list to a single data frame
    df_result = do.call(rbind, lapply(res_list, as.data.frame))
    rownames(df_result) = NULL
    
    # re-order columns
    df_result = df_result[, c("test_model", "base_model", "model_type", "k", "h2_total", "prop_add", 
                              "prop_dom", "mean_r2", "sd_r2")]
    
    return(df_result)
  })
  
  # Combine all results
  final_df = do.call(rbind, result_dfs)
  
  # sort
  final_df = final_df %>%
    mutate(
      h2_total = as.numeric(as.character(h2_total)),  # factor-safe numeric
      k        = as.numeric(k)
    ) %>%
    arrange(h2_total, desc(k), model_type)
  
  return(final_df)
}

# plot it from table 5
plot_devdf_summary_r2 <- function(df_summary, h2 = 0.5, title = NULL, ylims = NULL, brks = NULL, title_size = 15, tick_size = 15) {
  
  model_colors = c(
    "ADD" = "#584053",  
    "DOM" = "#8DC6BF",  
    "XGB" = "#FCBC66", 
    "DNN" = "#F97B4F"   
  )
  
  # subset 
  df_summary = df_summary %>%
    filter(h2_total == h2)
  
  # orders
  model_order <- c("ADD", "DOM", "XGB", "DNN")
  dom_order <- sort(unique(df_summary$k), decreasing = TRUE)
  
  # factor
  df_plot = df_summary %>%
    mutate(
      model_type = factor(model_type, levels = model_order),
      k = factor(k, levels = dom_order)
    )
  
  # Create alternating shading using annotate (must use numeric indices of factor levels)
  rect_layers <- lapply(seq_along(dom_order), function(i) {
    if (i %% 2 == 1) {
      annotate(
        "rect",
        xmin = i - 0.5, xmax = i + 0.5,
        ymin = -Inf, ymax = Inf,
        fill = "grey95", alpha = 0.5
      )
    }
  })
  
  # set breaks and lines
  y_scale <- if (!is.null(brks)) {
    scale_y_continuous(
      breaks = brks,
      labels = scales::number_format(accuracy = 0.01)
    )
  } else {
    scale_y_continuous()
  }
  
  pos <- position_dodge(width = 0.75)
  
  ggplot(df_plot, aes(x = k, y = mean_r2, color = model_type)) +
    rect_layers + 
    # geom_line(aes(group = model_type), position = pos, linetype = "dotted", linewidth = 0.6) +
    geom_errorbar(
      aes(
        ymin = mean_r2 - sd_r2,
        ymax = mean_r2 + sd_r2,
        group = model_type     # <- this restores dodging
      ),
      position = position_dodge(width = 0.75),
      width = 0.2,
      linetype = "dashed",
      color = "black"
    ) + 
    geom_point(
      aes(fill = model_type),            
      color = "black",                    
      shape = 21,
      position = position_dodge(width = 0.75),
      size = 2.5,
      stroke = 0.8
    ) +
    scale_fill_manual(values = model_colors) + 
    labs(
      x = expression(italic(k)),
      y = expression(R^2),
      fill = "Model",
      title = title
    ) + 
    y_scale + 
    theme_minimal(base_size = 14) +
    theme(
      panel.background = element_blank(),     # removes panel (inside plot area) background
      plot.background = element_blank(),      # removes entire plot canvas background
      panel.grid.major.x = element_blank(),
      panel.grid.minor = element_blank(),
      panel.grid.major.y = element_line(color = "grey80", linetype = "dashed"),
      axis.line = element_line(color = "black", linewidth = 0.6),
      axis.ticks = element_line(color = "black"),
      # ---- add some text size args here babbyyyy
      axis.title.x = element_text(size = title_size),
      axis.title.y = element_text(size = title_size),
      axis.text.x  = element_text(angle = 45, hjust = 1, size = tick_size),
      axis.text.y  = element_text(size = tick_size),
      # ---- continue with regular programming
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      legend.position = "right", 
      legend.title = element_text(size = 15), 
      legend.text  = element_text(size = 15)
    ) +
    (if (!is.null(ylims)) coord_cartesian(ylim = ylims) else NULL)
}

# END FUNCTION DEFINITION ----

# load data
devdf = read_excel("/Users/nyb/phd_code/dl-prs-v2/output/man_sims/figs/june2025/dev_sims_fig_data_aug2025.xlsx", sheet = "Sheet1")

# add cv split numbers to the data frame
devdf = devdf %>% 
  rename(nsnps = Nsnps) %>%
  group_by(Model, Phenotype, nsnps, h2_total) %>% 
  mutate(cv_split = row_number()) %>% 
  ungroup()

# make prop dom percent
devdf$prop_dom_percent = devdf$prop_dom_percent = round((devdf$prop_dom / devdf$h2_total), 2)

# FORMAT FOR EACH POLYGENICITY LEVEL ----

# 100 SNP phenotypes
devdf = devdf.within.pheno.results(devdf, h2_vals = c(0.5, 0.2, 0.1))

# ---- Get y-axis ranges

# Step 1: Group by nsnps, h2_total, and prop_dom
r2_ranges = devdf %>%
  group_by(h2_total) %>%
  summarise(
    r2_min = min(mean_r2, na.rm = TRUE),
    r2_max = max(mean_r2, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(r2_range = r2_max - r2_min)

# Step 2: Get the largest range and add a small buffer
y_span = max(r2_ranges$r2_range, na.rm = TRUE) + 0.01
y_span = 0.09

# Step 3: make max and min ylims per phenotype

# h = 10%
max_h01 = 0.11
min_h01 = max_h01 - 0.09

# h = 20%
max_h02 = 0.21
min_h02 = max_h02 - 0.09

# h = 50%: #2 
max_h05 = 0.51
min_h05 = max_h05 - 0.1
min_h05_3b = 0.43

# make y ticks
brks01 = c(0.1, 0.08, 0.06, 0.04, 0.02)
brks02 = c(0.2, 0.18, 0.16, 0.14, 0.12)
brks05 = c(0.5, 0.48, 0.46, 0.44, 0.42)

# ---- make figures

# set tick and title size
titleSize = 15
tickSize = 15

# figure 1: 1 column (h2 = 50%) x 3 rows
# rows = number of causal SNPs
plot.h01 = plot_devdf_summary_r2(devdf, h2 = 0.1, title = NULL, 
                           title_size = titleSize, tick_size = tickSize, 
                           ylims = c(min_h01, max_h01), brks = brks01) + theme(axis.title.x = element_blank())
plot.h02 = plot_devdf_summary_r2(devdf, h2 = 0.2, title = NULL, 
                           title_size = titleSize, tick_size = tickSize,
                           ylims = c(min_h02, max_h02), brks = brks02) + theme(axis.title.x = element_blank())

plot.h05 = plot_devdf_summary_r2(devdf, h2 = 0.5, title = NULL,  
                           title_size = titleSize, tick_size = tickSize,
                           ylims = c(min_h05_3b, max_h05), brks = brks05)

# MAKE FIGURE 
fig3b_plot =
  (plot.h01) /
  (plot.h02) /
  (plot.h05) +
  plot_layout(guides = "collect") &
  theme(legend.position = "bottom")

# View in R
print(fig3b_plot)

# save and view
ggplot2::ggsave(
  filename = "~/Desktop/fig3b_plot.png",
  plot = fig3b_plot,
  dpi = 300
)






