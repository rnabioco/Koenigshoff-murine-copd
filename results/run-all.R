library(rmarkdown)
library(here)
rmds <- c("00_other/x_atlas.Rmd",
          "01_tdtomato_egfp_expt1/1_preprocess.Rmd",
          "01_tdtomato_egfp_expt1/2_indepth.Rmd",
          "03_tdtomato_egfp_expt2/00_processing.Rmd",
          "04_combined_analyses/combined_2.Rmd",
          "04_combined_analyses/combined_2.Rmd",
          "N_cellbrowser.Rmd")

for (i in seq_along(rmds)){
  render(file.path(here("results", rmds[i])))
}
