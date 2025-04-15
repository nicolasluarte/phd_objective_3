# lib load ----
pacman::p_load(
    tidyverse,
    ggplot2,
    patchwork
)
setwd(this.path::here())

# Plot stuff ----
theme_uncertainty <- ggpubr::theme_pubr() +
    update_geom_defaults("point", list(size = 5, alpha = 0.5, shape = 21)) +
    theme(
        text = element_text(size = 24),
        axis.text=element_text(size=14),
        plot.margin = unit(c(0.5,0.5,0.5,0.5), "cm"),
        legend.position = "none"
    )

boxplot_sig_bracket <- function(group1, group2){
    ggsignif::geom_signif(
        comparisons = list(c(group1, group2)),
        map_signif_level = TRUE,
        textsize = 0,
        tip_length = 0
    )
}

# Behavior ----
bh_d <- read_csv("../datasets/")

## mdls ----

### licks and events ----
bh_mdl <- lme4::glmer.nb(
    data = bh_d,
    Licks ~ TRT_CTRL * tipo_recompensa + n_sesion + (1|ID),
    control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(bh_mdl)


    
    
    
    
    
    
    
    
    
    