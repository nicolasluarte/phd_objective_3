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
bh_d <- read_csv("../datasets/lickometer_data_proc.csv") %>% 
    arrange(ID, group, rel_date)
bh_d

## mdls ----

### licks and events ----
bh_licks_mdl <- glmmTMB::glmmTMB(
    data = bh_d,
    licks ~ tipo_recompensa * group + rel_date + (1|ID),
    family = glmmTMB::nbinom2
)
summary(bh_licks_mdl)

bh_licks_emm <- emmeans::emmeans(
    bh_licks_mdl,
    pairwise ~ tipo_recompensa * group | rel_date,
    type = "response"
)
bh_licks_emm

bh_events_mdl <- glmmTMB::glmmTMB(
    data = bh_d,
    events ~ tipo_recompensa * group + rel_date + (1|ID),
    family = glmmTMB::nbinom2
)
summary(bh_events_mdl)

bh_events_emm <- emmeans::emmeans(
    bh_events_mdl,
    pairwise ~ tipo_recompensa * group | rel_date,
    type = "response"
)
bh_events_emm


# Neurons ----

groups <- bh_d %>% 
    select(ID, group) %>% 
    distinct(ID, group)
nd <- read_csv("../datasets/cell_count_matrix.csv")

oxa_positive_neurons <- nd %>% 
    filter(oxa_positive==1) %>% 
    group_by(animal_id) %>% 
    summarise(
        oxa_plus = n()
    ) %>% 
    left_join(., groups, by = c("animal_id"="ID"))
oxa_positive_neurons

oxa_rb_cfos <- nd %>% 
    left_join(., groups, by = c("animal_id"="ID")) %>% 
    group_by(animal_id, oxa_positive, rb_positive, group) %>% 
    summarise(
        cnt = n(),
        cfos_positive = sum(cfos_positive)
    ) %>% 
    filter(!(oxa_positive==0 & rb_positive==0)) %>% 
    mutate(
        percent_of_total = (cfos_positive/cnt),
        class = case_when(
            oxa_positive == 1 & rb_positive == 1 ~ "oxa+_rb+",
            oxa_positive == 1 & rb_positive == 0 ~ "oxa+_rb-",
            oxa_positive == 0 & rb_positive == 1 ~ "oxa-_rb+"
        ),
        successes = cfos_positive,
        failures = cnt - cfos_positive 
    )
oxa_rb_cfos

oxa_rb_cfos_mdl <- lme4::glmer(
    data = oxa_rb_cfos,
    cbind(successes, failures) ~ class * group + (1|animal_id),
    family = binomial(link = "logit"),
    control = lme4::glmerControl(
        optimizer = "bobyqa",
        optCtrl = list(maxfun = 2e5)
    )
)
summary(oxa_rb_cfos_mdl)

oxa_rb_cfos_emm <- emmeans::emmeans(
    oxa_rb_cfos_mdl,
    pairwise ~ group | class,
    type = "response"
)
oxa_rb_cfos_emm



## mdls ----
### oxa positive neurons ----

oxa_pos_mdl <- lm(
    data = oxa_positive_neurons,
    oxa_plus ~ group
)
summary(oxa_pos_mdl)

# p::oxa_positive_neurons ----

np1 <- oxa_positive_neurons %>% 
    ggplot(aes(
        group, oxa_plus
    )) +
    geom_boxplot(outlier.shape = NA, aes(color=group), width = 0.5) + 
    geom_point(aes(fill=group)) +
    theme_uncertainty +
    scale_x_discrete(labels = c("Ctrl.", "Trt.")) +
    scale_y_continuous(breaks = seq(0,500,100), 
                       limits = c(0,500), 
                       expand = c(0,0)) +
    ylab("OXA+ neurons") +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
np1

# p::oxa_rb_cfos ----

np2 <- oxa_rb_cfos %>% 
    mutate(class_plot=factor(as.factor(interaction(class, group)),
                             levels = c("oxa+_rb-.control",
                                        "oxa+_rb-.treatment",
                                        "oxa-_rb+.control",
                                        "oxa-_rb+.treatment",
                                        "oxa+_rb+.control",
                                        "oxa+_rb+.treatment"
                                        ))) %>% 
    ggplot(aes(
        class_plot, percent_of_total
    )) +
    geom_boxplot(outlier.shape = NA, aes(color = group), width = 0.5) +
    geom_point(aes(fill = group), position = position_dodge(0.75)) +
    ggsignif::geom_signif(
        comparisons = list(c(5, 6)),
        textsize = 0,
        y_position = 0.47
    ) +
    theme_uncertainty +
    scale_y_continuous(breaks = seq(0,0.6,0.1), 
                       limits = c(0,0.6), 
                       expand = c(0,0)) +
    ylab("cFos+") +
    xlab("") +
    scale_x_discrete(labels = c(
        latex2exp::TeX(r"($OXA^{+}_{Rb^{-}}$)"),
        "",
        latex2exp::TeX(r"($OXA^{-}_{Rb^{+}}$)"),
        "",
        latex2exp::TeX(r"($OXA^{+}_{Rb^{+}}$)"),
        ""
    )) + 
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) +
    theme(
        axis.text.x = element_text(hjust = -0.1),
        axis.ticks.x = element_blank()
    )
np2
    
# p::licks ----

bhp1 <- bh_d %>% 
    ungroup() %>% 
    group_by(ID, tipo_recompensa, group) %>% 
    summarise(
        licks = mean(licks)
    ) %>% 
    mutate(
        x_lab = factor(as.factor(interaction(group, tipo_recompensa)), levels =
                           c("treatment.sac", "treatment.wat", "control.wat"))
    ) %>% 
    ggplot(aes(
        x_lab, licks
    )) +
    geom_boxplot(outlier.shape = NA, aes(color = group)) +
    geom_point(aes(fill=group)) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 3), c(1, 2)),
        textsize = 0,
        y_position = c(1250, 1200)
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c(latex2exp::TeX(r"($Suc_{trt}$)"),
                                latex2exp::TeX(r"($Wat_{trt}$)"),
                                latex2exp::TeX(r"($Wat_{ctrl}$)"))) +
    scale_y_continuous(breaks = seq(0,1400,200), 
                       limits = c(0,1400), 
                       expand = c(0,0)) +
    ylab("# Licks") +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
bhp1

# p::rewards ----

bhp2 <- bh_d %>% 
    ungroup() %>% 
    group_by(ID, tipo_recompensa, group) %>% 
    summarise(
        licks = mean(events)
    ) %>% 
    mutate(
        x_lab = factor(as.factor(interaction(group, tipo_recompensa)), levels =
                           c("treatment.sac", "treatment.wat", "control.wat"))
    ) %>% 
    ggplot(aes(
        x_lab, licks
    )) +
    geom_boxplot(outlier.shape = NA, aes(color = group)) +
    geom_point(aes(fill=group)) +
    ggsignif::geom_signif(
        comparisons = list(c(1, 3), c(1, 2)),
        textsize = 0,
        y_position = c(50, 45)
    ) +
    theme_uncertainty +
    scale_x_discrete(labels = c(latex2exp::TeX(r"($Suc_{trt}$)"),
                                latex2exp::TeX(r"($Wat_{trt}$)"),
                                latex2exp::TeX(r"($Wat_{ctrl}$)"))) +
    scale_y_continuous(breaks = seq(0,60,10), 
                       limits = c(0,60), 
                       expand = c(0,0)) +
    ylab("# Rewards") +
    xlab("") +
    scale_fill_manual(values = c("black", "orange")) +
    scale_color_manual(values = c("black", "orange")) 
bhp2

# figure ----

bhp1 | bhp2 | np1 | np2
    
    