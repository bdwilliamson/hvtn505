# create the forest plot of CV-AUC for each of the top learner and SL from each assay combination
auc_plot_assays <- function(avg_aucs, main_font_size_forest, main_font_size_lab, sl_only = TRUE, immunoassay = TRUE, colors = NULL) {
  ## get SL and the top individual learner/screen combo from each assay
  top_learners_init <- avg_aucs %>% 
    group_by(assay) %>% 
    arrange(desc(AUC), .by_group = TRUE) %>% # arrange in descending order by assay group
    filter((Learner == "SL" & Screen == "All") | row_number() == 1) %>%  # select only the SL row and the top performing learner 
    mutate(learner_nm = make_nice_learner_name(Learner), screen_nm = make_nice_screen_name(Screen)) %>% 
    ungroup() %>% 
    arrange(desc(AUC))
  
  ## if sl_only, only use the SL algorithm for each assay
  if (sl_only) {
    top_learners <- top_learners_init %>% 
      filter(Learner == "SL")
  } else {
    top_learners <- top_learners_init
  }
  
  ## forest plot with the groups/learners from the analysis plan
  top_learner_plot <- top_learners %>% 
    ggplot(aes(x = AUC, y = factor(paste0(Screen, "_", Learner, "_", assay), levels = paste0(Screen, "_", Learner, "_", assay)[order(AUC)],
                                   labels = paste0(varset_label, " ", learner_nm, " ", screen_nm)[order(AUC)]))) + 
    geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
    geom_point() +
    xlab("CV-AUC") +
    ylab("") + 
    scale_x_continuous(breaks = seq(0.4, 1, 0.1),
                     labels = as.character(seq(0.4, 1, 0.1)),
                     limits = c(0.4, 1)) +
    theme(legend.position = "", 
          axis.text.y = element_blank(),
          text = element_text(size = main_font_size_forest),
          axis.title = element_text(size = main_font_size_forest), 
          axis.text.x = element_text(size = main_font_size_forest),
          axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0), size = main_font_size_forest),
          plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
  if (immunoassay) {
    ## add on a legend to top_learner_plot, color based on immunoassay set
    top_learner_plot <- top_learner_plot +
      geom_errorbarh(aes(xmin = top_learners$ci_ll, xmax = top_learners$ci_ul, color = top_learners$immunoassay_set)) +
      scale_color_manual(values = colors) +
      labs(color = "Immunoassay set") +
      theme(legend.position = c(0.63, 0.2), 
            axis.text.y = element_blank(),
            text = element_text(size = main_font_size_forest),
            axis.title = element_text(size = main_font_size_forest), 
            axis.text.x = element_text(size = main_font_size_forest),
            axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0), size = main_font_size_forest),
            plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
    ## separate plot with nice names, printed values of AUCs, based on immunoassays only
    top_learners_labels <- top_learners %>% 
      ungroup() %>% 
      mutate(lab_auc = paste0(format(round(AUC, 3), nsmall = 3), " [", 
                              format(round(ci_ll, 3), nsmall = 3), ", ", 
                              format(round(ci_ul, 3), nsmall = 3), "]")) %>%
      select(varset_label, lab_auc)
    ## melt to make a single "value" column
    top_learners_labels$var <- 1
    top_learners_labs <- melt(top_learners_labels, id.var = "var")
    ## tack on x, y coordinates
    top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1, function(x) which(grepl(x, c("varset_label", "lab_auc"))) - 1 + c(0.2, 0.5)[which(grepl(x, c("varset_label", "lab_auc")))])
    top_learners_labs$y_coord <- rep(rev(as.numeric(rownames(top_learners))), dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
  } else {
    ## separate plot with nice names, printed values of the AUCs
    top_learners_labels <- top_learners %>% 
      ungroup() %>% 
      mutate(lab_auc = paste0(format(round(AUC, 3), nsmall = 3), " [", 
                              format(round(ci_ll, 3), nsmall = 3), ", ", 
                              format(round(ci_ul, 3), nsmall = 3), "]")) %>%
      select(varset_label, learner_nm, screen_nm, lab_auc)
    ## melt to make a single "value" column
    top_learners_labels$var <- 1
    top_learners_labs <- melt(top_learners_labels, id.var = "var")
    ## tack on x, y coordinates
    top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1, function(x) which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_auc"))) - 1 + c(0, 0.3, -0.3, 0.2)[which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_auc")))])
    top_learners_labs$y_coord <- rep(rev(as.numeric(rownames(top_learners))), dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
  }
  ## make the plot
  top_learner_nms_plot <- top_learners_labs %>% 
    ggplot(aes(x = x_coord, y = y_coord, label = value)) +
    geom_text(size = main_font_size_lab, hjust = 0, vjust = 0.5) +
    xlim(c(min(top_learners_labs$x_coord), max(top_learners_labs$x_coord) + 0.75)) +
    ylim(c(1, max(top_learners_labs$y_coord))) +
    theme(legend.position="", 
          axis.line=element_blank(),
          axis.text=element_blank(),
          text = element_text(size = main_font_size_lab),
          axis.title = element_blank(), 
          axis.ticks = element_blank(),
          plot.margin=unit(c(0,0,0,0),"cm"),
          panel.background=element_blank(),
          panel.border=element_blank(),
          panel.grid.major=element_blank(),
          panel.grid.minor=element_blank(),
          plot.background=element_blank())
  
  
  return(list(top_learner_plot = top_learner_plot, top_learner_nms_plot = top_learner_nms_plot))
}