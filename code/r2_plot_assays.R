# create the forest plot of CV-R2 for each of the top learner and SL from each assay combination
r2_plot_assays <- function(avg_r2s, main_font_size) {
  ## get SL and the top individual learner/screen combo from each assay
  top_learners <- avg_r2s %>% 
    group_by(assay) %>% 
    arrange(desc(R2), .by_group = TRUE) %>% # arrange in descending order by assay group
    filter((Learner == "SL" & Screen == "All") | row_number() == 1) %>%  # select only the SL row and the top performing learner 
    mutate(learner_nm = make_nice_learner_name(Learner), screen_nm = make_nice_screen_name(Screen)) %>% 
    ungroup() %>% 
    arrange(desc(R2))
  
  ## another forest plot with the groups/learners from the analysis plan
  top_learner_plot <- top_learners %>% 
    ggplot(aes(x = R2, y = factor(paste0(Screen, "_", Learner, "_", assay), levels = paste0(Screen, "_", Learner, "_", assay)[order(R2)],
                                   labels = paste0(varset_label, " ", learner_nm, " ", screen_nm)[order(R2)]))) + 
    geom_errorbarh(aes(xmin = ci_ll, xmax = ci_ul)) +
    geom_point() +
    xlab(expression(paste("CV-", R^2))) +
    ylab("") + 
    xlim(c(-0.6, 0.5)) +
    theme(legend.position="", 
          axis.text.y=element_blank(),
          text = element_text(size=35),
          axis.title = element_text(size=25), 
          axis.title.x = element_text(margin = ggplot2::margin(t = 20, r = 0, b = 0, l = 0)),
          plot.margin=unit(c(1,0.5,0,0),"cm")) # top, right, bottom, left
  ## separate plot with nice names, printed values of the r2s
  top_learners_labels <- top_learners %>% 
    ungroup() %>% 
    mutate(lab_r2 = paste0(format(round(R2, 3), nsmall = 3), " [", 
                            format(round(ci_ll, 3), nsmall = 3), ", ", 
                            format(round(ci_ul, 3), nsmall = 3), "]")) %>%
    select(varset_label, learner_nm, screen_nm, lab_r2)
  ## melt to make a single "value" column
  top_learners_labels$var <- 1
  top_learners_labs <- melt(top_learners_labels, id.var = "var")
  ## tack on x, y coordinates
  top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1, function(x) which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_r2"))) - 1 + c(0, 0.3, -0.1, 0.2)[which(grepl(x, c("varset_label", "learner_nm", "screen_nm", "lab_r2")))])
  top_learners_labs$y_coord <- rep(rev(as.numeric(rownames(top_learners))), dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
  top_learner_nms_plot <- top_learners_labs %>% 
    ggplot(aes(x = x_coord, y = y_coord, label = value)) +
    geom_text(size = main_font_size, hjust = 0, vjust = 0.5) +
    xlim(c(min(top_learners_labs$x_coord), max(top_learners_labs$x_coord) + 0.75)) +
    ylim(c(1, max(top_learners_labs$y_coord))) +
    theme(legend.position="", 
          axis.line=element_blank(),
          axis.text=element_blank(),
          text = element_text(size = main_font_size),
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