# create the forest plot of the specified measure for each of the
# top learner and SL from each assay combination
# @param avgs the average predictiveness measures
# @param type the predictiveness measure ("auc" or "r_squared")
# @param main_font_size_forest the main font size for text in the forest plot
# @param main_font_size_lab the font size for axis labels
# @param sl_only should we only plot the SL estimates, or all learners?
# @param immunoassay should we color by immunoassay?
# @param colors the colors to use
# @param point_size what size should the points be?
plot_assays <- function(avgs = NULL, type = "auc", main_font_size_forest = 5,
                        main_font_size_lab = 5, sl_only = TRUE,
                        immunoassay = TRUE,
                        colors = NULL, point_size = 3, x_lim = c(0.4, 1), 
                        lgnd_pos = c(0.56, 0.2)) {
  # if type is AUC, make correct labels
  if (type == "auc") {
    x_lab <- "CV-AUC"
    avgs <- avgs %>%
      mutate(measure = AUC)
  } else {
    x_lab <- expression(paste("CV-", R^2))
    # x_lim <- c(-1, 1)
    avgs <- avgs %>%
      mutate(measure = R2)
    # lgnd_pos <- c(0.65, 0.15)
  }
  # get SL and the top individual learner/screen combo from each assay
  top_learners_init <- avgs %>%
    group_by(assay) %>%
    arrange(desc(measure), .by_group = TRUE) %>%
    filter((Learner == "SL" & Screen == "All") | row_number() == 1) %>%
    mutate(learner_nm = make_nice_learner_name(Learner), 
           screen_nm = make_nice_screen_name(Screen)) %>%
    ungroup() %>%
    arrange(desc(measure)) %>%
    mutate(ci_ll_truncated = pmax(ci_ll, x_lim[1]),
           ci_ul_truncated = pmin(ci_ul, x_lim[2]))

  # if sl_only, only use the SL algorithm for each assay
  if (sl_only) {
    top_learners <- top_learners_init %>%
      filter(Learner == "SL")
  } else {
    top_learners <- top_learners_init
  }

  # forest plot with the groups/learners from the analysis plan
  top_learner_plot <- top_learners %>%
    ggplot(aes(y = measure,
        x = factor(paste0(Screen, "_", Learner, "_", assay),
        levels = paste0(Screen, "_", Learner, "_", assay)[order(measure)],
        labels = paste0(varset_label, " ", learner_nm, " ", 
                        screen_nm)[order(measure)]),
        xend = factor(paste0(Screen, "_", Learner, "_", assay),
        levels = paste0(Screen, "_", Learner, "_", assay)[order(measure)],
        labels = paste0(varset_label, " ", learner_nm, " ", 
                        screen_nm)[order(measure)]))) +
    geom_errorbar(aes(ymin = ci_ll, ymax = ci_ul), size = point_size) +
    geom_pointrange(aes(ymin = ci_ll_truncated, ymax = ci_ul_truncated),
    size = point_size) +
    ylab(x_lab) +
    xlab("") +
    scale_y_continuous(breaks = round(seq(x_lim[1], x_lim[2], 0.1), 1),
                       labels = as.character(round(seq(x_lim[1], 
                                                       x_lim[2], 0.1), 1)),
                       limits = x_lim) +
    coord_flip() +
    theme(legend.position = "",
          axis.text.y = element_blank(),
          text = element_text(size = main_font_size_forest),
          axis.title = element_text(size = main_font_size_forest),
          axis.text.x = element_text(size = main_font_size_forest),
          axis.ticks = element_line(size = 1),
          axis.ticks.length = unit(0.5, "cm"),
          axis.title.x = element_text(margin = ggplot2::margin(t = 5,
              r = 0, b = 0, l = 0), size = main_font_size_forest),
          plot.margin=unit(c(1.25,0.5,0.5,0),"cm")) # top, right, bottom, left
  if (immunoassay) {
    # add on a legend to top_learner_plot, color based on immunoassay set
    top_learner_plot <- top_learner_plot +
      geom_errorbar(aes(ymin = ci_ll,
                        ymax = ci_ul,
                        color = immunoassay_set),
                    size = point_size) +
      geom_pointrange(aes(ymin = ci_ll_truncated,
                         ymax = ci_ul_truncated,
                         color = immunoassay_set),
                     size = point_size) +
      scale_color_manual(values = colors) +
      labs(color = "Assay set") +
      theme(legend.position = lgnd_pos,
            axis.text.y = element_blank(),
            text = element_text(size = main_font_size_forest),
            axis.title = element_text(size = main_font_size_forest),
            axis.text.x = element_text(size = main_font_size_forest),
            axis.title.x = element_text(margin = ggplot2::margin(t = 5,
                 r = 0, b = 0, l = 0), size = main_font_size_forest),
            legend.text = element_text(size = main_font_size_forest),
            plot.margin=unit(c(1.25,0.5,0.5,0),"cm")) 
    
    # separate plot with nice names,
    # printed values of measures, based on immunoassays only
    top_learners_labels <- top_learners %>%
      ungroup() %>%
      mutate(lab_measure = paste0(format(round(measure, 3), nsmall = 3), " [",
                              format(round(ci_ll, 3), nsmall = 3), ", ",
                              format(round(ci_ul, 3), nsmall = 3), "]")) %>%
      select(varset_label, lab_measure)
    # melt to make a single "value" column
    top_learners_labels$var <- 1
    top_learners_labs <- reshape2::melt(top_learners_labels, id.var = "var")
    # tack on x, y coordinates
    top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1,
    function(x) which(grepl(x, c("varset_label", "lab_measure"))) -
    1 + c(0.2, 0.5)[which(grepl(x, c("varset_label", "lab_measure")))])
    top_learners_labs$y_coord <- rep(rev(as.numeric(rownames(top_learners))),
     dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
  } else {
    # separate plot with nice names, printed values of the measures
    top_learners_labels <- top_learners %>%
      ungroup() %>%
      mutate(lab_measure = paste0(format(round(measure, 3), nsmall = 3), " [",
                              format(round(ci_ll, 3), nsmall = 3), ", ",
                              format(round(ci_ul, 3), nsmall = 3), "]")) %>%
      select(varset_label, learner_nm, screen_nm, lab_measure)
    # melt to make a single "value" column
    top_learners_labels$var <- 1
    top_learners_labs <- reshape2::melt(top_learners_labels, id.var = "var")
    # tack on x, y coordinates
    top_learners_labs$x_coord <- apply(matrix(top_learners_labs$variable), 1,
    function(x) which(grepl(x,
        c("varset_label", "learner_nm", "screen_nm", "lab_measure"))) -
        1 + c(0, 0.3, -0.3, 0.2)[which(grepl(x,
            c("varset_label", "learner_nm", "screen_nm", "lab_measure")))])
    top_learners_labs$y_coord <- rep(rev(as.numeric(rownames(top_learners))),
     dim(top_learners_labs)[1]/length(as.numeric(rownames(top_learners))))
  }
  # make the plot
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


  return(list(top_learner_plot = top_learner_plot, 
              top_learner_nms_plot = top_learner_nms_plot))
}
