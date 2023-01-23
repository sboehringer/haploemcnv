#' @title Nested-loop plot
#' 
#' @description Creating of a nested-loop plot to visualize the estimated haplotype frequencies and effect sizes for a large number of simulation scenarios. The nested-loop plot is based on \cite{}.
#'
#' @details ... For the cut-off value of the outliers, \code{outlier_prob} is used as quantile for the distribution of all values to plot in \code{EM_betas} + three times the interquantile range. All values above this cut-off value are replaced by empty triangles and plotted with the cut-off value as their \emph{true} value. Another option would be to manually choose the cut-off value via \code{maxY}. But this is less elegant, because then the value is not based on the data and 
#'
#' @param Scenario_overview A matrix. Overview of the characteristics of all simulation scenarios
#' @param EM_betas A list. Each list element are the values to plot obtained via a different EM-algorithm analysis method (\emph{e.g.,} with vs without outcome analysis). 
#' @param naive_betas An optional list. Each list element are the values to plot obtained via a different naive analysis method (\emph{e.g.,} dictorial vs democratic analysis).
#' @param list_variable An optional character. Variable name that explains the different list elements of \code{EM_betas} and \code{naive_betas}.
#' @param variables A vector. The variables where the simulation scenarios differ from each other. Their values must be stated in \cite{Scenario_overview}.
#' @param mean_plot A logical scalar. Whether or not the mean value of the values to plot for each scenario must be plotted as a horizontal line.
#' @param title An optional vector. Title of the Nested-loop plot.
#' @param yname A character. The name of the y-axis, how the values to plot are calculated. 
#' @param outlier_prob An optional numeric value \in \{0, 1\}. Indicating which part of the distribution is discarded and placed as empty triangles on top of the graph.
#' @param maxY An optional numeric value. A crude cut-off value, replacing the one calculated with \code{outlier_prob}. Less sophisticated and not data-driven. 
#' @param extra_lines An optional list. Additional values to plot as colored horizontal lines for each scenario. Predominantly used for the individual haplotype frequencies.
#' @param second_Yaxis An optional vector. If supplied, the values of the y-axis will be transformed with the specified value and also plotted on the right. Useful if the values in \code{extra_lines} have been scaled. The name of \code{second_Yaxis} will be used as the name for the second y-axis.
#' @param altered_replications An optional vector. Indicates per scenario how many replications have been replaced due to the Euclidean distance estimates.  
#' 
#' @return The nested-loop plot and the order in which the scenarios are plotted
#' 
Nested_loop_plot = function(Scenario_overview, EM_betas, naive_betas = NULL, list_variable = "outcome", HTR_name = NULL, variables = c("naive", "LD", "Amb"), mean_plot = FALSE, 
                            mean_name = "mean", title = NULL, yname = "RMSE", outlier_prob = 0.99, maxY = NULL, extra_lines = NULL, second_Yaxis = NULL, altered_replications = NULL){
  

  if("naive" %in% variables & !is.null(naive_betas)){
    naive = TRUE
  } else {
    naive = FALSE
    if("naive" %in% variables){
      variables = variables[variables != "naive"]
    }
  }
  
  
  

##### Setting-up Data #####
  if(is.list(EM_betas)){
    EM_betas_list = EM_betas
    EM_betas = EM_betas_list[[1]]
    
    nr_lists = length(EM_betas_list); lists = TRUE
    
    if(naive){
      naive_betas_list = naive_betas
      naive_betas = naive_betas_list[[1]]
      
      if(length(naive_betas_list) != nr_lists){
        stop("Please supply the same number of lists for both the EM as naive analysis.")
      }
    }
  } else {
    lists = FALSE
  }
  
  
  Scenario_overview_sub = Scenario_overview[!duplicated(Scenario_overview[, colnames(Scenario_overview) %in% c("Exp_name", variables)]), ]
  data_lines = Scenario_overview_sub[Scenario_overview_sub[, "Exp_name"] %in% colnames(EM_betas), ]
  
  EM_betas = EM_betas[, colnames(EM_betas) %in% Scenario_overview_sub[, "Exp_name"]]

  
  
  if(lists){
    if(naive){
      naive_betas = naive_betas[, colnames(naive_betas) %in% Scenario_overview_sub[, "Exp_name"]]

      scenario_rep = nr_lists * 2
    } else {
      scenario_rep = nr_lists
      
    }
  } else {
    if(naive){
      naive_betas = naive_betas[, colnames(naive_betas) %in% Scenario_overview_sub[, "Exp_name"]]

      scenario_rep = 2
    } else {
      scenario_rep = 1
      
    }
  }
  nr_scenarios = ncol(EM_betas) * scenario_rep
  
  
  
  nr_variables = length(variables)
  for(i in nr_variables:1){
    if(variables[i] != "naive"){
      data_lines = data_lines[order(data_lines[, variables[i]]), ]
    }
  }
  data_lines$Model = (1:nrow(data_lines) * scenario_rep) - (scenario_rep - 1)
  
  
  
  
##### Variables #####
  variables = variables[variables != "Exp_name"]
  
  index1 = variables %in% colnames(data_lines)[lengths(apply(data_lines, 2, unique)) != 1]
  index2 = variables == "naive"
  
  variables = variables[index1 | index2]; nr_variables = length(variables)
  
  if(lists){
    if(variables[1] == "naive"){
      variables = c(variables[1], list_variable, variables[2:nr_variables])
    } else {
      variables = c(list_variable, variables)
    }
    
    nr_variables = length(variables)
    
    list_variable_labels = names(EM_betas_list)
  }
  
  
  
  variable_labels = variable_values = NULL
  for(i in 1:nr_variables){
    variable = variables[i]
    
    if(variable == "naive"){
      variable_labels[i] = c("Model\n(EM, MV)")
    } else if(variable == list_variable){
      variable_labels[i] = paste(list_variable, "\n(", paste(names(EM_betas_list), collapse = ", "),")", sep = "")
    } else {  
      variable_labels[i] = paste(variable, "\n(", paste(unique(sort(data_lines[, variable])), collapse = ", "), ")", sep = "")
    }
  }
  
  variable_labels1 = unlist(lapply(strsplit(variable_labels, "\\\n"), function(x){x[1]}))
  variable_labels2 = unlist(lapply(strsplit(variable_labels, "\\\n"), function(x){x[2]}))
  
  
  
##### Dataset #####
  if(!lists){
    dat_EM = as.data.frame(as.table(EM_betas))
    colnames(dat_EM) = c("Haplo", "Scenario", "RMSE")
    dat_EM$naive = 0
    
    if(naive){
      dat_naive = as.data.frame(as.table(naive_betas))
      colnames(dat_naive) = c("Haplo", "Scenario", "RMSE")
      dat_naive$naive = 1  
      
      plot_dat = rbind(dat_EM, dat_naive)
    } else {
      plot_dat = dat_EM
    }
    
  } else {
    dat_EM = NULL
    for(i in 1:nr_lists){
      dat_EM_i = as.data.frame(as.table(EM_betas_list[[i]]))
      colnames(dat_EM_i) = c("Haplo", "Scenario", "RMSE")
      dat_EM_i$naive = 0
      
      ##
      if(!is.logical(list_variable_labels[i])){
        dat_EM_i = cbind(dat_EM_i, "NEW" = as.character(list_variable_labels[i]))
      } else {
        dat_EM_i = cbind(dat_EM_i, "NEW" = as.integer(as.logical(list_variable_labels[i])))
      }
      ##
      colnames(dat_EM_i)[colnames(dat_EM_i) == "NEW"] = list_variable
      
      dat_EM = rbind(dat_EM, dat_EM_i)
    }
    
    if(naive){
      dat_naive = NULL
      for(i in 1:nr_lists){
        dat_naive_i = as.data.frame(as.table(naive_betas_list[[i]]))
        colnames(dat_naive_i) = c("Haplo", "Scenario", "RMSE")
        dat_naive_i$naive = 1
        
        ##
        if(!is.logical(list_variable_labels[i])){
          dat_naive_i = cbind(dat_naive_i, "NEW" = as.character(list_variable_labels[i]))
        } else {
          dat_naive_i = cbind(dat_naive_i, "NEW" = as.integer(as.logical(list_variable_labels[i])))
        }
        ##
        colnames(dat_naive_i)[colnames(dat_naive_i) == "NEW"] = list_variable
        
        dat_naive = rbind(dat_naive, dat_naive_i)
      }
      
      plot_dat = rbind(dat_EM, dat_naive)
    } else {
      
      plot_dat = dat_EM
    }
  }
  plot_dat = plot_dat[order(plot_dat$Scenario), ]
  plot_dat_var = c("Model", variables[!variables %in% colnames(plot_dat)]); len_plot_dat_var = length(plot_dat_var)
  
  
  ##
  #...$Model contains the order in which the data is shown, need to make that good:
  #     Now have: EM_noY; EM_Y; naive_noY; naive_Y
  #
  #     This code can be made much faster
  ##
  
  
  plot_dat = cbind(plot_dat, matrix(NA, nrow = nrow(plot_dat), ncol = len_plot_dat_var, dimnames = list(NULL, plot_dat_var)))
  if(lists){
    for(i in 1:nrow(plot_dat)){
      index = which(data_lines[, "Exp_name"] == plot_dat$Scenario[i])
      
      for(j in 1:len_plot_dat_var){
        plot_dat[i, plot_dat_var[j]] = data_lines[index, plot_dat_var[j]]
      }
      
      if(is.logical(list_variable_labels[1])){
        plot_dat[i, "Model"] = data_lines[index, "Model"] + 0.5 + plot_dat[i, list_variable] * 1 + plot_dat[i, "naive"] * 2
      } else {
        plot_dat[i, "Model"] = data_lines[index, "Model"] + 0.5 + (which(plot_dat[i, list_variable] == list_variable_labels) - 1) * 1 + plot_dat[i, "naive"] * 2
      }
    }
    
  } else {
    for(i in 1:nrow(plot_dat)){
      index = which(data_lines[, "Exp_name"] == plot_dat$Scenario[i])
      
      for(j in 1:len_plot_dat_var){
        plot_dat[i, plot_dat_var[j]] = data_lines[index, plot_dat_var[j]]
      }
      
      plot_dat[i, "Model"] = data_lines[index, "Model"] + 0.5 + plot_dat[i, "naive"]
    }
  }

  plot_order = sort(tapply(plot_dat$Model, plot_dat$Scenario, mean))
  scenarios = names(plot_order)
  
  
#  if(TRUE %in% (plot_dat$RMSE > 15)){
#    plot_dat$RMSE = sapply(plot_dat$RMSE, function(x){ifelse(x > 3, (log(x, 10) + 2), x)})  # The plus 2 is to make sure that all values above 3 (taken the log of) are higher than 3
#    log_scale = TRUE
#  } else {
#    log_scale = FALSE
#  }

  

    
##### Outliers #####
  # if(is.null(maxY)){
  outlier_value = quantile(plot_dat$RMSE, probs = outlier_prob, na.rm = TRUE) + 3 * IQR(plot_dat$RMSE, na.rm = TRUE)
  # } else {
  #   outlier_value = maxY
  # }
  
  y_lim_max = max(plot_dat$RMSE[plot_dat$RMSE < outlier_value], na.rm = TRUE) * 1.05  # outlier_value * 1.05
  y_lim_min = -(y_lim_max * 0.25)
  cat("With this value", sum(plot_dat$RMSE > y_lim_max, na.rm = TRUE), "values have been collapsed\n")
    
    
  plot_dat2 = plot_dat[plot_dat$RMSE > outlier_value & !is.na(plot_dat$RMSE), ]
  plot_dat1 = plot_dat[!(rownames(plot_dat) %in% rownames(plot_dat2)), ]

  
#  y_lim_max = ceiling(max(plot_dat$RMSE, na.rm = TRUE) / 0.5) * 0.5
#  y_lim_min = min(unlist(lapply(line_coord, function(x){min(x$y)})))
  
#  continue = TRUE; nr_y_lines = nr_y_lines - 1
#  while(continue & nr_y_lines < (y_lim_max + 1)){
#    nr_y_lines = nr_y_lines + 1
#    
#    y_range = seq(from = 0, to = y_lim_max, length.out = nr_y_lines)
#    
#    if(!(FALSE %in% (y_range %% 0.5 == 0))){
#      continue = FALSE
#    }
#  }
#  
#  if(continue){  # If above was correct (continue has become false) this is skipped
#    y_range = seq(from = 0, to = ceiling(y_lim_max), by = 1)
#  }
  
  

##### Extra lines #####
  if(mean_plot){
    if(lists){
      EM_cmeans = NULL
      for(i in 1:nr_lists){
        EM_cmeans[[i]] = colMeans(EM_betas_list[[i]], na.rm = TRUE)[scenarios]
      }
      
      if(naive){
        naive_cmeans = NULL
        for(i in 1:nr_lists){
          naive_cmeans[[i]] = colMeans(EM_betas_list[[i]], na.rm = TRUE)[scenarios]
        }
        
        yline = NULL
        for(i in 1:(nr_scenarios / scenario_rep)){
          yline = c(yline, EM_cmeans[[1]][i], EM_cmeans[[2]][i], naive_cmeans[[1]][i], naive_cmeans[[2]][i])
        }
        hline = data.frame(xline = 1:nr_scenarios, yline = yline)
        
        if(TRUE %in% (hline$yline > y_lim_max)){
          hline$yline[which(hline$yline > y_lim_max)] = y_lim_max
        }
        
      } else {
        yline = NULL
        for(i in 1:(nr_scenarios / scenario_rep)){
#          yline = c(yline, EM_cmeans[[1]][i], EM_cmeans[[2]][i])
          yline = c(yline, unlist(lapply(EM_cmeans, function(z){z[i]})))
        }
        hline = data.frame(xline = 1:nr_scenarios, yline = yline)
        
        if(TRUE %in% (hline$yline > y_lim_max)){
          hline$yline[which(hline$yline > y_lim_max)] = y_lim_max
        }
      }
      
    } else {
      
      EM_cmeans = colMeans(EM_betas, na.rm = TRUE)[scenarios]
      
      if(naive){
        naive_cmeans = colMeans(naive_betas, na.rm = TRUE)[scenarios]
        
        yline = NULL
        for(i in 1:(nr_scenarios / scenario_rep)){
          yline = c(yline, EM_cmeans[i], naive_cmeans[i])
        }
        hline = data.frame(xline = 1:nr_scenarios, yline = yline)
        
        if(TRUE %in% (hline$yline > y_lim_max)){
          hline$yline[which(hline$yline > y_lim_max)] = y_lim_max
        }
        
      } else {
        yline = NULL
        for(i in 1:nr_scenarios){
          yline = c(yline, EM_cmeans[i])
        }
        hline = data.frame(xline = 1:nr_scenarios, yline = yline)
        
        if(TRUE %in% (hline$yline > y_lim_max)){
          hline$yline[which(hline$yline > y_lim_max)] = y_lim_max
        }
      }
    }
  }
  
  
  # if(!is.null(extra_lines)){
  #   values = extra_lines[[1]]
  #   names_values = names(values[[1]]); len_values = length(names_values)
  # 
  #   scenarios = names(plot_order)
  # 
  #   if(naive){
  #     stop("Sorry, this has not been implemented yet\n")
  # 
  #   } else {
  #     hline_extra = data.frame(xline = 1:nr_scenarios, 0, row.names = scenarios)
  # 
  #     if(len_values != 1){
  #       for(i in 2:len_values){
  #         hline_extra = cbind(hline_extra, 0)
  #       }
  #     }
  #     colnames(hline_extra) = c("xline", names_values)
  # 
  # 
  #     for(i in 1:len_values){
  #       name_value = names_values[i]
  # 
  #       for(j in scenarios){
  #         hline_extra[j, name_value] = values[[j]][name_value]
  #       }
  #     }
  #   }
  # }
  ## Above extra lines is tailer made for what I inputted, but the input has changed because it was unlogical.
  ## But what do we want to do now? list in list? probably the easiest for extensions...
  if(!is.null(extra_lines)){
    ##
    for(i in 1:length(extra_lines)){
      EL_name = colnames(extra_lines[[i]])
      extra_lines[[i]] = matrix(extra_lines[[i]][scenarios, ], nrow = nrow(extra_lines[[i]]), dimnames = list(scenarios, EL_name))
    }
    ##
    
    hline_extra = data.frame(xline = 1:nr_scenarios)
    hline_extra$rownames = rep(scenarios, each = scenario_rep)
    
    names_values = colnames(extra_lines[[1]]); len_values = length(names_values)
    # if(len_values == 1){
    #   names_values = unlist(strsplit(names_values, "noY_"))[[2]]
    # }
    
    for(i in 1:len_values){
      hline_extra = cbind(hline_extra, 0)
    }
    names(hline_extra) = c("xline", "rownames", names_values)
    
    for(i in 1:scenario_rep){
      for(j in 1:nrow(extra_lines[[1]])){
        hline_extra[i + scenario_rep * (j - 1), names_values] = extra_lines[[i]][j, ]  
      }
    }
  } else {
    second_Yaxis = NULL
  }
  
  
  

##### Indication lines #####
  distance_one_line = abs(y_lim_min / nr_variables)
  distance = distance_one_line * 0.8
  IndLine_height = distance_one_line * 0.2
  
#  max_vals = -IndLine_height - (distance * 1.5) * 0:(nr_variables - 1)
  max_vals = -IndLine_height - distance_one_line * 0:(nr_variables - 1)
  min_vals = max_vals - distance
  
  text_coord = max_vals - distance * (1 / 4)
  
  line_coord = rep(list(NULL), nr_variables)
  for(i in 1:nr_variables){
    variable = variables[i]
    
    if(variable == "naive"){
      if(lists){
        obs = rep(c(0, 0, 1, 1), nrow(data_lines))
      } else {
        obs = rep(c(0, 1), nrow(data_lines))  # added the /2 because twice as many x-axis
      }

    } else if(variable %in% list_variable){
      if(naive){
        obs = rep(0:(nr_lists - 1) / (nr_lists - 1), nrow(data_lines) * 2)
        
      } else {
        obs = rep(0:(nr_lists - 1) / (nr_lists - 1), nrow(data_lines))
      
    #   if(naive){
    #     obs = rep(c(0, 1, 0, 1), nrow(data_lines))
    #   } else {
    #     obs = rep(c(0, 1), nrow(data_lines))
      }
      
    } else {
      obs = NULL
      for(j in 1:nrow(data_lines)){
        
        obs = c(obs, rep(data_lines[j, variable], scenario_rep))
        
        # if(lists){
        #   if(naive){
        #     obs = c(obs, rep(data_lines[j, variable], 4))
        #   } else {
        #     obs = c(obs, rep(data_lines[j, variable], 2))
        #   }
        # } else {
        #   if(naive){
        #     obs = c(obs, rep(data_lines[j, variable], 2))
        #   } else {
        #     obs = c(obs, rep(data_lines[j, variable], 1))
        #   }
        # }
      }
    }
    len_obs = length(obs)
    
    un_obs = sort(unique(obs))
    len_un_obs = length(un_obs)
    
    # Calculation x-coordinates
    x = 1
    for(j in 2:len_obs){
      x = c(x, j)
      
      if(obs[(j - 1)] != obs[j]){
        x = c(x, j)
      }
    }
    
    x = c(x, (len_obs + 1))
    len_x = length(x)
    
    
    # Calculation y-coordinates
    diff = (max_vals[i] - min_vals[i]) / (len_un_obs - 1)
    
    y = curr_val = min_vals[i]
    for(j in 1:(len_x - 1)){
      if(x[j] == x[(j + 1)]){
        curr_val = curr_val + diff * (which(un_obs == obs[x[j]]) - which(un_obs == obs[x[(j - 1)]]))
      }
      
      y = c(y, curr_val)
    }
    
    line_coord[[i]] = list(x = x, y = y)
  }
  
  
  
  
  

##### ggplot #####
  p = ggplot(data = plot_dat1, mapping = aes(x = Model, y = RMSE)) +
    
    theme(
      panel.grid.minor.x = element_blank(), # get rid of minor grid
      panel.grid.major.x = element_blank(),

      panel.grid.minor.y = element_blank(), # get rid of minor grid
      panel.grid.major.y = element_blank(),
      
      legend.title = element_blank(),
      legend.text = element_text(size = 8),
      legend.position = "right",
      
      axis.title.x = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank()
    ) +
    
    scale_color_discrete(guide = guide_legend(ncol = 1))
  
  
## Plot title
  if(!is.null(title)){
    p = p + 
      labs(title = title)
  }

  

  ## Horizontal lines  
  gg_breaks = as.list(ggplot_build(p)$layout$panel_params[[1]]$y.sec)$breaks
  gg_breaks = gg_breaks[!is.na(gg_breaks)]
  gg_labels = gg_breaks
  
  y_range = range(gg_breaks, na.rm = TRUE)
  
  if(is.null(second_Yaxis)){
    p = p + 
      scale_y_continuous(breaks = gg_breaks, labels = gg_labels, name = yname)
  } else {
    p = p + 
      scale_y_continuous(breaks = gg_breaks, labels = gg_labels, name = yname,
                         sec.axis = sec_axis(trans = ~.*second_Yaxis, name = names(second_Yaxis), breaks = gg_breaks * second_Yaxis))
  }
  

  gg_labels = gg_breaks[!is.na(gg_breaks)]
  len_gg_labels = length(gg_labels)
  for(i in 1:len_gg_labels){
    p = p + geom_hline(yintercept = gg_labels[i], col = "white")
  }
  
  
  
## Vertical lines  
  len_obs_p1 = len_obs + 1
  segment_data = data.frame(x = 1:len_obs_p1, xend = 1:len_obs_p1, y = rep(0, len_obs_p1), yend = rep(Inf, len_obs_p1))
  p = p + geom_segment(data = segment_data, mapping = aes(x = x, y = y, xend = xend, yend = yend), col = "white")
  
  

    
## Indication lines
  for(i in 1:nr_variables){
    p = p +   
      geom_line(data = data.frame(line_coord[[i]]), mapping = aes(x = x, y = y)) +
#      geom_text(label = variable_labels[i], x = -2.5, y = text_coord[i], hjust = "left", size = 3)
      
      geom_text(label = variable_labels1[i], x = max(-2.5, 0.5 - nr_scenarios * 0.025), y = text_coord[i], hjust = "left", size = 3) +
      geom_text(label = variable_labels2[i], x = max(-2.4, 0.6 - nr_scenarios * 0.025), y = text_coord[i] - (distance * (2 / 5)), hjust = "left", size = 3) 
  }
  
  
  
  
## RMSE Values
  nr_haplos = nr_haplos_remain = length(unique(c(unique(plot_dat1$Haplo), unique(plot_dat2$Haplo))))
  nr_cols = 1
  
  while(nr_haplos_remain > 35){
    nr_haplos_remain = nr_haplos_remain - 35
    nr_cols = nr_cols + 1
  }
  
  if(nr_cols == 2 & nr_haplos_remain <= 4){
    nr_cols = 1
  }

  
  
  p = p + geom_jitter(aes(col = Haplo), width = 0.25, size = 2.5) +
    scale_color_manual(values = hue_pal()(nr_haplos),
                       guide = guide_legend(ncol = nr_cols))
  
  if(nrow(plot_dat2) != 0){
    plot_dat2$RMSE = y_lim_max
    
    #plot_dat2$col = unique(ggplot_build(p)$data[[13]]["colour"])
    
    p = p + geom_jitter(data = plot_dat2, mapping = aes(x = Model, y = RMSE, colour = Haplo), size = 2.5, shape = 2, show.legend = FALSE, height = 0, width = 0.25) +  # If show.legend = TRUE, empty triangles will be added to legend
      scale_color_manual(values = hue_pal()(nr_haplos),
                         guide = guide_legend(ncol = nr_cols))
  } 
  
  
  
  
    


## Horizontal mean/extra lines
  if(mean_plot | !is.null(extra_lines)){
    df_all = line_colours = line_colours_names = NULL; nr_lines = 0; EL_legend = FALSE
    
#    cbbPalette <- c("#000000", "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
    palette_names = c("noY_Al", "Y_alleles_Al", "Y_backward_Al", "Y_penalized_Al", "naive_MV_noY_Al", "naive_Dict_Al", "naive_Demo_Al", "naive_MV_Y_Al")
    palette = brewer.pal(length(palette_names), "Paired")  #; display.brewer.pal(length(palette_names), "Paired")
    
    if(mean_plot){
      df_all = cbind(hline, "group" = mean_name)
      
      nr_lines = nr_lines + 1
      line_colours = "#000000"
      line_colours_names = mean_name
    }
    
    if(!is.null(extra_lines)){
      if(len_values == 1 & !is.null(HTR_name)){
        df_all = rbind(df_all, cbind("xline" = hline_extra$xline, "yline" = hline_extra[, names_values[1]], "group" = names(HTR_name)))
        
        line_colours = c(line_colours, as.character(HTR_name)); line_colours_names = c(line_colours_names, names(HTR_name))
        nr_lines = nr_lines + 1
        
      } else {
        for(i in 1:len_values){
          df_all = rbind(df_all, cbind("xline" = hline_extra$xline, "yline" = hline_extra[, names_values[i]], "group" = names_values[i]))
          
          line_colours = c(line_colours, palette[str_detect(palette_names, names_values[i]) & ifelse(rep(naive, length(palette_names)), str_detect(palette_names, "naive"), !str_detect(palette_names, "naive"))])
          
          ##
          # if(length(palette[palette_names == names_values[i]]) == 0){
          #   if(names_values[i] == "Al"){
          #     line_colours = c(line_colours, palette[6])
          #   }
          # }
          ##
          
          nr_lines = nr_lines + 1
        }
        
        # line_colours = c(line_colours, palette[1:len_values])
        line_colours_names = c(line_colours_names, names_values)
      }
      

      
      EL_legend = NA 
    }
    
    df_all = as.data.frame(df_all)
    df_all$xline = as.numeric(df_all$xline)
    df_all$yline = as.numeric(df_all$yline)
    
    names(line_colours) = line_colours_names
    
    
    p = p +
      
      new_scale_color() +
#      new_scale("shape") +
      
      geom_segment(data = df_all, aes(x = xline - 0.02, xend = xline + 1.02, y = yline, yend = yline, linetype = group, col = group), size = 1.5, show.legend = EL_legend) + 
      
      scale_linetype_manual(values = rep(1, nr_lines), labels = line_colours_names) +
      scale_color_manual(values = line_colours, labels = line_colours_names) +
      
#      scale_shape_manual(values = c(0, 0, 0)) +
        
      theme(legend.title = element_blank())
    
  }
  
  
  ##
  if(!is.null(altered_replications)){
    altered_replications = altered_replications[scenarios]

    p = p +
      geom_text(label = "Alt. Repl.:", x = max(-2.5, 0.5 - nr_scenarios * 0.025), y = (max_vals[1] - distance) * (1 / 16), size = 3, hjust = "left")
    
    for(i in 1:nrow(data_lines)){
#      if(altered_replications[i] != 0){
        p = p + 
          geom_text(label = altered_replications[i], x = data_lines[i, "Model"] + nr_lists - 0.5, y = (max_vals[1] - distance) * (1 / 16), size = 3)
#      }
    }
  }
  ##
  

  #  p$layers = c(geom_boxplot(), p$layers)  Can change layers; not needed here
  
  # attributes(ggplot_build(p))
  # as_ggplot(get_legend(p))

    
  legend_width = as.numeric(get_legend(p)$widths[3])
  
  return(list("p" = p, "plot_order" = scenarios, "legend_width" = legend_width))
}
