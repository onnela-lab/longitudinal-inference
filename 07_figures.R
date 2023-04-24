#################################################################################################### 
# This script generates plots for our paper
#################################################################################################### 
library(ggplot2)
library(gridExtra)
library(cowplot)
library(stringr)
library(reshape2)

## Network Plots
network_plots <- function(){
  
  # Parameters for the Network
  mu <- 0.000
  rho <- 0.3
  sigma <- 0.1
  w_0 <- 0.4
  w_1 <- 0.2
  g_union <- c(1, 6, 12)
  
  param_all <- c("mu", "rho", "sigma",  "w_0", "w_1", "g_union", "lag") 
  
  # Creating the graphs
  graphs <- list()
  edge_tracker <- list()
  edge_color <- list()
  edge_type <- list()
  
  for(i in 1:1){
    temp_graphs <- list()
    temp_edge_tracker <- list()
    temp_edge_color <- list()
    temp_edge_type<- list()
    for(j in 1:3){
      print(c(i,j))
      results <- Sweden_MSM(num_nodes = 250, 
                            mu = mu[i],
                            rho = rho[i],
                            sigma = sigma[i], 
                            w_0 = w_0[i],
                            w_1 = w_1[i],
                            iter = 500,
                            g_union = g_union[j],
                            distribution = "discrete",
                            lag = 15)
      
      temp_graphs[[j]] <- results[[2]]
      
      temp_graphs[[j]] <- simplify(temp_graphs[[j]])
      
      temp_edge_tracker[[j]] <- results[[3]]
      
      last_iteration <- max(as.numeric(temp_edge_tracker[[j]][,"iteration"]))
      
      cas_edges <- which(as.numeric(temp_edge_tracker[[j]][, "start"]) >= last_iteration - g_union[j] + 1 &
                           temp_edge_tracker[[j]][, "type"] == "c")
      
      steady_edges <- which(as.numeric(temp_edge_tracker[[j]][, "end"]) > last_iteration - g_union[j] + 1  &
                              temp_edge_tracker[[j]][, "type"] == "s")
      
      
      e_color <- NULL
      e_type <- NULL
      
      ce <- c(t(temp_edge_tracker[[j]][cas_edges, c("e1", "e2")]))
      ce_i <- edge_index(temp_graphs[[j]], ce)
      e_color[ce_i] <- "blue"
      e_type[ce_i] <- 1
      
      # All repeats default to steady
      se <- c(t(temp_edge_tracker[[j]][steady_edges, c("e1", "e2")]))
      se_i <- edge_index(temp_graphs[[j]], se)
      e_color[se_i] <- "red"
      e_type[se_i] <- 2
      
      if(any(is.na(e_type))){
        nn <- which(is.na(e_type))
        e_color[nn] <- "red"
        e_type[nn] <- 2
        print(length(nn))
      }
      
      temp_edge_color[[j]] <- e_color
      temp_edge_type[[j]] <- e_type
    }
    
    graphs[[i]] <- temp_graphs
    edge_tracker[[i]] <- temp_edge_tracker
    edge_color[[i]] <- temp_edge_color
    edge_type[[i]] <- temp_edge_type
  }
  
  # Plotting the graphs 
  
  par(mfrow = c(1,3))
  for(i in 1:1){
    for(j in 1:3){
      title <- paste("mu = ", mu[i],
                     ", rho = ", rho[i],
                     ", sigma = ", sigma[i],
                     ",\n w_0 = ", w_0[i],
                     ", w_1 = ", w_1[i],
                     ", aggregation = ", g_union[j], sep = "")
      
      plot(graphs[[i]][[j]], 
           vertex.size = 0.5,
           vertex.label = NA,
           edge.color = edge_color[[i]][[j]],
           edge.lty = edge_type[[i]][[j]],
           edge.width = .15)
      
      #title(title, cex.main=1)
    }
  }
  
  
}

################################################################################
################################################################################

## Mapping Functions 
plots_csv <- function(){
  
  population <- c("Sweden")
  request <- c("num_singletons", "steady_str_end_len","two_concurrency", "steady_num")
  param_all <- c("rho", "sigma",  "w_0", "w_1")
  d_c <- c("discrete_free", "discrete_fixed")
  
  for(d_c_i in d_c){
    all_graphs <- list()
    for(request_i in request){
      for(z in population){
        for(i in param_all){
          name_i <- i
          print(c(z, request_i, i))
          wd_all <- file.path("/Volumes/SanDisk/Research/JP/HIV/Data/figures")
          setwd(wd_all)
          
          # Reading in Data
          #dim_all <- as.data.frame(read.csv(paste(i, request_i ,"mapping_data.csv", sep = "_")))
          dim_all <- as.data.frame(read.csv(paste(i, request_i ,d_c_i, "mapping_data.csv", sep = "_")))
          
          # Correcting column names
          c_names <- str_replace(colnames(dim_all), "X", "")
          for(k in 1:length(c_names)){
            if(nchar(c_names[k]) < 5){
              if(c_names[k] == 1){
                c_names[k] <- "1.000"
              }else{
                c_names[k] <- paste(c_names[k], 
                                    rep(0, 5 - nchar(c_names[k]) ),
                                    sep = "")
              }
            }
          }
          
          colnames(dim_all) <- substr(c_names, 1, 4)
          dim_all <- dim_all[, c(seq(2, 30, 3))]
          
          
          # Creating ggplot
          #title_1 <- paste(z, "Summary Statistic Plot")
          title_1 <- NULL
          if(d_c_i == "discrete_free"){
            if(request_i == "num_singletons"){
              #ylab_1 <- "Proportion of Nodes \n with 0 Partners"
              ylab_1 <- parse(text = "s[1]")
              y_lim <- c(0,1)
            }else{
              if(request_i == "steady_str_end_len"){
                #ylab_1 <- "Average Duration of \n Relationships During \n Study"
                ylab_1 <- parse(text = "s[2]")
                y_lim <- c(0,30)
              }else{
                if(request_i == "two_concurrency"){
                  #ylab_1 <- "Proportion of \n CasualRelationshiops \n aslo in a Steady \n Relationship"
                  ylab_1 <- parse(text = "s[3]")
                  y_lim <- c(0,1)
                }else{
                  if(request_i == "steady_num"){
                    #ylab_1 <- "Proportion of \n Relationships that \n are Steady"
                    ylab_1 <- parse(text = "s[4]")
                    y_lim <- c(0,1)
                  }
                }
              }
            }
          }else{
            if(request_i == "num_singletons"){
              #ylab_1 <- "Proportion of Nodes \n with 0 Partners"
              ylab_1 <- parse(text = "s[1]")
              y_lim <- c(0,0.8)
            }else{
              if(request_i == "steady_str_end_len"){
                #ylab_1 <- "Average Duration of \n Relationships During \n Study"
                ylab_1 <- parse(text = "s[2]")
                y_lim <- c(0,15)
              }else{
                if(request_i == "two_concurrency"){
                  #ylab_1 <- "Proportion of \n CasualRelationshiops \n aslo in a Steady \n Relationship"
                  ylab_1 <- parse(text = "s[3]")
                  y_lim <- c(0.1,1)
                }else{
                  if(request_i == "steady_num"){
                    #ylab_1 <- "Proportion of \n Relationships that \n are Steady"
                    ylab_1 <- parse(text = "s[4]")
                    y_lim <- c(0.2,1)
                  }
                }
              }
            }
          }
          
          if(i == "rho"){
            ylab_1 <- ylab_1
          }else{
            ylab_1 <- ""
          }
          
          if(i == "rho"){
            x_lab <-  parse(text = i)
          }else{
            if(i == "sigma"){
              x_lab <- parse(text = i)
            }else{
              if(i == "w_0"){
                x_lab <- parse(text = "omega[0]")
              }else{
                if(i == "w_1"){
                  x_lab <- parse(text = "omega[1]")
                }
              }
            }
          }
          
          if(request_i == "steady_num"){
            if(i == "rho"){
              p1 <- ggplot(stack(dim_all), aes(x = ind, y = values)) +
                theme_bw() +
                geom_boxplot() + 
                xlab(x_lab) +
                ylab(ylab_1) + 
                ylim(y_lim[1], y_lim[2]) +
                theme(axis.text.x = element_text(angle = 90, size = 18),
                      axis.text.y = element_text(size = 18), 
                      text = element_text(size = 24))
            }else{
              p1 <- ggplot(stack(dim_all), aes(x = ind, y = values)) +
                theme_bw() +
                geom_boxplot() + 
                xlab(x_lab) +
                ylab(ylab_1) + 
                ylim(y_lim[1], y_lim[2]) +
                theme(axis.text.x = element_text(angle = 90, size = 18),
                      axis.text.y = element_blank(), 
                      text = element_text(size = 24))
            }
          }else{
            if(i == "rho"){
              p1 <- ggplot(stack(dim_all), aes(x = ind, y = values)) +
                theme_bw() +
                geom_boxplot() + 
                xlab("") +
                ylab(ylab_1) + 
                ylim(y_lim[1], y_lim[2]) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_text(size = 18), text = element_text(size = 24))
            }else{
              p1 <- ggplot(stack(dim_all), aes(x = ind, y = values)) +
                theme_bw() +
                geom_boxplot() + 
                xlab("") +
                ylab(ylab_1) + 
                ylim(y_lim[1], y_lim[2]) +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank())
            }
          }
          
          
          #plot(p1)
          
          all_graphs[[request_i]][[name_i]] <- p1
        }
        
      }
      
    }
    grid.arrange(grobs = align_plots(plotlist = list(all_graphs[[1]][[1]], all_graphs[[1]][[2]], all_graphs[[1]][[3]], all_graphs[[1]][[4]],
                                                     all_graphs[[2]][[1]], all_graphs[[2]][[2]], all_graphs[[2]][[3]], all_graphs[[2]][[4]],
                                                     all_graphs[[3]][[1]], all_graphs[[3]][[2]], all_graphs[[3]][[3]], all_graphs[[3]][[4]], 
                                                     all_graphs[[4]][[1]], all_graphs[[4]][[2]], all_graphs[[4]][[3]], all_graphs[[4]][[4]]),
                                     align = "hv"), nrow = 4, ncol = length(request))
  }
  
  #return(all_graphs)
  }
  
plots_csv()
  
################################################################################
################################################################################

## Posterior plots
final_function <- function(param_matrix,
                           reg_sum_data,
                           dist_types,
                           g_union_check,
                           lag, 
                           z){
  
  nn <- which(complete.cases(reg_sum_data))
  reg_sum_data <- reg_sum_data[nn, ]
  param_matrix <- param_matrix[nn, ]
  
  
  
  # Standardized summary data
  std_sum_data <- scale(reg_sum_data)
  
  
  # Creating the distance columns (different metrics will be used)
  dist_matrix <- matrix(NA, dim(std_sum_data)[1], length(dist_types))
  colnames(dist_matrix) <- dist_types
  
  # Contribution looks at which variables contributes to distance 
  # (percentage explained as well in projections)
  
  L2_calc <- function(x){
    # calculating L2 norm
    
    # x
    # row of the summary matrix
    # 1st row is observed data
    x <- x - std_sum_data[1, ]
    L2 <- norm(x, type = "2")
    
    return(L2)
    
  }
  
  dist_matrix[, dist_types[k]] <- apply(std_sum_data, # use standardized data for L2
                                        MARGIN = 1, 
                                        FUN = L2_calc) 
  
  contribution <- request
  
  final <- cbind(param_matrix, 
                 reg_sum_data)
  final <- cbind(final, dist_matrix)
  
  return(list(final, 
              request))
  
}

ABC_ref <- function(request,
                    population,
                    param_all,
                    num_samples,
                    dist_types,
                    quant,
                    g_union_check,
                    stat_table = FALSE,
                    reg_adj = FALSE,
                    return_grid = FALSE,
                    cheat = TRUE,
                    three = FALSE,
                    lag = 0,
                    prior = FALSE,
                    final_lag = FALSE){
  # request:
  # Summary statistics we are requesting
  # num_edges: number of edges
  # edges_m: proportion of edges maintained between time points
  # num_tri: number of triangles in the graph
  # quant_25: 25th degree distribution percentile
  # quant_50: 50th degree distribution percentile
  # quant_75: 75th degree distribution percentile
  # quant_max: maximum of the degree distribution 
  # quant_mean: mean of the degree distribution 
  # median_local_cc: median of the local clustering coefficient
  # mean_local_cc: mean of the local clustering coefficient 
  # num_con_comp: number of components with 2 or more nodes
  # size_LCC: size of the largest connected component
  # num_singletons: number of nodes with degree 0
  # num_exclusive_pairs: number of monogamous pairs
  # prolonged_single:
  # entering
  # population:
  # Country/City we are looking at
  # Sweden
  # Amsterdam
  # param_all
  # parameters to be considered (each model will be different)
  # num_samples
  # number of samples we are considering from the prior distribution
  # dist_types:
  # Distances to be included in the ABC reference table
  # L2
  # L2 norm
  # quant:
  # percentile to threshold the distance metric
  # g_union_check:
  # which union iteration to check
  # stat_table:
  # TRUE: we already calcualted the stat table
  # FALSE: we still need to calculate stat table
  # reg_adj:
  # Regression adjustment
  # return_grid
  # Return ggplots from posterior
  # three:
  # If true calculate 3 different observed versions of the graph
  # results
  # prior:
  # Wheter you investigate the prior instead of the posterior
  
  for(z in population){
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    # Generating parameter matrix
    file <- paste(z, "param_table.csv", sep = "_")
    param_matrix <- read.csv(file)
    
    
    # Generating reference table
    # Saving time and reading long portion of data in
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    file_name <- paste(z, "ABC_stat_table.csv", sep = "_")
    
    reg_sum_data <- read.csv(file_name)
    reg_sum_data <- reg_sum_data[, request, drop = FALSE] # removing extra column
    
    
    nn <- sample(which(param_matrix[, "g_union"] == g_union_check &
                         param_matrix[, "lag"] == lag))
    
    reg_sum_data <- reg_sum_data[nn, , drop = FALSE]
    param_matrix <- param_matrix[nn, , drop = FALSE]
    
    nw <- which(param_matrix[, "rho"] > .13 & param_matrix[, "rho"] < .14 &
                  param_matrix[, "sigma"] > .24 & 
                  param_matrix[, "w_0"] > .2 & 
                  param_matrix[, "w_1"] > .1)
    
    ni <- nw
    
    nn <- c(ni, 1:I(ni-1), I(ni+1):dim(reg_sum_data)[1])
    reg_sum_data <- reg_sum_data[nn, ,drop = FALSE]
    param_matrix <- param_matrix[nn, , drop = FALSE]
    
    # observed data
    o_data <- cbind(param_matrix[1,,drop = FALSE], reg_sum_data[1,, drop = FALSE])
    
    
    for(i_o in 1:dim(o_data)[1]){
      final_full <- final_function(param_matrix = param_matrix,
                                   reg_sum_data = reg_sum_data,
                                   dist_types = dist_types,
                                   g_union_check = g_union_check,
                                   lag = lag,
                                   z = z)
      
      if(prior){
        final <- cbind(param_matrix, reg_sum_data)
        final <- cbind(final, NA)
      }else{
        final <- final_full[[1]]
        contribution <- final_full[[2]]
      }
      
      
      
      # Creating Posterior Plots
      post_plots <- function(full_data, dist_k, quant_x, prior = FALSE){
        # full_data:
        # data to be used to obtain the posterior
        # dist_k:
        # particular distances to use when creating the posterior 
        # plot
        # quant_x:  
        # percentile to check when selecting samples
        
        
        # identifying data
        full_sim_data <- full_data
        full_sim_data <- full_sim_data[c(1, sample(1:dim(full_sim_data)[1], 10000)), ] # reference table
        full_dist_data <- full_sim_data[-1, dist_k]
        quant_threshold <- quantile(full_dist_data, probs = quant_x, na.rm = TRUE)
        print(quant_threshold)
        if(prior){
          post_data <- full_sim_data
        }else{
          post_data <- full_sim_data[which(full_dist_data <= quant_threshold) + 1, ,drop = FALSE ] #+1 resets order
        }
        
        # storing the plots
        p_length <- 4
        
        
        
        # Regression adjustment
        if(reg_adj){
          
          # Adjusting posterior
          for(k in 1:length(param_all)){
            #print(param_all[k])
            formula <- paste(param_all[k], "~", paste(request, collapse = "+"))
            lm1 <- lm(formula = formula, data = as.data.frame(post_data))
            
            
            #define weights to use
            #wt <- 1 / lm(abs(lm1$residuals) ~ lm1$fitted.values)$fitted.values^2
            
            #perform weighted least squares regression
            #wls_model <- lm(formula = formula, data = post_data, weights = wt)
            
            # observed estimate
            obs_est <- predict(lm1, data.frame(full_data[1, request, drop = FALSE]))
            
            # posterior residuals (neighborhood within s_obs)
            post_resid <- resid(lm1)
            
            # Adjustment
            theata_c <- obs_est + post_resid
            
            print(c(sum(theata_c < 0), param_all[k]))
            
            # Making negative predicted values 0
            theata_c <- ifelse(theata_c < 0, 0, theata_c)
            
            # sub
            post_data[, param_all[k]] <- theata_c
          }
        }
        
        all_ggplots = NULL
        # Making plots
        for(k in 1:length(param_all[1:4])){
          if(param_all[k] == "gamma"){
            x_lim = c(1.5, 4)
            xlab = expression(gamma)
          }else{
            if(param_all[k] == "p_0"){
              x_lim = c(0, 0.8)
              xlab = param_all[k] 
            }else{
              if(param_all[k] == "p_cas"){
                x_lim = c(0, 1)
                xlab = param_all[k]
              }else{
                if(param_all[k] == "k_max"){
                  x_lim = c(5, 100)
                  xlab = param_all[k]
                }else{
                  if(param_all[k] == "migration"){
                    x_lim = c(0, 0.2)
                    xlab = param_all[k]
                  }else{
                    if(param_all[k] == "mu"){
                      x_lim = c(0, 0.5)
                      xlab = expression(mu)
                    }else{
                      if(param_all[k] == "rho"){
                        x_lim = c(0, 1)
                        xlab = parse(text = "rho")
                      }else{
                        if(param_all[k] == "sigma"){
                          x_lim = c(0, 1)
                          xlab = parse(text = "sigma")
                        }else{
                          if(param_all[k] == "w_0"){
                            x_lim = c(0, 1)
                            xlab = parse(text = "omega[0]")
                          }else{
                            if(param_all[k] == "w_1"){
                              x_lim = c(0, 1)
                              xlab = parse(text = "omega[1]")
                            }else{
                              if(param_all[k] == "theta"){
                                x_lim <- c(0, 1)
                                xlab = param_all[k]
                              }
                            }
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
          
          
          if(k == 1){
            y_lab = "Density"
          }else{
            y_lab = ""
          }
          
          # Making the plots
          
          if(param_all[k] == "rho"){
            if(final_lag){
              all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
                geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
                #scale_y_continuous(labels = scales::percent) + 
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, size = 18),
                      axis.text.y = element_text(size = 18), 
                      text = element_text(size = 24)) +
                xlab(xlab)+
                ylab(y_lab) +
                ylim(0, 30) +
                xlim(0, 1) +
                #ylim(c(0, dim(post_data)[1])) + 
                geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.05, 
                                                 na.rm = TRUE),
                           col="blue") + 
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.95, 
                                                 na.rm = TRUE),
                           col="blue")
            }else{
              all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
                geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
                #scale_y_continuous(labels = scales::percent) + 
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_text(size = 18), 
                      text = element_text(size = 24)) +
                xlab("")+
                ylab(y_lab) +
                ylim(0, 30) +
                xlim(0, 1) +
                #ylim(c(0, dim(post_data)[1])) + 
                geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.05, 
                                                 na.rm = TRUE),
                           col="blue") + 
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.95, 
                                                 na.rm = TRUE),
                           col="blue")
            }
          }else{
            if(final_lag){
              all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
                geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
                #scale_y_continuous(labels = scales::percent) + 
                theme_bw() +
                theme(axis.text.x = element_text(angle = 90, size = 18),
                      axis.text.y = element_blank(), 
                      text = element_text(size = 24)) +
                xlab(xlab)+
                ylab("") +
                ylim(0, 30) +
                xlim(0, 1) +
                #ylim(c(0, dim(post_data)[1])) + 
                geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.05, 
                                                 na.rm = TRUE),
                           col="blue") + 
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.95, 
                                                 na.rm = TRUE),
                           col="blue")
            }else{
              all_ggplots[[k]] <- p <- ggplot(as.data.frame(post_data), aes_string(x = param_all[k]))+
                geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
                #scale_y_continuous(labels = scales::percent) + 
                theme_bw() +
                theme(axis.text.x = element_blank(),
                      axis.text.y = element_blank(), 
                      text = element_text(size = 24)) +
                xlab("")+
                ylab("") +
                ylim(0, 30) +
                xlim(0, 1) +
                #ylim(c(0, dim(post_data)[1])) + 
                geom_vline(xintercept = full_data[1, param_all[k]], col = "red", linetype="dotted") + #observed
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.05, 
                                                 na.rm = TRUE),
                           col="blue") + 
                geom_vline(xintercept = quantile(post_data[, param_all[k]], 
                                                 probs = 0.95, 
                                                 na.rm = TRUE),
                           col="blue")
            }
          }
        }
        
        return(all_ggplots)
      
      }
      
      
      return(post_plots(full_data = final,
                        dist_k = dist_types,
                        quant_x = quant,
                        prior = prior))
      
    }
    
    
    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/ref_table")
    file_name <- paste(z, "ABC_ref_table.csv", sep = "_")
    
    write.csv(final, file = file_name)
    
  }
  
  #return(final)
  
}


results_prior <- ABC_ref(request = c("num_singletons", "steady_str_end_len", "two_concurrency", "steady_num"),
                         population = c("Sweden"),
                         param_all = c("rho", "sigma", "w_0", "w_1"),
                         dist_types = "No_Selection",
                         quant = 0.01,
                         g_union_check = 1,
                         stat_table = TRUE, 
                         reg_adj = FALSE,
                         return_grid = FALSE,
                         cheat = FALSE,
                         three = FALSE,
                         lag = -1,
                         prior = TRUE)

results_no_lag <- ABC_ref(request = c("num_singletons", "steady_str_end_len", "two_concurrency", "steady_num"),
                         population = c("Sweden"),
                         param_all = c("rho", "sigma", "w_0", "w_1"),
                         dist_types = "No_Selection",
                         quant = 0.01,
                         g_union_check = 1,
                         stat_table = TRUE, 
                         reg_adj = TRUE,
                         return_grid = FALSE,
                         cheat = FALSE,
                         three = FALSE,
                         lag = -1,
                         prior = FALSE)

results_50 <- ABC_ref(request = c("num_singletons", "steady_str_end_len", "two_concurrency", "steady_num"),
                          population = c("Sweden"),
                          param_all = c("rho", "sigma", "w_0", "w_1"),
                          dist_types = "No_Selection",
                          quant = 0.01,
                          g_union_check = 1,
                          stat_table = TRUE, 
                          reg_adj = TRUE,
                          return_grid = FALSE,
                          cheat = FALSE,
                          three = FALSE,
                          lag = 50,
                          prior = FALSE,
                      final_lag = TRUE)

plot_grid(plotlist = list(results_prior[[1]],results_prior[[2]],results_prior[[3]],results_prior[[4]],
                          results_no_lag[[1]],results_no_lag[[2]],results_no_lag[[3]],results_no_lag[[4]],
                          results_50[[1]],results_50[[2]],results_50[[3]],results_50[[4]]),
          align = "hv",
          ncol = 4,
          nrow = 3)



################################################################################
################################################################################

# Prior Plots
prior_plots <- function(){
  # Sampling from the prior
  s <- runif(10000, 1, 90)
  sigma <- 1/s

  r <- runif(10000, 1, 50)
  rho <- 1/r
  
  w1 <- runif(10000, 1, 61)
  w_1 <- 1/w1
  
  w0 <- runif(10000, 1, 40)
  w_0 <- 1/w0
  
  mu <- 0.000
  
  df = data.frame(s, sigma, r, rho, w1, w_1, w0, w_0)
  
  
  p1 <-ggplot(df, aes(x=r)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "1/rho")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_text(size = 18)) +
    theme(text=element_text(size=24)) + 
    ylab("Density") +
    ylim(0, 0.04) #+
    #xlim(0, 95)
  
  p2 <-ggplot(df, aes(x=s)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "1/sigma")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24)) + 
    ylab("") +
    ylim(0, 0.04) #+
    #xlim(0, 55)
  p3 <-ggplot(df, aes(x=w0)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "1/omega[0]")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24))+ 
    ylab("") +
    #xlim(c(0, 42)) #+
    ylim(0, 0.04)
  p4 <-ggplot(df, aes(x=w1)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "1/omega[1]")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24))+ 
    ylab("")+
    ylim(0, 0.04) #+
    #xlim(0, 65)
  
  p5 <-ggplot(df, aes(x=rho)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "rho")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_text(size = 18)) +
    theme(text=element_text(size=24)) + 
    ylab("Density")+
    ylim(0, 20)+
    xlim(0, 1) 
  
  p6 <-ggplot(df, aes(x=sigma)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "sigma")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24))+ 
    ylab("")+
    ylim(0, 20)+
    xlim(0, 1) 
  
  p7 <-ggplot(df, aes(x=w_0)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "omega[0]")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24))+ 
    ylab("")+
    ylim(0, 20) +
    xlim(0, 1) 
    
  
  p8 <-ggplot(df, aes(x=w_1)) +
    geom_histogram(aes(y=..density..), color="black", fill="white", bins = 20, boundary = 1) +
    xlab(parse(text = "omega[1]")) + 
    theme_bw() +
    theme(axis.text.x = element_text(angle = 90, size = 18),
          axis.text.y = element_blank()) +
    theme(text=element_text(size=24))+ 
    ylab("")+
    ylim(0, 20) +
    xlim(0, 1) 
  
  plot_grid(plotlist = list(p1, p2, p3, p4, p5, p6, p7, p8),
            align = "hv",
            nrow = 2,
            ncol = 4)
  

}

prior_plots()

################################################################################
################################################################################

# RMSE Plots
ABC_lag_plots <- function(request,
                          population,
                          param_all,
                          dist_types,
                          quant,
                          g_union_check,
                          stat_table = TRUE,
                          root,
                          euclidean){ # bins isnt actually needed 
  
  for(z in population){
    #keep

    setwd("/Volumes/SanDisk/Research/JP/HIV/Data/figures")
    # # Generating discrepancy from the prior
    prior_point <- 2.219386 #Non- reg_adjusted
    
    if(root){ #euclidean must be false if root is true
      if(!euclidean){
        lag_data <- as.data.frame(read.csv(paste("Root_Euclidean_", 
                                                 euclidean,
                                                 "_RMSE_data.csv", sep = "")))
      }else{
        print("error")
      }
    }else{
      lag_data <- as.data.frame(read.csv(paste("Euclidean_", 
                                               euclidean,
                                               "_RMSE_data.csv", sep = "")))
    }

    lag_data$group <- as.factor(lag_data$group)
    lag_data$lag <- lag_data$lag + 1
    
    if(euclidean){
        p <- ggplot(lag_data, aes(x=lag, y=mean, group = group, colour = group)) + 
        geom_point(aes(size = 5, shape = group)) +
        ylab("Average Discrepancy") +
        xlab("Observation Lag") +
        theme_bw() +
        geom_smooth(method = "loess", size = 1.5) + 
        theme(axis.text.x = element_text(angle = 90, size = 18),
              axis.text.y = element_text(size = 18),
              text = element_text(size = 18)) +
        theme(legend.position="top") +
        geom_hline(yintercept= prior_point, linetype="dashed", color = "red") + 
        #geom_errorbar(aes(ymin=mean-1.96*s.e, ymax=mean+1.96*s.e), width=.1) +
        scale_color_manual(labels = c("No Regression Adjustment", "Regression Adjustment"),
                           values = c("blue", "red"), name = NULL) + 
        annotate("text", 16, 2.0, vjust = -1, label = "Prior Average Error", color = "red") + 
          scale_shape_manual(values=c(16, 15))
        
        # Calculating percision gains
        
        #bb <- ggplot_build(p)
        
        #bb_d <- bb$data[[2]]
        
        #bb_d_1 <- bb_d[which(round(bb_d[, "x"]) == 50 & bb_d[, "group"] == 1), ]
        #bb_d_2 <- bb_d[which(round(bb_d[, "x"]) == 50 & bb_d[, "group"] == 2), ]
        
        #(1 - bb_d_2[, "y"]/bb_d_1[, "y"])*100
        #(1 - bb_d_2[, "y"]/prior_point)*100
        
      plot(p)
    }else{
      lag_data <- lag_data[, -which(colnames(lag_data) == "mean_mu")]
      lag_data <- lag_data[which(lag_data[, "group"] == 2), ]
      mdf <- melt(lag_data,id.vars= c("lag", "group"))
      
      p <- ggplot(mdf, aes( x=lag, y=value, colour=variable, group= variable)) + 
        geom_point(aes(size = 5, shape = variable))  + 
        geom_smooth(method = "loess", size = 0.9) +
        ylab("Average Discrepancy") +
        xlab("Observation Lag") +
        theme_bw() +
        theme(axis.text.x = element_text(angle = 90, size = 18),
              axis.text.y = element_text(size = 18),
              text = element_text(size = 18)) +
        theme(legend.position="top") +
        #geom_hline(yintercept= prior_point, linetype="dashed", color = "red") + 
        scale_color_manual(labels = expression(rho,sigma, omega[0], omega[1],
                                      "Total"), name = NULL,
                           values = c("blue", "red", "green", "grey","black")) + 
        scale_shape_manual(values=c(3, 7, 16, 15, 17))
      
      plot(p)
      
      
    }
  }
}


ABC_lag_plots(request,
              population,
              param_all,
              dist_types,
              quant,
              g_union_check,
              stat_table = TRUE,
              euclidean = FALSE,
              root = TRUE)


ABC_lag_plots(request,
              population,
              param_all,
              dist_types,
              quant,
              g_union_check,
              stat_table = TRUE,
              euclidean = TRUE,
              root = FALSE)
