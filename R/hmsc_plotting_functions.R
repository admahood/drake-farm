# hmsc_plotting functions

# model convergence ============================================================
ggplot_convergence <- function(Hm, beta = TRUE, V=FALSE, gamma = FALSE, 
                             omega=FALSE, title = "Model Convergence"){
  
  mpost <- convertToCodaObject(Hm)
  
  d <- 
    bind_rows(
  effectiveSize(mpost$Beta) %>%
    as_tibble() %>% mutate(fit_statistic = "ess", variable = "beta"),
  gelman.diag(mpost$Beta, multivariate=FALSE)$psrf%>%
    as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
    mutate(variable = "beta", fit_statistic = "psrf")
    )
  
  if(V) {
    d <- d %>%
      bind_rows(
        effectiveSize(mpost$V) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "V")) %>%
      bind_rows(gelman.diag(mpost$V, multivariate=FALSE)$psrf%>%
          as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
          mutate(variable = "V", fit_statistic = "psrf")
      )
    }
  
  if(gamma) {
    d <- d %>%
      bind_rows(
        effectiveSize(mpost$Gamma) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "gamma")) %>%
      bind_rows(gelman.diag(mpost$Gamma, multivariate=FALSE)$psrf%>%
                  as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
                  mutate(variable = "gamma", fit_statistic = "psrf")
      )
    }
  
  if(omega){
    sppairs = matrix(sample(x = 1:Hm$ns^2, size = 100))
    tmp = mpost$Omega[[1]]
    for (chain in 1:length(tmp)){
      tmp[[chain]] = tmp[[chain]][,sppairs]
    }
    
    d <- d %>%
      bind_rows(
        effectiveSize(tmp) %>%
          as_tibble() %>% mutate(fit_statistic = "ess", variable = "omega")) %>%
      bind_rows(gelman.diag(tmp, multivariate=FALSE)$psrf%>%
                  as_tibble() %>% dplyr::rename(value = `Point est.`) %>%
                  mutate(variable = "omega", fit_statistic = "psrf")
      )
  }
  
  vline_df <- data.frame(fit_statistic = c("ess", "psrf"),
                         xintercept = c(length(mpost$Beta)*nrow(mpost$Beta[[1]]),
                                        1.01))

  
  ggplot(d, aes(x=value)) + 
    geom_histogram(bins=70) + 
    geom_vline(data = vline_df, aes(xintercept = xintercept), color="red", lty=2)+
    facet_grid(variable~fit_statistic, scales='free') +
    ggtitle(title)

}

trace_plot <- function(Hm, which){
  co <- convertToCodaObject(Hm)
  if(which == "beta") return(plot(co$Beta))
  if(which == "gamma") return(plot(co$Gamma))
  if(which == "v") return(plot(co$V))
}
# model fit ====================================================================

ggplot_fit <- function(Hm, which = "r2", sp_names = "none"){
  mpost <- convertToCodaObject(Hm)
  preds <- computePredictedValues(Hm)
  MF <- evaluateModelFit(hM=Hm, predY=preds)
  df <- MF %>% 
    as_tibble() %>%
    pivot_longer(cols = names(.)) 
  spp <- colnames(Hm$Y)
    
  if(which == "named"){
    means <- df%>%
      filter(name == "TjurR2") %>%
      mutate(species = spp)
    if(sp_names != "none") means <- mutate(species = sp_names[species])
    return(ggplot(means, aes(x=value, y=species)) +
      geom_bar(stat = "identity") +
      xlab("Tjur R<sup>2</sup>") +
      theme(axis.title.x = element_markdown()))
  }
  
  if(which == "all"){
    means <- df %>%
      group_by(name) %>%
      summarise(mean = mean(value)) %>%
      ungroup()
    
  return(ggplot(df) +
    geom_histogram(aes(x=value)) +
    facet_wrap(~name) +
    geom_text(data = means, aes(label = paste("Avg =", round(mean, 2))), x=.75, y=4))}
  
  if(which == "r2"){
    means <- df%>%
      filter(name == "TjurR2") %>%
      group_by(name) %>%
      summarise(mean = mean(value),
                max = max(value), 
                min = min(value)) %>%
      ungroup()
    
    return(ggplot(df%>% filter(name == "TjurR2") ) +
      geom_histogram(aes(x=value),bins = 15) +
      ggtitle(paste("Tjur R<sup>2</sup>, Avg:", round(means$mean, 2),
                    ", Range: ", round(means$min, 2)," - ",round(means$max, 2))) +
      theme(plot.title = element_markdown()))
    }
}

# variance partitioning ========================================================

ggplot_vp <- function(Hm, title = "Variance Explained", 
                      cols = NULL, 
                      lut_varnames = NULL, 
                      lut_sppnames = NULL){
  VP <- computeVariancePartitioning(Hm)
  mpost <- convertToCodaObject(Hm)
  
  prevalence <-   colSums(Hm$Y) %>%
    as_tibble(rownames = "Species") %>%
    dplyr::rename(prevalence = value) %>%
    arrange(desc(prevalence))

  
  mf_df <- data.frame(Species = colnames(Hm$Y)) %>%
    left_join(prevalence)
  
  vp_df <- VP$vals%>%
    as_tibble(rownames = "variable") %>%
    pivot_longer(cols=names(.)[2:ncol(.)], 
                 names_to = "Species", 
                 values_to = "value") %>%
    left_join(prevalence) %>%
    na.omit()
  
  # vp_summary <- vp_df %>%
  #   group_by(origin, variable) %>%
  #   summarise(value = mean(value)) %>%
  #   ungroup() %>%
  #   pivot_wider(id_cols = "variable", 
  #               names_from = "origin", 
  #               values_from = "value")
  
  if(is.vector(lut_varnames)) vp_df <- vp_df %>% mutate(variable = lut_varnames[variable])
  if(is.vector(lut_sppnames)) vp_df <- vp_df %>% mutate(Speices = lut_varnames[Species])

  vp_summary <- vp_df %>%
    group_by(variable) %>%
    summarise(value_pct = mean(value) * 100) %>%
    ungroup() %>%
    mutate(variable_pct = paste0(variable, " (", round(value_pct,1), "%)")) %>%
    dplyr::select(-value_pct)
  
  vp_order <- vp_df %>%
    filter(variable == first(vp_df$variable %>% unique())) %>%
    arrange(prevalence) %>%
    mutate(Species_f = factor(Species, levels = .$Species)) %>%
    dplyr::select(Species, Species_f) 
  
  
  p_vp <- left_join(vp_df, vp_order) %>% 
    left_join(vp_summary) %>%
    mutate(variable = factor(variable),
           value = value) %>%
    ggplot(aes(x=value,y=Species_f, fill = variable_pct)) +
    geom_bar(stat="identity", color = "black")+
    theme_classic() +
    scale_fill_discrete(name = "Variable\n (Avg Variance Explained)") +
    geom_hline(yintercept = table(vp_order$origin)[1]+.5) +
    geom_hline(yintercept = nrow(vp_order)+.5) +
    annotate("text", x = 1.2, y=1, label="Introduced", angle=90, vjust="bottom",
             hjust="left", size=8)+
    annotate("text", x = 1.2, y=nrow(vp_order), label="Native", angle=90, vjust="top",
             hjust="right", size=8)+
    ylab("Species") +
    xlab("Proportion of Variance Explained") +
    theme(legend.position = "right",
          legend.text = element_markdown(),
          legend.title = element_blank(),
          legend.justification = c(1,0),
          legend.background = element_rect(color="black")) 
  
  if(is.vector(cols)) p_vp <- p_vp + scale_fill_manual(values = cols)
  if(is.vector(title)) p_vp <- p_vp + ggtitle(title)
     
  return(p_vp)

  }

# betas ========================================================================


# add species exlusion?, add fill behind require(ggnewscale)
ggplot_beta <- function(Hm,
                        grouping_var = NA,
                        grouping_labels = c("Introduced", "Native"),
                        support_level = 0.89,
                        lut_varnames = NULL, 
                        lut_sppnames = NULL, no_intercept = TRUE, title = NA){
  postBeta <- getPostEstimate(Hm, parName = "Beta")
  
  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }
  

  
  means <- postBeta$mean %>%
    as_tibble() %>%
    rowid_to_column("env_var") %>%
    mutate(env_var = c(covNames)) %>%
    pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Species", values_to = "Mean")
  

  supported <- postBeta$support %>% 
    as_tibble() %>%
    rowid_to_column("env_var") %>%
    mutate(env_var = covNames) %>%
    pivot_longer(cols=names(.)[2:ncol(.)], 
                 names_to = "Species", 
                 values_to = "Support") %>%
    filter(Support > support_level | Support < (1-support_level),
           env_var != "(Intercept)") %>%
    left_join(means, by = c("env_var", "Species"))%>%
    mutate(sign = ifelse(Mean>0, "+", "-"))
  
  vp_order <-   colSums(Hm$Y) %>%
    as_tibble(rownames = "Species") %>%
    dplyr::rename(prevalence = value) %>%
    left_join(Hm$TrData %>% tibble::rownames_to_column("Species")) %>%
    arrange(desc(introduced), (prevalence)) %>%
    mutate(Species_f = factor(Species, levels = .$Species)) %>%
    filter(Species %in% supported$Species)
  
  supported <- supported %>%
    left_join(vp_order)#
  
  if(!is.na(grouping_var)){
    tab <- vp_order %>%
      pull(grouping_var) %>%
      table()
    
    iline_position <- tab[2] + 0.5  
    txt <- data.frame(label = grouping_labels,
                      x = c(1.02,1.02), 
                      y = c(iline_position - 5, iline_position + 5))
    y_position <- length(unique(supported$env_var)) + .8
  }
  
  # length(vp_order$origin[vp_order$origin=="I"]) + .5
  
  if(is.vector(lut_varnames)) supported <- supported %>% mutate(env_var = lut_varnames[env_var])
  if(is.vector(lut_sppnames)) supported <- supported %>% mutate(Speices = lut_varnames[Species])
  if(no_intercept) supported <- supported %>% filter(env_var != "(Intercept)")
  
  p_beta <- supported %>%
    ggplot(aes(x=env_var,y=reorder(Species_f,Species))) +
    geom_tile(lwd=.5, aes(fill = Mean, color = sign)) +
    theme_pubclean()+
    scale_fill_steps2() +
    scale_color_manual(values = c(("red"), ("blue"))) +
    guides(color = "none")+
    scale_x_discrete(expand = c(0,1)) +
    theme(axis.text.x = element_markdown(angle=45, vjust=1,hjust = 1),
          legend.position = "right",
          plot.background = element_rect(color="black"),
          plot.title = element_text(hjust = 1, face = "bold")) +
    xlab("Environmental Filters")+
    ylab("Species")
  
  if(!is.na(grouping_var)) {
    p_beta <- p_beta +
    geom_hline(data = txt, yintercept = iline_position) +
    geom_text(data = txt, aes(label = label, x=c(y_position,y_position),y=y),
              angle=90, fontface="bold")
  }

  if(!is.na(title)) p_beta <- p_beta + ggtitle(title)
  
  return(p_beta)
}

# traits =======================================================================

ggplot_gamma <- function(Hm, support_level = 0.89, no_intercept = TRUE, title = "Effects on Traits"){
  
  covNamesNumbers <- c(TRUE, FALSE)
  covNames = character(Hm$nc)
  for (i in 1:Hm$nc) {
    sep = ""
    if (covNamesNumbers[1]) {
      covNames[i] = paste(covNames[i], Hm$covNames[i], sep = sep)
      sep = " "
    }
    if (covNamesNumbers[2]) {
      covNames[i] = paste(covNames[i], sprintf("(C%d)", i), sep = sep)
    }
  }
  
  trNames = character(Hm$nt)
  trNamesNumbers = c(T,F)
  for (i in 1:Hm$nt) {
    sep = ""
    if (trNamesNumbers[1]) {
      trNames[i] = paste(trNames[i], Hm$trNames[i], sep = sep)
      sep = " "
    }
    if (trNamesNumbers[2]) {
      trNames[i] = paste(trNames[i], sprintf("(T%d)", i), 
                         sep = sep)
    }
  }
  
  trNames <- str_remove_all(trNames, "yes") %>%
    str_remove("pp")
  
  postGamma = getPostEstimate(Hm, parName="Gamma")

  means_gamma <- postGamma$mean %>%
    as_tibble() %>%
    rowid_to_column("env_var") %>%
    mutate(env_var = c(covNames)) %>%
    pivot_longer(cols=names(.)[2:ncol(.)], names_to = "Trait", values_to = "Mean")
  
  lut_gamma <- trNames
  names(lut_gamma) <- unique(means_gamma$Trait)
  
  supported_gamma <- postGamma$support %>% 
    as_tibble() %>%
    rowid_to_column("env_var") %>%
    mutate(env_var = covNames) %>%
    pivot_longer(cols=names(.)[2:ncol(.)], 
                 names_to = "Trait", 
                 values_to = "Support") %>%
    filter(Support > support_level |Support< (1-support_level),
           env_var != "(Iintercept)") %>%
    left_join(means_gamma, by = c("env_var", "Trait"))%>%
    mutate(sign = ifelse(Mean>0, "+", "-"),
           Trait = lut_gamma[Trait])%>%
    filter(env_var != "(Intercept)")
  
  p_gamma <- supported_gamma %>%
    ggplot(aes(x=env_var,y=(Trait), fill = Mean, color = sign)) +
    geom_tile(lwd=.5) +
    theme_pubclean()+
    scale_fill_steps2() +
    scale_color_manual(values = c(("red"), ("blue"))) +
    guides(color = "none")+
    theme(axis.text.x = element_text(angle=45, vjust=1,hjust = 1),
          # axis.title = element_blank(),
          legend.position = "right",
          plot.background = element_rect(color="black"),
          plot.title = element_text(hjust = 1, face = "bold")) +
    ggtitle(title) +
    xlab("Environmental Filters") +
    ylab("Traits")
  
  return(p_gamma)
}

# omegas =======================================================================
omega_table <- function(Hm) {
  # make a thing to include the mean and abs value for species 2 origin
  OmegaCor = computeAssociations(Hm)
  
  hmdf_mean <- OmegaCor[[1]]$mean %>%
    as.matrix
  hmdf_support <- OmegaCor[[1]]$support %>%
    as.matrix
  
  # avg association strengths
  tab<- OmegaCor[[1]]$mean %>%
    abs() %>%
    rowSums() %>%
    as_tibble(rownames = "Species") %>%
    arrange(desc(value)) %>%
    left_join(prevalence) %>%
    left_join(OmegaCor[[1]]$mean %>%
                rowSums() %>%
                as_tibble(rownames = "Species") %>%
                dplyr::rename(mean = value)) %>%
    dplyr::rename(abs_value = value)
  return(tab)
}


ggplot_omega <- function(Hm, 
                         support_level = 0.89, 
                         hc_method = "single", 
                         axis_text_colors = "black",
                         title = "Species Associations",
                         dots = TRUE){
  OmegaCor = computeAssociations(Hm)

  hmdf_mean <- OmegaCor[[1]]$mean %>%
    as.matrix
  hmdf_support <- OmegaCor[[1]]$support %>%
    as.matrix

  pmat <-OmegaCor[[1]]$support 
  pmat[pmat < (1- support_level)] <- 0.99
  
  get_upper_tri <- function(cormat){
    cormat[upper.tri(cormat)]<- NA
    return(cormat)
  }
  reorder_cormat <- function(cormat){
    # Use correlation between variables as distance
    dd <- as.dist((1-cormat)/2)
    hc <- hclust(dd, "ave")
    cormat <-cormat[hc$order, hc$order]
    return(cormat)
  }
  
  # if(bar){
  # 
  # hmdf_mean %>%
  #   # reorder_cormat() %>%
  #   get_upper_tri() %>%
  #   as_tibble(rownames = "species1") %>%
  #   pivot_longer(names_to = "species2", values_to = "rescor", -species1) %>%
  #   mutate(rescor = ifelse(rescor ==1, NA, rescor)) %>%
  #   left_join(traits %>% as_tibble(rownames = "species1")) %>%
  #   left_join(traits %>% as_tibble(rownames = "species2") %>% dplyr::select(species2, origin2=origin)) %>%
  #   mutate(origin = ifelse(origin == "I", "Introduced", "Native")) %>%
  #     mutate(origin2 = ifelse(origin2 == "I", "Introduced", "Native")) %>%
  #   filter(!is.na(rescor)) %>%
  #   ggplot(aes(x=species1, y = species2, fill = rescor)) +
  #   # geom_tile(aes(color = origin2),key_glyph = "path", lwd=.25) +
  #   # scale_color_manual(values = c("red", "black")) +
  #     geom_tile(color="black") +
  #   scale_fill_gradient2(mid = "grey90") +
  #   facet_grid(origin2~origin, scales = "free_x")+
  #   theme(axis.text.x = element_text(angle =90),
  #         panel.grid.major.y = element_line(color = "grey", linetype = 2))
  #  
  #   
  #   
  #   hmdf_mean %>%
  #     as_tibble(rownames = "species1") %>%
  #     pivot_longer(names_to = "species2", values_to = "rescor", -species1) %>%
  #     left_join(traits %>% as_tibble(rownames = "species1")) %>%
  #     left_join(traits %>% as_tibble(rownames = "species2") %>% dplyr::select(species2, origin2=origin)) %>%
  #     mutate(origin = ifelse(origin == "I", "Introduced", "Native")) %>%
  #     mutate(origin2 = ifelse(origin2 == "I", "Introduced", "Native")) %>%
  #     filter(!is.na(rescor)) %>%
  #     ggplot(aes(x=abs(rescor), color = paste0(origin, ":", origin2),
  #                fill = paste0(origin, ":", origin2))) +
  #     geom_density(alpha=0.25)
  #     # group_by(origin, origin2) %>%
  #     # summarise(meanabs = mean(abs(rescor)),
  #     #           mean = mean(rescor)) %>%
  #     # ungroup() %>%
  #     # arrange(desc((meanabs)))
  # }
  if(dots){
    pcor1<- ggcorrplot::ggcorrplot(hmdf_mean,
                                 type = "lower",
                                 hc.order = TRUE,
                                 hc.method = hc_method,
                                 colors = c("red", "white", "blue"),
                                 p.mat = pmat,
                                 pch = 20,
                                 sig.level = support_level,
                                 title = title) +
    theme(plot.background = element_rect(color="black"),
          legend.position = c(0,1),
          legend.justification = c(0,1),
          legend.background = element_rect(color="black"),
          plot.title = element_text(hjust = 1, face = "bold"))
    }else{
      pcor1<- ggcorrplot::ggcorrplot(hmdf_mean,
                                     type = "lower",
                                     hc.order = TRUE,
                                     hc.method = hc_method,
                                     colors = c("red", "grey90", "blue"),
                                     title = title, 
                                     legend.title = "Residual\nCorrelation") +
        theme(plot.background = element_rect(color="black"),
              legend.position = c(0,1),
              axis.text = element_text(color = axis_text_colors),
              legend.justification = c(0,1),
              legend.background = element_rect(color="black"),
              plot.title = element_text(hjust = 1, face = "bold"))
          }

  return(pcor1)
}
