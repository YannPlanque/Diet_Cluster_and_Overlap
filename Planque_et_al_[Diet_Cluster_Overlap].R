##--------------------------------------------------------------------------------------------------------
## SCRIPT : Identification of harbour and grey seals' diet in the baie de Somme (France), and interspecific dietary overlap
## Specific content : - Identification of diet clusters
##                    - dietary overlap (Pianka index)
##
## As part of : 
##    Planque et al. "Trophic niche overlap between sympatric harbour seals (Phoca vitulina) and grey seals (Halichoerus grypus)
##     at their Southern European limit range (Eastern English Channel)". Submitted in Ecology & Behaviour journal
##
## Author : Yann Planque(1)*
## Affiliation : 
##    (1) Centre d'Etudes Biologiques de Chizé (CEBC, UMR 7372 CNRS - La Rochelle Université), La Rochelle, France## Contact : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr
##
## Contact* : yann.planque@univ-lr.fr ; yann.planque@hotmail.fr##
## Last update : 2020-11-05
##
##
## R version 4.0.2 (2020-06-22) -- "Taking Off Again"
## Copyright (C) 2020 The R Foundation for Statistical Computing
## Platform: x86_64-w64-mingw32/x64 (64-bit)
##-------------------------------------------------.-------------------------------------------------------

### 0 // Useful functions #################################################################################

  ## Bootstrap (with replacement) on diet matrix
  ## Input data: diet matrix (in mass). Row: sample x Column: prey species/type 
  ## Output: diet composition (% in mass by prey)
  ## n = number of rows in the matrix; m = number of bootstrap iterations
    # without CI
    bootproie<-function(n,m,data){
    
    result<-matrix(NA,ncol=ncol(data),nrow=m)
    matrix<-matrix(NA,ncol=ncol(data),nrow=n)
    for (k in 1:m){
      
      est<-seq(1,nrow(data),1)
      estboot<-sample(est,n,replace=T)
      for (j in 1:n){
        matrix[j,]<-as.vector(data[estboot[j],])}
      
      freqboot<-rep(NA,ncol(matrix))
      for (i in 1:ncol(matrix)){
        freqboot[i]<-sum(matrix[,i])}
      freqboot2<-freqboot/sum(freqboot)
      
      result[k,]<-freqboot2}
    result<-as.data.frame(result)
    result
  }
  
    # With CI (at 95% by default)
    bootproie_CI <- function(n,m,data,proba=0.95){
      NAMES <- colnames(data)
      result_prov<-matrix(NA,ncol=ncol(data),nrow=m)
      matrix<-matrix(NA,ncol=ncol(data),nrow=n)
      for (k in 1:m){
        
        est<-seq(1,nrow(data),1)
        estboot<-sample(est,n,replace=T)
        for (j in 1:n){
          matrix[j,]<-as.vector(data[estboot[j],])}
        
        freqboot<-rep(NA,ncol(matrix))
        for (i in 1:ncol(matrix)){
          freqboot[i]<-sum(matrix[,i])}
        freqboot2<-freqboot/sum(freqboot)
        
        result_prov[k,]<-freqboot2}
      
      result_prov
      result_CI <- apply(result_prov,2,quantile, c((1-proba)/2, (1+proba)/2))
      colnames(result_CI) <- 1:(ncol(result_CI))
      result_CI <- as.data.frame(t(result_CI))
      
      result_CI$Prop <- as.vector(apply(data, 2, sum)/sum(data))
      result_CI$Prey <- NAMES
      result_CI <- result_CI[,c(4,3,1,2)]
      colnames(result_CI) <- c("Prey", "Prop", "CI_inf", "CI_sup")
      result <- result_CI
      result
    }
    
    
    ## Bootstrap to calculate number of samples by parameter (default parameter : "Cluster" column in data)
    ## Input data: a data.frame
    ## Output: Proportion of samples by parameter with associated CI
    ## n = number of rows in the matrix; m = number of bootstrap iterations
    boot_N <- function(data, m=1000, proba=0.95, Parameter="Cluster"){
      
      options(dplyr.summarise.inform = FALSE)
      
      test0 <- as.data.frame(
        data %>% group_by_at(Parameter) %>%
          summarise(N = n()))
      test0$Prop_N <- test0$N / sum(test0$N)
      
      Param_raw <- test0[,1]
      
      Prop_param_CI <- NULL
      
      for (k in 1:m){
        
        
        est1<-seq(1,nrow(data),1)
        estboot1<-sample(est1,nrow(data),replace=T)
        
        test1 <- data[estboot1,]
        
        test1b <- as.data.frame(
          test1 %>% group_by_at(Parameter) %>%
            summarise(N = n()))
        
        Param_null <- Param_raw[(Param_raw %in% test1b[,1]) == FALSE] 
        
        if (length(Param_null) == 0){
        }else{
          test1bis <- data.frame(Parameter=Param_null, N=rep(0, length(Param_null)))
          colnames(test1bis) <- c(Parameter, "N")
          test1b <- rbind(test1b, test1bis)
          rm(test1bis)
        }
        rm(Param_null)
        
        test1b$Prop_N <- test1b$N / sum(test1b$N)
        test1b$n_boot <- k
        
        Prop_param_CI <- rbind(Prop_param_CI, test1b)
        
      }
      
      
      test_CI <-
        as.data.frame(
          Prop_param_CI %>% group_by_at(Parameter) %>% 
            summarise(Prop_N_inf = quantile(Prop_N, (1-proba)/2),
                      Prop_N_sup = quantile(Prop_N, 1-(1-proba)/2)))
      
      result <- left_join(test0, test_CI)
      
    }

    
  ## Pianka index (Pianka 1974) / Diet overlap
  ## Input data: two diet matrix (in mass) of the same structure (in columns; 1 column = 1 prey species/type)
  ## Output: Pianka value (0 = no overlap, 1 = complete overlap)
    # Without CI
    pianka_index<-function(data1, data2){
      Test <- colnames(data1) == colnames(data2)
      
      if (sum(Test) == ncol(data1)&sum(Test) == ncol(data2)){
        Tot_1 <- sum(data1)
        Tot_col_1 <- apply(data1,2,sum)
        Data_frame_1 <- data.frame(prey=names(Tot_col_1), mass=as.data.frame(Tot_col_1)[,1])
        Data_frame_1$Prop_tot <- Data_frame_1$mass/Tot_1
        Data_frame_1$Prop_sqr <- (Data_frame_1$Prop_tot)^2
        Sqr_1 <- sum(Data_frame_1$Prop_sqr)
        
        Tot_2 <- sum(data2)
        Tot_col_2 <- apply(data2,2,sum)
        Data_frame_2 <- data.frame(prey=names(Tot_col_2), mass=as.data.frame(Tot_col_2)[,1])
        Data_frame_2$Prop_tot <- Data_frame_2$mass/Tot_2
        Data_frame_2$Prop_sqr <- (Data_frame_2$Prop_tot)^2
        Sqr_2 <- sum(Data_frame_2$Prop_sqr)
        
        result <- sum(Data_frame_1$Prop_tot * Data_frame_2$Prop_tot)/sqrt(Sqr_1*Sqr_2)
        
      }else{
        
        result <- FALSE}
      
    }
    
    # With CI (i.e. with bootstrap on diet data) (at 95% by default)
    pianka_index_CI<-function(n1,n2,m,data1,data2, proba=0.95){
      
      proba_1 <- 1-proba
      CI_inf_prob <- 0+proba_1/2
      CI_sup_prob <- 1-proba_1/2
      
      
      Test <- colnames(data1) == colnames(data2)
      
      if (sum(Test) == ncol(data1)&sum(Test) == ncol(data2)){
        Tot_1 <- sum(data1)
        Tot_col_1 <- apply(data1,2,sum)
        Data_frame_1 <- data.frame(prey=names(Tot_col_1), mass=as.data.frame(Tot_col_1)[,1])
        Data_frame_1$Prop_tot <- Data_frame_1$mass/Tot_1
        Data_frame_1$Prop_sqr <- (Data_frame_1$Prop_tot)^2
        Sqr_1 <- sum(Data_frame_1$Prop_sqr)
        
        Tot_2 <- sum(data2)
        Tot_col_2 <- apply(data2,2,sum)
        Data_frame_2 <- data.frame(prey=names(Tot_col_2), mass=as.data.frame(Tot_col_2)[,1])
        Data_frame_2$Prop_tot <- Data_frame_2$mass/Tot_2
        Data_frame_2$Prop_sqr <- (Data_frame_2$Prop_tot)^2
        Sqr_2 <- sum(Data_frame_2$Prop_sqr)
        
        Pianka_all <- sum(Data_frame_1$Prop_tot * Data_frame_2$Prop_tot)/sqrt(Sqr_1*Sqr_2)
        
        
        colnames(data1)== colnames(data2)
        
        result1<-matrix(NA,ncol=ncol(data1),nrow=m)
        matrix1<-matrix(NA,ncol=ncol(data1),nrow=n1)
        
        for (k in 1:m){
          
          est<-seq(1,nrow(data1),1)
          estboot<-sample(est,n1,replace=T)
          for (j in 1:n1){
            matrix1[j,]<-as.vector(data1[estboot[j],])}
          
          freqboot<-rep(NA,ncol(matrix1))
          for (i in 1:ncol(matrix1)){
            freqboot[i]<-sum(matrix1[,i])}
          freqboot2<-freqboot/sum(freqboot)
          
          result1[k,]<-freqboot2}
        result1<-as.data.frame(result1)
        result1
        
        
        
        result2<-matrix(NA,ncol=ncol(data2),nrow=m)
        matrix2<-matrix(NA,ncol=ncol(data2),nrow=n2)
        
        for (l in 1:m){
          
          est<-seq(1,nrow(data2),1)
          estboot<-sample(est,n2,replace=T)
          for (j in 1:n2){
            matrix2[j,]<-as.vector(data2[estboot[j],])}
          
          freqboot<-rep(NA,ncol(matrix2))
          for (i in 1:ncol(matrix2)){
            freqboot[i]<-sum(matrix2[,i])}
          freqboot2<-freqboot/sum(freqboot)
          
          result2[l,]<-freqboot2}
        result2<-as.data.frame(result2)
        result2
        
        
        Prod1x2 <- apply(result1*result2,1,sum)
        
        result1_sqr <- result1^2
        result2_sqr <- result2^2
        
        sqr1 <- apply(result1_sqr,1,sum)
        sqr2 <- apply(result2_sqr,1,sum)
        
        Pianka <- Prod1x2/sqrt(sqr1*sqr2)
        
        resultbis <- data.frame(Pianka=Pianka_all, CI_inf=quantile(Pianka,CI_inf_prob), CI_sup=quantile(Pianka,CI_sup_prob))
        row.names(resultbis) <- "1"
        result <- resultbis
        
      }else{
        
        result <- FALSE}
      
    }
    
    # To obtain raw pianka values (i.e. with no summarised results)
    pianka_index_CI_RAWDATA<-function(n1,n2,m,data1,data2, proba=0.95){
      
      proba_1 <- 1-proba
      CI_inf_prob <- 0+proba_1/2
      CI_sup_prob <- 1-proba_1/2
      
      
      Test <- colnames(data1) == colnames(data2)
      
      if (sum(Test) == ncol(data1)&sum(Test) == ncol(data2)){
        Tot_1 <- sum(data1)
        Tot_col_1 <- apply(data1,2,sum)
        Data_frame_1 <- data.frame(prey=names(Tot_col_1), mass=as.data.frame(Tot_col_1)[,1])
        Data_frame_1$Prop_tot <- Data_frame_1$mass/Tot_1
        Data_frame_1$Prop_sqr <- (Data_frame_1$Prop_tot)^2
        Sqr_1 <- sum(Data_frame_1$Prop_sqr)
        
        Tot_2 <- sum(data2)
        Tot_col_2 <- apply(data2,2,sum)
        Data_frame_2 <- data.frame(prey=names(Tot_col_2), mass=as.data.frame(Tot_col_2)[,1])
        Data_frame_2$Prop_tot <- Data_frame_2$mass/Tot_2
        Data_frame_2$Prop_sqr <- (Data_frame_2$Prop_tot)^2
        Sqr_2 <- sum(Data_frame_2$Prop_sqr)
        
        Pianka_all <- sum(Data_frame_1$Prop_tot * Data_frame_2$Prop_tot)/sqrt(Sqr_1*Sqr_2)
        
        
        colnames(data1)== colnames(data2)
        
        result1<-matrix(NA,ncol=ncol(data1),nrow=m)
        matrix1<-matrix(NA,ncol=ncol(data1),nrow=n1)
        
        for (k in 1:m){
          
          est<-seq(1,nrow(data1),1)
          estboot<-sample(est,n1,replace=T)
          for (j in 1:n1){
            matrix1[j,]<-as.vector(data1[estboot[j],])}
          
          freqboot<-rep(NA,ncol(matrix1))
          for (i in 1:ncol(matrix1)){
            freqboot[i]<-sum(matrix1[,i])}
          freqboot2<-freqboot/sum(freqboot)
          
          result1[k,]<-freqboot2}
        result1<-as.data.frame(result1)
        result1
        
        
        
        result2<-matrix(NA,ncol=ncol(data2),nrow=m)
        matrix2<-matrix(NA,ncol=ncol(data2),nrow=n2)
        
        for (l in 1:m){
          
          est<-seq(1,nrow(data2),1)
          estboot<-sample(est,n2,replace=T)
          for (j in 1:n2){
            matrix2[j,]<-as.vector(data2[estboot[j],])}
          
          freqboot<-rep(NA,ncol(matrix2))
          for (i in 1:ncol(matrix2)){
            freqboot[i]<-sum(matrix2[,i])}
          freqboot2<-freqboot/sum(freqboot)
          
          result2[l,]<-freqboot2}
        result2<-as.data.frame(result2)
        result2
        
        
        Prod1x2 <- apply(result1*result2,1,sum)
        
        result1_sqr <- result1^2
        result2_sqr <- result2^2
        
        sqr1 <- apply(result1_sqr,1,sum)
        sqr2 <- apply(result2_sqr,1,sum)
        
        Pianka <- Prod1x2/sqrt(sqr1*sqr2)
        
        resultbis <- Pianka
        result <- resultbis
        
      }else{
        
        result <- FALSE}
      
    }
    
  ## Apply a linear translation to make a scaled matrix positive
  Linear_positiv_translation <- function(matrix){
      result <- matrix + abs(min(matrix)) + 0.1
    }
    
##############################################################################################################
  
  
### 0 // Packages ############################################################################################

lapply(c("ggplot2", "dplyr", "lubridate", "NbClust",
         "factoextra", "ggplotify", "ggpubr"), library, character.only=TRUE)

##############################################################################################################
  

### I // Data ################################################################################################
  ## 1 / Direction 
  Direction <- ".../Script_Planque et al_[Diet_Cluster_Overlap]"

  ## DATA VAILABILITY ##
  #  Data used in this script are freely downloadable from SEANOE website: seanoe.org/XXXXXXXXX
  #  Please include data in "Input" file
  ##  
    
  ## 2 / Import data
  #  Data available on SEANOE: 
  #   Planque Yann, Vincent Cécile, Caurant Florence, Spitz Jérôme (2020).
  #   Content of harbour seal (Phoca vitulina) and grey seal (Halichoerus grypus) scats
  #   in prey classified by functional groups (samples collected in the baie de Somme, France,
  #   from 2002 to 2019). SEANOE. https://doi.org/10.17882/76780
  
  #   Diet matrices in reconstructed mass of identified prey (g)
  #   Pv : Phoca vitulina (harbour seals)  //  Hg : Halichoerus grypus (grey seals)
  All_data_raw <- read.table(paste(Direction,"Input","Diet_Seals_BDS_data.csv", sep="/"),
                        dec=".", sep=";", h=T)
  
  Pv_data <- All_data_raw[All_data_raw$Species == "Pv",]
  Hg_data <- All_data_raw[All_data_raw$Species == "Hg",]

  head(Hg_data) # Columns 6 to 11 :  functionnal groups of prey
  
  # Prey names data frame
  Prey_names <- read.table(paste(Direction,"Input","Diet_Seals_BDS_preynames.csv", sep="/"),
                           dec=".", sep=";", h=T)
  Prey_names
  
  
  ## 3 / Prepare data
  Hg_data$Date <- as.Date(Hg_data$Date, "%d/%m/%Y")
  Pv_data$Date <- as.Date(Pv_data$Date, "%d/%m/%Y")
  
  Hg_data$Species <- "Hg"
  Pv_data$Species <- "Pv"
  
  
    # a / Determine seasons according to months
    Hg_data$Season <- ifelse(month(Hg_data$Date) >= 4 & month(Hg_data$Date) <= 8, "Spring_Summer", "Autumn_Winter")
    Hg_data_SS <- Hg_data[Hg_data$Season == "Spring_Summer",]
    Hg_data_AW <- Hg_data[Hg_data$Season == "Autumn_Winter",]
    
    Pv_data$Season <- ifelse(month(Pv_data$Date) >= 4 & month(Pv_data$Date) <= 8, "Spring_Summer", "Autumn_Winter")
    Pv_data_SS <- Pv_data[Pv_data$Season == "Spring_Summer",]
    Pv_data_AW <- Pv_data[Pv_data$Season == "Autumn_Winter",]
    
    All_data <- rbind(Hg_data, Pv_data)

    
    # b / Prepare matrix for diet analyses
    Hg_mtx <- as.matrix(Hg_data[,colnames(Hg_data) %in% Prey_names$Prey])
    row.names(Hg_mtx) <- Hg_data$Num_collection
    
    # Select columns with diet data
    Hg_mtx_SS <- as.matrix(Hg_data_SS[,colnames(Hg_data_SS) %in% Prey_names$Prey]) 
    Hg_mtx_AW <- as.matrix(Hg_data_AW[,colnames(Hg_data_AW) %in% Prey_names$Prey])
    
    Pv_mtx <- as.matrix(Pv_data[,colnames(Pv_data) %in% Prey_names$Prey])
    row.names(Pv_mtx) <- Pv_data$Num_collection
    
    Pv_mtx_SS <- as.matrix(Pv_data_SS[,colnames(Pv_data_SS) %in% Prey_names$Prey])
    Pv_mtx_AW <- as.matrix(Pv_data_AW[,colnames(Pv_data_AW) %in% Prey_names$Prey])

##############################################################################################################

    
### II // Diet composition and interspecific overlap #########################################################
  ## 1 / Diet composition
    # All periods
    Hg_diet <- bootproie_CI(n=nrow(Hg_mtx), m=1000, Hg_mtx, proba = 0.95)
    Hg_diet
    
    Pv_diet <- bootproie_CI(n=nrow(Pv_mtx), m=1000, Pv_mtx, proba = 0.95)
    Pv_diet
    
    # By season
    #Spring/Summer
    Hg_diet_SS <- bootproie_CI(n=nrow(Hg_mtx_SS), m=1000, Hg_mtx_SS, proba = 0.95)
    Pv_diet_SS <- bootproie_CI(n=nrow(Pv_mtx_SS), m=1000, Pv_mtx_SS, proba = 0.95)
    Hg_diet_SS
    Pv_diet_SS
    
    #Autumn/Winter
    Hg_diet_AW <- bootproie_CI(n=nrow(Hg_mtx_AW), m=1000, Hg_mtx_AW, proba = 0.95)
    Pv_diet_AW <- bootproie_CI(n=nrow(Pv_mtx_AW), m=1000, Pv_mtx_AW, proba = 0.95)
    Hg_diet_AW
    Pv_diet_AW
    
    # Plot it...
    Hg_diet <- merge(Hg_diet, Prey_names, by="Prey")
    Pv_diet <- merge(Pv_diet, Prey_names, by="Prey")
    
    Hg_diet$Seal <- "Hg"
    Pv_diet$Seal <- "Pv"
    
    Seal_diet <- rbind(Pv_diet, Hg_diet)
    Seal_diet$Seal <- factor(Seal_diet$Seal, levels=c("Pv", "Hg"))
    
    Plot_raw_diet <- 
    ggplot(data = Seal_diet, 
           mapping = aes(x = reorder(Functionnal_groups, -Num_group_functionnal), fill = Seal, 
                         y = ifelse(test = Seal == "Pv", yes = -Prop*100, no = Prop*100)))+
      geom_bar(stat = "identity", width = 0.5) +
      scale_y_continuous(labels = abs, limits = max(Seal_diet$CI_sup*100) * c(-1,1), breaks=seq(-100, 100, by=20)) +
      labs(y = "Mass (%)", x="") +
      coord_flip()+
      geom_errorbar(aes(ymin=ifelse(test = Seal == "Pv", yes = -CI_inf*100, no = CI_inf*100), 
                        ymax=ifelse(test = Seal == "Pv", yes = -CI_sup*100, no = CI_sup*100)), width=0)+
      scale_fill_manual(labels=c(paste0("Harbour seals\nN=", nrow(Pv_mtx)," scats / 2002-19               "),
                                 paste0("Grey seals\nN=", nrow(Hg_mtx)," scats / 2016-19")),
                        values = c("Pv" = "#4DE600","Hg" = "#0000DE" ), name="")+
      theme_bw(base_size = 16) +
      theme(legend.position="top",
            legend.text=element_text(size=15),
            axis.text=element_text(size=15,color="black") 
      )+
      geom_hline(yintercept=0,color = "black", size=0.5)+
      guides(fill=guide_legend(""))
    
    Plot_raw_diet
    
    ggsave(Plot_raw_diet, filename = paste(Direction,"Plot","00_RAW_Diet_Pv_Hg.png", sep = "/"), dpi = 600,
            width = 9.5, height = 5)
    
    
  ## 2 / Interspecific overlap
    # All periods
    Pianka_Pv_Hg <- pianka_index_CI(n1=nrow(Hg_mtx), data1 = Hg_mtx,
                                    n2=nrow(Pv_mtx), data2 = Pv_mtx,
                                    m=1000)
    Pianka_Pv_Hg  # result
    
    # By season
    Pianka_Pv_Hg_SS <- pianka_index_CI(n1=nrow(Hg_mtx_SS), data1 = Hg_mtx_SS,
                                       n2=nrow(Pv_mtx_SS), data2 = Pv_mtx_SS,
                                       m=1000)
    
    Pianka_Pv_Hg_AW <- pianka_index_CI(n1=nrow(Hg_mtx_AW), data1 = Hg_mtx_AW,
                                       n2=nrow(Pv_mtx_AW), data2 = Pv_mtx_AW,
                                       m=1000)
    
    Pianka_Pv_Hg_SS  # results
    Pianka_Pv_Hg_AW
    
    
    # Probability to have a higher overlap during Spring / Summer 
    # Raw Pianka values
    Pianka_SS_RAW <- pianka_index_CI_RAWDATA(n1=nrow(Hg_mtx_SS), data1 = Hg_mtx_SS,
                                             n2=nrow(Pv_mtx_SS), data2 = Pv_mtx_SS,
                                             m=1000)
    
    Pianka_AW_RAW <- pianka_index_CI_RAWDATA(n1=nrow(Hg_mtx_AW), data1 = Hg_mtx_AW,
                                             n2=nrow(Pv_mtx_AW), data2 = Pv_mtx_AW,
                                             m=1000)
    Hist_pianka_season <- data.frame(Season = c(rep("Spring / Summer", 1000), rep("Autumn / Winter", 1000)),
                                     Pianka = c(Pianka_SS_RAW, Pianka_AW_RAW))
    
      # Histogram of Pianka values (generated 1000 times from bootstrap) by season
      Plot_hist_pianka <- 
      Hist_pianka_season %>%
      ggplot( aes(x=Pianka, fill=Season)) +
        geom_histogram( color="#e9ecef", alpha=0.6, position = 'identity') +
        scale_fill_manual(values=c("#69b3a2", "#404080")) +
        theme_bw(base_size = 16) +
        theme(axis.text=element_text(color="black"))+
        labs(fill="")
      
      Plot_hist_pianka
      
      ggsave(Plot_hist_pianka, filename = paste(Direction,"Plot", "Diet_cluster_00_Hist.png", sep = "/"), dpi = 600,
             width = 7.5, height = 5)
    
    # Probability :
    mean(ifelse(Pianka_SS_RAW > Pianka_AW_RAW, 1, 0))
    
    
    # Plot of Pianka results
    Pianka_results <- 
    rbind(Pianka_Pv_Hg, rbind(Pianka_Pv_Hg_AW, Pianka_Pv_Hg_SS))
    Pianka_results$Period <- c("All periods", "Autumn / Winter", "Spring / Summer")
    
    Plot_pianka <- 
    Pianka_results %>%
    ggplot(aes(x=Period, y=Pianka)) +
      geom_errorbar(aes(ymin=CI_inf, ymax=CI_sup), width=0, size=0.8) +
      geom_point(aes(fill=Period), shape=21, size=5) +
      scale_fill_manual(values=c("#BBBBBB", "#E69F00", "#56B4E9")) +
      ylim(0,1) +
      labs(title = "Diet overlap between harbour and grey seals",
           x="", y="Pianka value [CI95%]") +
      theme_bw(base_size = 16) +
      theme(axis.text.y = element_text(color="black", size=14),
            axis.text.x = element_blank(),
            axis.ticks.x = element_blank(),
            legend.title = element_blank(),
            legend.text = element_text(color="black", size=14), 
            axis.title.y = element_text(color="black", size=16),
            plot.title = element_text(color="black", size=18)) #+
      # geom_hline(yintercept = 0,color = "#999999", size = 0.5) +
      # geom_hline(yintercept = 1,color = "#999999", size = 0.5)
    
    Plot_pianka
    
    ggsave(Plot_pianka, filename = paste(Direction,"Plot", "Diet_cluster_00_Pianka.png", sep = "/"), dpi = 600,
           width = 6, height = 6)
    

# Export results 
write.table(Pianka_results, paste(Direction,"Output","Pianka_results.csv", sep="/"),dec=".", sep=";", row.names = F)

##############################################################################################################


### III // Identification of diet clusters (hierarchical clustering) #########################################
  ## 1 / Prepare a matrix in % mass in prey by sample (i.e. by row)
  Hg_mtx0 <- Hg_mtx
  nrow(Hg_mtx0)

  Pv_mtx0 <- Pv_mtx
  nrow(Pv_mtx0)
  
  
  # Merge both species mass matrices
  All_data0 <- All_data[,colnames(All_data) %in% c("Num_collection", Prey_names$Prey, "Species")]
  All_mtx0 <- as.matrix(All_data0[,colnames(All_data0) %in% Prey_names$Prey])
  
  
  #ID
  row.names(All_mtx0) <- All_data0$Num_collection
  
  # Proportion matrix
  All_mtx <- All_mtx0/rowSums(All_mtx0)
  
  All_mtx
  
  


  ## 2 / Hierarchical clustering on this matrix
  ## Distance: euclidean ; Method = ward.D2
  
    # a / Determine the accurate number of clusters (NbClust package)
    # scale matrix
    Scale_mtx0 <- All_mtx %>% scale() 
    Scale_mtx <- Scale_mtx0[,c(1:ncol(All_mtx))]
    Scale_mtx <- Linear_positiv_translation(Scale_mtx) # need to be positive for NbClust (in R > 4.0)
    
    # number of clusters
    res.nbclust <- 
      NbClust(Scale_mtx, distance = "euclidean",
              min.nc = 2, max.nc = 10, 
              method = "ward.D2", index = "all")
    
    Clust_Nb <- fviz_nbclust(res.nbclust, ggtheme = theme_minimal())
    Clust_Nb

    
    # Best number of clusters : 6
  
    ggsave(Clust_Nb, filename = paste(Direction,"Plot", "Diet_cluster_00_Nb.png", sep = "/"), dpi = 600,
          width = 7.5, height = 5)
    
    K <- as.numeric(as.character(Clust_Nb$data[Clust_Nb$data$freq == max(Clust_Nb$data$freq),][,1])) # Extract nb cluster (K)
    K
    
    # according to plot
    
    # b / Compute clustering (factoextra package)
    # K = 6
    
    # heatmap.2(All_mtx %>% scale(), 
    #           distfun = function(y) dist(y, "euclidean"),
    #           hclustfun=function(x) hclust(x, method="ward.D2") )
    
    res.hc <- All_mtx %>%
      scale() %>%  # Scale data
      eclust("hclust", k = K, graph = FALSE, stand = FALSE,
             hc_method = "ward.D2", hc_metric = "euclidean", nboot = 100) 
    
    Plot_dendro <- fviz_dend(res.hc) # Plot
    Plot_dendro

    # c / Reorder manualy the ID of cluster (to organise it according to plot mapping)
    Plot_dendro_data <- ggplot_build(Plot_dendro)  # ggplot data
    
    Order_col <- data.frame(colour = unique(Plot_dendro_data$data[[2]]["colour"]), Cluster_expected = c(1:K)) # Order of mapping
    
    # Nb ind by colour
    Order_col <- left_join(Order_col, as.data.frame(Plot_dendro_data$data[[2]] %>% group_by(colour) %>% summarise(Nb_ind = n())))
    Order_col # With expected cluster numbers
    
    
    # Order in cluster output
    N_by_clust0 <- data.frame(ID = names(res.hc$cluster), Cluster_output = res.hc$cluster)
    N_by_clust <- as.data.frame(N_by_clust0 %>% group_by(Cluster_output) %>% summarise(Nb_ind=n()))
    N_by_clust <- left_join(N_by_clust, Order_col)  # Merge with expected cluster numbers
    
    N_by_clust$Cluster <- N_by_clust$Cluster_expected
    
    # Distance to map cluster number on dendrogram
    N_by_clust <- N_by_clust %>% arrange(Cluster) %>% mutate(N_mid = Nb_ind/2) %>% mutate(N_cum = cumsum(Nb_ind)) %>% mutate (N_posit = N_cum - N_mid)
    
    
    # d / Final dendrogram
      # Cluster colors
    Col_clusters <- c("firebrick2", "goldenrod", "green3", "deepskyblue", "royalblue4", "darkorchid2") # 6 colors for plot
    
    # Plot A - Dendrogram with k clusters
    Plot_A <- 
      fviz_dend(res.hc, k = res.nbclust, # Cut in four groups
                cex = 0.5, # label size
                color_labels_by_k = TRUE, # color labels by groups
                rect = F, # Add rectangle around groups
                show_labels = F,
                k_colors = Col_clusters, 
                main = "",
                ggtheme = theme_classic()
      )+ theme( axis.line.y = element_line(colour = "black", 
                                           size = 0.65, linetype = "solid"),
                axis.text=element_text(size=13, colour="black"),
                axis.title=element_text(size=14))+
      annotate("text", x = N_by_clust$N_posit, y = -2, label = N_by_clust$Cluster, size=5,
               color=Col_clusters)
    
    Plot_A
    
    ggsave(Plot_A, filename = paste(Direction, "Plot", "Diet_cluster_01_Dendro.png", sep = "/"), dpi = 600,
           width = 6.5, height = 5)
    ###
    
    
    # Clustering validation statistics (silhouette method)
      # values range from 1 to - 1:
          # -> close to 1 indicates that the object is well clustered.
          # -> close to -1 indicates that the object is poorly clustered
    # Inspect the silhouette plot:
    Silh <- fviz_silhouette(res.hc) 
    Silh
    # obs: most of samples are well clustered
    
    # ggsave(Silh, filename = paste(Direction, "Plot", "Diet_cluster_02_Silhouette.png", sep = "/"), dpi = 600,
    #        width = 6.5, height = 5)
    
    
    # e / Assign clusters to data

    Cluster_assign <- data.frame(Num_collection = names(res.hc$cluster), Cluster_output = as.numeric(res.hc$cluster))
    nrow(Cluster_assign) == nrow(All_data)
    
    # New object with raw data (in mass by prey) and cluster result
    All_data_clust <- All_data %>% mutate(Num_collection = as.character(Num_collection))
    All_data_clust <- left_join(All_data, Cluster_assign)
    All_data_clust$Cluster_output == Cluster_assign$Cluster_output
    
    # Assign new cluster id
    All_data_clust <- left_join(All_data_clust, N_by_clust[,c(1, 5)], by="Cluster_output")
    
    # Take seasons and sum_biomass
    All_data_clust$Sum_biomass <- rowSums(All_data_clust[,colnames(All_data_clust) %in% Prey_names$Prey])
    
    
    ## 3 / Distribution of scats in each cluster & associated scat content
      # A / Prepare
      SPP_unique <- unique(All_data_clust$Species)
      Season_unique <- unique(All_data_clust$Season)
      Clusters_unique <- 1:K
      
      All_data_clust$Clust_FFstat <- ifelse(All_data_clust$Cluster %in% c(1,2), 1, 2) # Flatfish clusters, 1 and 2...
      #1 = Flatfish dominance
      #2 = other types of prey
      
      
      # B / All calculations in loops
      Cluster_diet_content <- NULL
      Prop_N_ALL <- NULL
      Prop_N_ALL_season <- NULL
      
      Prop_N_ALL_FFstats <- NULL
      Prop_N_ALL_season_FFstats <- NULL

      for (i in 1:length(SPP_unique)){ # Species level
        
        Data_spp <- All_data_clust[All_data_clust$Species == SPP_unique[i],]
        
        # Proportion of scats in each cluster (by species)
        Prop_N <- boot_N(Data_spp, m=1000, proba=0.95, Parameter="Cluster")
        
        
        # Clusters with zero (to add to data)
        Clust_zero <- unique(All_data_clust$Cluster)[(unique(All_data_clust$Cluster) %in% unique(Prop_N$Cluster)) == FALSE]
        
        if(length(Clust_zero) == 0){
        }else{
          Prop_N_zero <- data.frame(Cluster=Clust_zero, N=rep(0, length(Clust_zero)),
                                    Prop_N = rep(0, length(Clust_zero)), Prop_N_inf = rep(0, length(Clust_zero)), Prop_N_sup = rep(0, length(Clust_zero)))
          Prop_N <- rbind(Prop_N, Prop_N_zero)
        }
        
        Prop_N$Species <- SPP_unique[i]
        
        Prop_N_ALL <- rbind(Prop_N_ALL, Prop_N)
        
        
            ### Same for flatfish clusters statistics...### ///
            Prop_N_FF <- boot_N(Data_spp, m=1000, proba=0.95, Parameter="Clust_FFstat")
            
            # FF_cClusters with zero (to add to data)
            Clust_zero_FF <- unique(All_data_clust$Clust_FFstat)[(unique(All_data_clust$Clust_FFstat) %in% c(1, 2)) == FALSE]
            
            if(length(Clust_zero_FF) == 0){
            }else{
              Prop_N_zero_FF <- data.frame(Cluster=Clust_zero_FF, N=rep(0, length(Clust_zero_FF)),
                                        Prop_N = rep(0, length(Clust_zero_FF)), Prop_N_inf = rep(0, length(Clust_zero_FF)), Prop_N_sup = rep(0, length(Clust_zero_FF)))
              Prop_N_FF <- rbind(Prop_N_FF, Prop_N_zero)
            }
            
            Prop_N_FF$Species <- SPP_unique[i]
            
            Prop_N_ALL_FFstats <- rbind(Prop_N_ALL_FFstats, Prop_N_FF)
            ### ///
        
        rm(Clust_zero, Prop_N_zero, Clust_zero_FF, Prop_N_zero_FF)
        
        for (j in 1:length(Season_unique)){ # Season level
          
          # Proportion of scats in each cluster (by species and seasons)
          Data_spp_season <- Data_spp[Data_spp$Season == Season_unique[j],]
          
          Prop_N_season <- boot_N(Data_spp_season, m=1000, proba=0.95, Parameter="Cluster")
          
          # Clusters with zero (to add to data)
          Clust_zero <- unique(All_data_clust$Cluster)[(unique(All_data_clust$Cluster) %in% unique(Prop_N_season$Cluster)) == FALSE]
          
          if(length(Clust_zero) == 0){
          }else{
            Prop_N_zero <- data.frame(Cluster=Clust_zero, N=rep(0, length(Clust_zero)),
                                      Prop_N = rep(0, length(Clust_zero)), Prop_N_inf = rep(0, length(Clust_zero)), Prop_N_sup = rep(0, length(Clust_zero)))
            Prop_N_season <- rbind(Prop_N_season, Prop_N_zero)
          }
          
          Prop_N_season$Species <- SPP_unique[i]
          Prop_N_season$Season <- Season_unique[j]
          
          Prop_N_ALL_season <- rbind(Prop_N_ALL_season, Prop_N_season)
          
          
              ### Same for flatfish clusters statistics...### ///
              Prop_N_Seaon_FF <- boot_N(Data_spp_season, m=1000, proba=0.95, Parameter="Clust_FFstat")
              
              # FF_cClusters with zero (to add to data)
              Clust_zero_FF <- unique(All_data_clust$Clust_FFstat)[(unique(All_data_clust$Clust_FFstat) %in% c(1, 2)) == FALSE]
              
              if(length(Clust_zero_FF) == 0){
              }else{
                Prop_N_zero_FF <- data.frame(Cluster=Clust_zero_FF, N=rep(0, length(Clust_zero_FF)),
                                             Prop_N = rep(0, length(Clust_zero_FF)), Prop_N_inf = rep(0, length(Clust_zero_FF)), Prop_N_sup = rep(0, length(Clust_zero_FF)))
                Prop_N_Seaon_FF <- rbind(Prop_N_Seaon_FF, Prop_N_zero_FF)
              }
              
              Prop_N_Seaon_FF$Species <- SPP_unique[i]
              Prop_N_Seaon_FF$Season <- Season_unique[j]
              
              Prop_N_ALL_season_FFstats <- rbind(Prop_N_ALL_season_FFstats, Prop_N_Seaon_FF)
              ### ///
          
        }
        
        # Content by cluster
        for (k in 1:K){
          Data_spp_clust <- Data_spp[Data_spp$Cluster == k,]
          
          if(nrow(Data_spp_clust) == 0){
            
            Diet_spp_clust <- data.frame(Prey = Prey_names$Prey[1:6], 
                                         Prop = rep(0, length(Prey_names$Prey[1:6])),
                                         CI_inf = rep(0, length(Prey_names$Prey[1:6])),
                                         CI_sup = rep(0, length(Prey_names$Prey[1:6])),
                                         Species = rep(SPP_unique[i], length(Prey_names$Prey[1:6])),
                                         Cluster = rep(k, length(Prey_names$Prey[1:6])))
              
             
          }else{
          
            Mtx_spp_clust <- as.matrix(Data_spp_clust[,colnames(Data_spp_clust) %in% Prey_names$Prey])
            Diet_spp_clust <- bootproie_CI(n=nrow(Mtx_spp_clust), m=1000, Mtx_spp_clust, proba = 0.95)
          
            Diet_spp_clust$Species <- SPP_unique[i]
            Diet_spp_clust$Cluster <- k
            
          }
          Cluster_diet_content <- rbind(Cluster_diet_content, Diet_spp_clust)
        }
      }
      
      # Prepare results for plotting and exporting
      Prop_N_ALL$Species <- factor(Prop_N_ALL$Species, levels=c("Pv", "Hg"))
      Prop_N_ALL$Cluster <- factor(Prop_N_ALL$Cluster, levels=c(1:K))
      
      Prop_N_ALL_season$Species <- factor(Prop_N_ALL_season$Species, levels=c("Pv", "Hg"))
      Prop_N_ALL_season$Cluster <- factor(Prop_N_ALL_season$Cluster, levels=c(1:K))
      
      Prop_N_ALL_season$Season.labs <- ifelse(Prop_N_ALL_season$Season == "Autumn_Winter", "Autumn / Winter",
                                              ifelse(Prop_N_ALL_season$Season == "Spring_Summer", "Spring / Summer", FALSE))
      
      Cluster_diet_content <- 
        left_join(Cluster_diet_content, Prey_names %>% dplyr::select(Prey, Functionnal_groups, Num_group_functionnal))
      
      Cluster_names <- as.vector(apply(expand.grid("Cluster ", (1:K)), 1, paste, collapse=""))
      names(Cluster_names) <- (1:K)
      
      
      # Numeric results...
      Prop_N_ALL <- Prop_N_ALL %>% dplyr::select(Species, Cluster, N,  Prop_N, Prop_N_inf, Prop_N_sup) 
      Prop_N_ALL_season <- Prop_N_ALL_season %>% dplyr::select(Species, Cluster, Season, Season.labs,  N,  Prop_N, Prop_N_inf, Prop_N_sup)
      Cluster_diet_content <- Cluster_diet_content %>% dplyr::select(Species, Cluster, Prey, Functionnal_groups,
                                                              Num_group_functionnal, Prop, CI_inf, CI_sup)
      
      Prop_N_ALL_FFstats <- Prop_N_ALL_FFstats %>% dplyr::select(Species, Clust_FFstat, N, Prop_N, Prop_N_inf, Prop_N_sup)
      Prop_N_ALL_season_FFstats <- Prop_N_ALL_season_FFstats %>% dplyr::select(Species, Clust_FFstat, Season, N, Prop_N, Prop_N_inf, Prop_N_sup)
      
      Prop_N_ALL
      Prop_N_ALL_season
      Cluster_diet_content
      
      Prop_N_ALL_FFstats
      Prop_N_ALL_season_FFstats
      
      
      #... to be exported
      write.table(Prop_N_ALL, paste(Direction,"Output","Result_prop_N_prey_all.csv", sep="/"),dec=".", sep=";", row.names = F)
      write.table(Prop_N_ALL_season, paste(Direction,"Output","Result_prop_N_prey_seasons.csv", sep="/"),dec=".", sep=";", row.names = F)
      
      write.table(Prop_N_ALL_FFstats, paste(Direction,"Output","Result_prop_N_FLATFISH.csv", sep="/"),dec=".", sep=";", row.names = F)
      write.table(Prop_N_ALL_season_FFstats, paste(Direction,"Output","Result_prop_N_FLATFISH_all.csv", sep="/"),dec=".", sep=";", row.names = F)
      
      write.table(Cluster_diet_content, paste(Direction,"Output","Result_clust_diet_content.csv", sep="/"),dec=".", sep=";", row.names = F)
      
      
      
      # C / Plot that
      
      # Plot B - Proportions of scats in each cluster
      theme_set(theme_bw(base_size = 14))
      MAX_Y <- ceiling(max(Prop_N_ALL$Prop_N_sup)*20)*5
      Plot_B <- 
      ggplot(data = Prop_N_ALL, 
             mapping = aes(x = Cluster, fill = Species, 
                           y = Prop_N*100))+
        geom_bar(stat = "identity", width = 0.6,
                 position=position_dodge()) +
        geom_errorbar(aes(ymin=Prop_N_inf*100, 
                          ymax=Prop_N_sup*100), width=0, stat = "identity", position=position_dodge(0.6))+
        geom_hline(yintercept=0,color = "black", size=0.5)+
        scale_y_continuous(limits=c(0,MAX_Y), breaks=seq(0, MAX_Y, by=20))+
        labs(y = "% scat samples", x="Diet cluster")+
        scale_fill_manual(breaks=c("Pv", "Hg"), labels=c("Harbour seals   ", "Grey seals"), values = c("Pv" = "gray82","Hg" = "Gray28" ))+
        theme(legend.position="none", 
              legend.text=element_text(size=14),
              axis.title = element_text(size=14,color="black"),
              axis.text.y=element_text(size=13,color="black"),
              axis.text.x = element_text(size=13,colour = Col_clusters))
      
      Plot_B
      
      ggsave(Plot_B, filename = paste(Direction,"Plot", "Diet_cluster_03_Distrib_all.png", sep = "/"), dpi = 600,
             width = 5.25, height = 4.75)
      ###
      
      
      
      

      
      # Plot B_season - Proportions of scats in each cluster, by season                                      
      Plot_B_season <- 
      ggplot(data = Prop_N_ALL_season, 
             mapping = aes(x = as.factor(Cluster), fill = Species, 
                           y = Prop_N*100))+
        facet_grid(. ~ Season.labs)+
        geom_bar(stat = "identity", width = 0.6,
                 position=position_dodge()) +
        geom_errorbar(aes(ymin=Prop_N_inf*100, 
                          ymax=Prop_N_sup*100), width=0, stat = "identity", position=position_dodge(0.6))+
        geom_hline(yintercept=0,color = "black", size=0.5)+
        scale_y_continuous(limits=c(0,100), breaks=seq(0, 100, by=20))+
        labs(y = "% scat samples", x="Diet cluster")+
        scale_fill_manual(breaks=c("Pv", "Hg"), labels=c("Harbour seals   ", "Grey seals"), values = c("Pv" = "gray82","Hg" = "Gray28" ))+
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              axis.text.y=element_text(size=13,color="black"),
              axis.text.x = element_text(size=13,colour = "black"), # Col_clusters
              strip.text.x = element_text(
                size = 12, color = "black", face = "bold"))
      
      Plot_B_season
      
      ggsave(Plot_B_season, filename = paste(Direction,"Plot", "Diet_cluster_04_Distrib_season.png", sep = "/"), dpi = 600,
             width = 8.5, height = 4.75)
      ###

      # Plot C - Diet content in scats of each cluster 
      # (prepare the plot)
      Plot_C0 <- 
      Cluster_diet_content %>% 
      ggplot(mapping = aes(x = reorder(Functionnal_groups, -Num_group_functionnal), fill = Species, 
                           y = ifelse(test = Species == "Pv", yes = -Prop*100, no = Prop*100)))+
        facet_wrap(~Cluster, labeller = labeller(Cluster = Cluster_names), ncol=2)+
        geom_bar(stat = "identity", width = 0.6) +
        scale_y_continuous(labels = abs,  limits=c(-100,100), breaks=seq(-100, 100, by=50)) +
        labs(y = "Percentage by mass (%)", x="") +
        coord_flip()+
        geom_errorbar(aes(ymin=ifelse(test = Species == "Pv", yes = -CI_inf*100, no = CI_inf*100), 
                          ymax=ifelse(test = Species == "Pv", yes = -CI_sup*100, no = CI_sup*100)), width=0)+
        scale_fill_manual(breaks=c("Pv", "Hg"), labels=c("Harbour seals   ", "Grey seals"), values = c("Pv" = "gray82","Hg" = "Gray28" ))+
        theme(legend.position="bottom",
              legend.text=element_text(size=14),
              axis.title = element_text(size=14,color="black"),
              axis.text.x =  element_text(size=13,color="black"),
              axis.text.y = element_text(size=12,color="black"),
              strip.background = element_rect(fill=c("green")),
              strip.text.x = element_text(
                size = 12, color = "black", face = "bold"))+
        guides(fill=guide_legend(""))+
        geom_hline(yintercept=0,color = "black", size=0.6)
      
      Plot_C0 # temporary plot
      
      # (plot with colors)
      g <- ggplot_gtable(ggplot_build(Plot_C0))
      stripr <- which(grepl('strip-t', g$layout$name))
      Col_clustersV2 <- Col_clusters[c(5, 6, 3, 4, 1, 2)]
      
      w <- 1
      for (i in stripr) {
        j <- which(grepl('rect', g$grobs[[i]]$grobs[[1]]$childrenOrder))
        g$grobs[[i]]$grobs[[1]]$children[[j]]$gp$fill <- Col_clustersV2[w]
        w <- w+1
      }
      rm(i, j, w)
      
      Plot_C <- ggplotify::as.ggplot(g)
      Plot_C

      ggsave(Plot_C, filename = paste(Direction,"Plot", "Diet_cluster_05_Diet_content.png", sep = "/"), dpi = 600,
             width = 8.5, height = 7.5)
      ###


      # Final plot (A + B + C)
      Plot_diet_clust <-
        ggarrange(
          ggarrange(Plot_A, Plot_B, labels = c("A", "B"), ncol = 2, nrow = 1,
                    font.label = list(size = 15, color = "black", face = "plain")),
          Plot_C, labels = c("", "C"), ncol = 1, nrow = 2,
          font.label = list(size = 15, color = "black", face = "plain"),
          heights = c(1, 2.25)
        )
      Plot_diet_clust

      ggsave(Plot_diet_clust, filename = paste(Direction,"Plot", "Diet_cluster_06_Summary_plot.png", sep = "/"), dpi = 600,
             width = 7.5, height = 9.5)
      
save.image(paste(Direction, "Output", "Diet_analyses.RData", sep = "/"), safe = TRUE)
#load(paste(Direction, "Output", "Diet_analyses.RData", sep = "/"))

##############################################################################################################


####################################################################################################################################
######################################################### END ######################################################################
####################################################################################################################################