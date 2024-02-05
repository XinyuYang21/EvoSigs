

#' Build CCF Matrix for samples
#' @name ccfMatBuild
#' @param samplelist samplename list, eg., TCGA-05-1425
#' @param ccf_upper ccf upper bound
#' @param ccf_folder ccf files folder path
#' @param binwidth bin width
#' @param random output randomized ccf matrix
#' @return ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat and ccfBandFractionMat.random
#' @export
#' @importFrom NMF randomize
ccfMatBuild <- function(samplenamelist,ccf_upper=1,ccf_folder,binwidth=0.01,random=FALSE){
  
  n_sample <- length(samplelist)
  rows <- ccf_upper/binwidth + 1
  ccfBand <- seq(0,upper,length.out = rows)
  ccfBandCountsMat <- matrix(nrow = rows,ncol = n_sample)
  
  if (n_sample == 0) stop("The number of samples is 0")
  
  for (i in 1:n_sample){
    
    sample_name <- as.character(samplelist[i])
    ssm <- load_ccf(sample_name,input = ccf_folder)
    matchBandCountDf <- data.frame(Var1=as.character(1:rows))
    matchBandCountDf <- suppressWarnings(left_join(matchBandCountDf,as.data.frame(table(findInterval(ssm$ccube_ccf,ccfBand))),stringAsFactors = F,by="Var1"))
    ccfBandCountsMat[,i] <- matchBandCountDf$Freq
    
  }
  
  if (random==TRUE){
    
    if (n_sample == 1) {
      ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
      ccfBandFractionMat.random <- ccfBandFractionMat[sample(1:(rows-1))]
      ccfBandCountsMat.random <- ccfBandCountsMat[sample(1:(rows-1))]
    } else {
      ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
      ccfBandCountsMat.random <- randomize(ccfBandCountsMat)
      ccfBandFractionMat.random <- apply(ccfBandCountsMat.random,2,function(x) x/sum(x))
    }
    outputlist <- list(ccfBandCountsMat,ccfBandCountsMat.random,ccfBandFractionMat,ccfBandFractionMat.random)
    
  } else {
    
    if (n_sample == 1) {
      ccfBandFractionMat <- ccfBandCountsMat/sum(ccfBandCountsMat)
    } else {
      ccfBandFractionMat <- apply(ccfBandCountsMat,2,function(x) x/sum(x))
    }
    outputlist <- list(ccfBandCountsMat,ccfBandFractionMat)
  }
  
  outputlist <- lapply(outputlist,function(x){x <- as.data.frame(t(x)) %>% mutate(samplename=samplelist)})
  
  return(outputlist)
}



#' Build count matrix for input samples per cancer type
#' @name ccfMatrixBuild
#' @param ccf_upper CCF upper bound
#' @param ccf_folder ccf files folder
#' @param output output folder path
#' @param RankEstimateType output rank estimate matrix format
#' @param minsample minimum samples required for each cancer type
#' @return CCF matrix for each cancer type
#' @export
#' @import dplyr
ccfMatBuild_output <- function(ccf_folder,output=NA,ccf_upper=1,RankEstimateType="fraction",minsample=30){
  
  post_summary = Build_post_summary_TCGA(ccf_folder=ccf_folder,minsample=minsample)
  samplelist_all <- post_summary$samplename
  type <- post_summary %>% group_by(.data$cancertype)%>% summarize(n=n())
  typelist <- as.character(subset(type,n>=minsample)$cancertype)
  ntype <- length(typelist)
  
  cat(paste0("\n Start constructing ccf count matrix for ",ntype," types (>",minsample," samples) \n"))
  
  if (!is.na(output)) {
    ccfOutput_path <- paste0(output,"ccfMat/")
    ccfOutput_all_path <- paste0(output,"ccfMat/All/")
    ccfOutput_rank_path <- paste0(output,"ccfMat/rank_estimate/")
    
    multi_dir_create(c(ccfOutput_all_path,ccfOutput_rank_path))
    
    for (i in 1:ntype){
      
      type <-  typelist[i]
      cat(paste0("\n >>>> loading ",i,'th type - ',type," <<<< \n"))
      
      samplelist_type <- subset(post_summary,cancertype==type)$samplename
      n_sample <- length(samplelist_type)
      
      ccfBandCountsMat <- suppressWarnings(ccfMatBuild(samplelist_type,ccf_folder = ccf_folder,upper=ccfupper))
      
      ccfCountMatrix <-ccfBandCountsMat[[1]]
      ccfCountsMatrix.random <- ccfBandCountsMat[[2]]
      ccfFractionMatrix <- ccfBandCountsMat[[3]]
      ccfFractionMatrix.random <- ccfBandCountsMat[[4]]
      
      # create directio
      multi_dir_create(paste0(ccfOutput_path,type,"/"))  
      
      # Output Matrix 
      output_format <- paste0(ccfOutput_path,type,"/",type,"_", n_sample,"_0-",ccfupper)
      output_format_rank <- paste0(ccfOutput_rank_path,type,"_", n_sample,"_0-",ccfupper)
      
      save(ccfCountMatrix,file=paste0(output_format,"_ccfCountMatrix_",Sys.Date(),".RData"))
      save(ccfCountsMatrix.random,file=paste0(output_format,"_ccfCountMatrix.random_",Sys.Date(),".RData"))
      save(ccfFractionMatrix,file=paste0(output_format,"_ccfFractionMatrix_",Sys.Date(),".RData"))
      save(ccfFractionMatrix.random,file=paste0(output_format,"_ccfFractionMatrix.random_",Sys.Date(),".RData"))
      
      # output to rank estimate folder
      if (RankEstimateType=="fraction") {
        write.csv(ccfFractionMatrix.random,file=paste0(output_format_rank,"_ccfFractionMatrix.random_",Sys.Date(),".csv"))
        write.csv(ccfFractionMatrix,file=paste0(output_format_rank,"_ccfFractionMatrix_",Sys.Date(),".csv"))
      }
      
      if (RankEstimateType=="count") {
        write.csv(ccfCountsMatrix.random,file=paste0(output_format_rank,"_ccfCountMatrix.random_",Sys.Date(),".csv"))
        write.csv(ccfCountMatrix,file=paste0(output_format_rank,"_ccfCountMatrix_",Sys.Date(),".csv"))
      }
    }
  }
  
  # ouput ccf matrix for all samples
  ccfBandCountsMat_all <- suppressWarnings(ccfMatBuild(samplelist_all,input_folder = input_folder,upper=ccfupper,add_samplename = add_samplename))
  
  ccfCountMatrix_all <- ccfBandCountsMat_all[[1]]
  ccfFractionMatrix_all <- ccfBandCountsMat_all[[3]]
  
  output_all_format <- paste0(ccfOutput_all_path,"All_",length(samplelist_all),"_0-",ccfupper)
  
  # Output Matrix 
  save(ccfCountMatrix_all,file=paste0(output_all_format,"_ccfCountMatrix_",Sys.Date(),".RData"))
  save( ccfFractionMatrix_all,file=paste0(output_all_format,"_ccfFractionMatrix_",Sys.Date(),".RData"))
  
  cat(paste0("\n Done, saved results to ",output))
}


#' Combine signatures from different cancer type together
#' @name combine_sig_nmf
#' @param sig_folder signatures folder
#' @param cancertype cancer type list
#' @param output_folder output_folder
#' @export
#' @return matrix combining all signature matrix
combine_sig_nmf <- function(typewise_sigs,output_folder=NA,cancertype){
  
  for (i in 1:length(cancertype)){
    tryCatch({
      type <- cancertype[i]
      cat(paste0("-> loading Evolutionary signatures for ",i,"th type : ",type),"\n")
      sig_file <- dir(paste0(input_folder,type,"/"))[grep("sig",dir(paste0(input_folder,type,"/")))]
      load(paste0(input_folder,type,"/",sig_file))
      colnames(sig) <- paste0(type,"_sig",1:ncol(sig))
      
      if (i==1) {
        combine_sigs <- sig
      } else {
        combine_sigs <- cbind(combine_sig,sig)
      }
      
    },error=function(e) print("Fail load this cancer type"))
  }  
  
  if (!is.na(output_folder)) save(combine_sigs,file=paste0(output_folder,"combine_sig_",Sys.Date(),".RData"))
  
  return(combine_sigs)
}


#' Perform NMF with consensus signature
#' @param ccfMat ccf matrix for all samples
#' @param consensus_sig consensus signatures(columns as signaguture,row as bins)
#' @param output output folder path
#' @return exposure
#' @export
Extract_sig <- function(ccfMat,consensus_sig,output=NA){
  
  n_sig <- ncol(consensus_sig)
  Mat <- t(ccfMat[,1:100])
  Mat[is.na(Mat)] <- 0
  
  EvoSig_exposure <- as.data.frame(t(YAPSA::LCD(Mat,consensus_sig))) %>%
  EvoSig_exposure <- cbind(EvoSig_exposure,apply(EvoSig_exposure[,1:n_sig],2,function(x) x/sum(x))) 
  colnames(EvoSig_exposure) <- c(paste0("Evo_sig_",1:n_sig),paste0("Evo_sig_",1:n_sig,"_proportion"))
  EvoSig_exposure$samplename = ccfMat[,101]
  EvoSig_exposure$samplename = substr(filename$samplename,1,12)
  EvoSig_exposure$samplename = gsub("[.]","-", EvoSig_exposure$samplename)
    
  if (!is.na(output)) {
    if (!dir.exists(output)) {dir.create(output)} 
    save(EvoDynamics_exposure,file=paste0(output,"evoDynamic_exposure.RData"))
  }
  
  return(EvoSig_exposure)
}

#' Choose the Number of Consensus Signatures using Hierarchical Clustering Number Estimate
#'
#' This function utilizes hierarchical clustering and various clustering indices to estimate
#' the optimal number of consensus signatures for a given set of type-wise signatures.
#'
#' @title chooseNumOfConsensusSigs
#' @param combine_sig Combined type-wise signatures across cancer types
#' @param min Minimum clustering number
#' @param max Maximum clustering number
#' @return A data frame with suggested clustering numbers for different methods and distances
#' @import NbClust
#' @export
chooseNumOfConsensusSigs <- function(combine_sig,min = 2,max = 10){
  
  # Define available clustering methods and distances
  methods = c("kl","ch","hartigan",
                    "cindex","db","silhouette","ratkowsky","ball",
                    "ptbiserial","gap", "frey", "mcclain",  "gamma", "gplus", "tau", "dunn", 
                    "sdindex", "sdbw") # "hubert","dindex"
  distance = c("metodo","euclidean", "maximum", "manhattan", "canberra")
  
  # Function to get the mode of a vector
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  # Initialize result table
  result_table <- as.data.frame(matrix(ncol = length(distances), nrow = length(methods)))
  colnames(result_table) <- distances
  
  # Loop through distances and methods to estimate clustering numbers
  for (j in 2:length(distance))
    for(i in 1:length(methods)){
      tryCatch({
        nb <- NbClust::NbClust(combine_sigs,distance = distance[j],
                     min.nc = min, max.nc = max, 
                     method = "complete", index =methods[i])
        result_table[i, j] <- nb$Best.nc[1]
        result_table[i, 1] <- methods[i]
      },error=function(e) print("error"))
    } 
  
  # Add a row indicating the most common clustering number for each distance
  result_table <- rbind(result_table,c("Most_common",apply(result_table[,2:5],2,getmode)))
  
  return(result_table)
}


#' Combine signatures and estimate number of clusters
#' @name hc_combine_sig Combine signatures and estimate number of clusters
#' @param combined_sigs combined signatures for all cancer types
#' @param distance distance function
#' @param output_folder output folder path
#' @param min.nc minimum cluster
#' @param max.nc maximum cluster
#' @return list(combine_sig,cluster)
#' @export
hc_combine_sig <- function(combine_sigs,output_folder,min.nc=3,max.nc=10,method="ward.D2",index="hubert",distance="euclidean"){
  cat(red("Start run hierarchical clustering to extract consensus signature. \n"))
  
  #combine_sig <- combine_sig_nmf(input_folder=input_folder,output_folder=output_folder,cancertype=cancertype) 
  
  if (index=="hubert") {
    
    pdf(paste0(output_folder,"hc_hubert_plot_",Sys.Date(),".pdf")) 
    
    res.nbclust  <- NbClust::NbClust(t(combine_sigs),min.nc = min.nc, max.nc = max.nc,method=method,index = index,distance=distance)
    second_diff = res.nbclust$All.index[2:(max.nc-min.nc+1)]-res.nbclust$All.index[1:(max.nc-min.nc)]
    cluster = as.numeric(which(second_diff==max(second_diff)))+min.nc
    print(paste0("The suggested number of clusters from the Hubert index was - ",cluster))
    cat(paste0("-> The Hubert index was saved as: ",paste0(output_folder,"hc_hubert_plot.pdf \n")))
    
    dev.off()
  }
  factoextra::fviz_nbclust(res.nbclust)
  return(list(combine_sigs,cluster))
}

#' Plotting for HC results 
#' @name hc_consensus ccfMat ccf matrix for all samples
#' @param combine_sig combined signatures for all cancer types
#' @param cluster clustering number
#' @param output_folder output folder path
#' @param distance distance function
#' @return consensus_sig
#' @import pheatmap
#' @importFrom cowplot save_plot
#' @import dplyr
#' @export
hc_consensus_sig <- function(combine_sig,cluster,output_folder,method="ward.D2",index="hubert",distance="euclidean"){
  
  combine_sig <- apply(combine_sig,2,function(x) x/sum(x))
  upper = quantile(as.matrix(combine_sig),0.95);breaksList = seq(0, upper, by = 0.01)
  col <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(n = cluster, name = "RdYlBu")))(length(breaksList))
  
  out <- pheatmap::pheatmap(combine_sig, cutree_cols = cluster, fontsize_col = 5,fontsize_row = 0.4,color = col, breaks = breaksList,clustering_distance_cols=distance,cluster_rows=F,filename=paste0(output_folder,distance,"_",cluster,"_hc_heatmap_",Sys.Date(),".pdf"),clustering_method = method)
  
  sig_label <- as.data.frame(cutree(out$tree_col,k=cluster)) 
  colnames(sig_label) = "cluster"
  sig_label$sig = rownames(sig_label)
  
  ## Compute consensu signatures for each cluster
  consensus_sig <- as.data.frame(t(combine_sig)) %>%
    mutate(sig=rownames(.)) %>%
    left_join(sig_label,by="sig") %>%
    group_by(cluster) %>%
    mutate(sig=NULL) %>%
    summarise_all(mean) 
  
  save(consensus_sig,file=paste0(output_folder,distance,"_",cluster,"_consensus_sig_",Sys.Date(),".RData"))
  
  consensus_sig <- apply(t(consensus_sig[,2:101]),2,as.numeric) 
  p1 <- plot_grid(sig_plot(consensus_sig))
  save_plot(paste0(output_folder,distance,"_",cluster,"_consensus_sig_",Sys.Date(),".pdf"),p1,base_asp = cluster)
  sig_plot(consensus_sig)
  cat(paste0("-> The consensus signature result was saved at: ",output_folder,"\n"))
  
  return(consensus_sig) 
}

#' sig_assignment 
#' @name hc_consensus ccfMat ccf matrix for all samples
#' @param signature signature
#' @param ccfMatrix ccfmatrix
#' @param output output folder path
#' @return exposure
#' @import dplyr
#' @importFrom cowplot save_plot
#' @importFrom YAPSA LCD
#' @export
sig_assignment <- function(signature,ccfMatrix,output=NA){
  
  ccfMatrix[is.na(ccfMatrix)] = 0
  
  exposure <-LCD(ccfMatrix,signature) 
  exposure <- as.data.frame(t(exposure)) %>% set_colnames(paste0("sig_",1:ncol(.))) 
  
  if (!is.na(output)) save(exposure,file=paste0(output,"/lcd_exposure",Sys.Date(),".RData"))
  
  return(exposure)
}



#' Plot signature matrix
#' @name sig_plot
#' @param sig signature matrix
#' @return mydata mutation data frame
#' @export
#' @importFrom reshape2 melt
#' @importFrom grid grid.draw
#' @importFrom RColorBrewer brewer.pal
#' @importFrom magrittr %>% set_colnames
#' @import dplyr
ccf_dist <- function(ccf_rows,
                     density_plot=FALSE,
                     tag=NULL,col=NA,theme="grey",fraction=FALSE,strip_title="Signature ",
                     title="Consensus Signature of evolutionary dynamics"){
  
  if (fraction) ccf_rows <- t(apply(ccf_rows,1,function(x) x/sum(x)))
  
  xx <- as.data.frame(ccf_rows) %>%
    magrittr::set_colnames(1:ncol(.)) %>%
    mutate(signature=paste0(strip_title,1:nrow(.))) %>%
    reshape2::melt(.,id=c("signature"))
  
  if (is.na(col)){
    if (nrow(ccf_rows)<=8) fills <-  RColorBrewer::brewer.pal(8, "Set3")[c(1,3:8,2)] else
      fills <- RColorBrewer::brewer.pal(nriw(ccf_rows), "Set3")[c(1,3:8,2,9:ncol(sig))]
  } else {fills <- col}
  
  p1 <- ggplot(xx,aes(y=value,x=variable)) + 
    geom_bar(aes(fill=signature),stat='identity') + 
    scale_fill_manual(values = fills)+ 
    theme_light()+
    #ggtitle(paste0("rank = ",rank,", cancertype = ",type, ", MatrixType = ",MatType )) + 
    theme(legend.title = element_blank(),
          legend.position = "none",
          strip.text.x = element_text(color= "white"),
          panel.grid= element_blank(),
          axis.title.x = element_text(color = "grey20"),
          axis.text.x = element_text(color = "grey20"),
          axis.text.y = element_text(color = "grey20")
          ) +
    facet_grid(rows  = vars(signature),scales="free")+ 
    scale_x_discrete(breaks=c("1","50","100") ,labels=c("0", "0.5", "1"))+
    labs(x="Cancer Cell Fraction",y="",title=title,tag=tag)
  
  g1 <- ggplot_gtable(ggplot_build(p1))
  
  strip_both <- which(grepl('strip-', g1$layout$name))
  
  k <- 1
  for (i in strip_both) {
    j1 <- which(grepl('rect', g1$grobs[[i]]$grobs[[1]]$childrenOrder))
    g1$grobs[[i]]$grobs[[1]]$children[[j1]]$gp$fill <- fills[k]
    #g1$layout$clip[j3] <- "off"
    k <- k+1
  }
  grid::grid.draw(g1)
  return(g1)
}


#' EvoSig_generate: Generate Evolutionary Signatures
#' 
#' This function performs NMF for a given cancer type and generates evolutionary signatures.
#' @name EvoSig_generate
#' @param ccfMat A matrix containing ccf distribution data
#' @param num_sigs The number of signatures suggested by rank estimate analysis
#' @param nrun Run times
#' @return A list containing results, signatures, and exposures
#' @import NMF
#' @importFrom magrittr %>% set_colnames
#' @import dplyr
EvoSig_generate <- function(ccfMat,rank,nrun){
  
  n_sample <- ncol(ccfMat)
  
  # Extracting the first 100 rows and transposing the matrix
  ccf <- t(apply(ccfMat[1:100,],1,as.numeric))
  
  # Preprocess rows with all 0
  index_p <- which(rowSums(ccf)>0)
  index_n <- which(!rowSums(ccf)>0)
  ccf<- ccf[index_p,]
  
  # Run NMFF
  res <- NMF::nmf(ccf,num_sigs,nrun=nrun,.opt='vp4')
  
  # Extracting signatures and exposures
  sig <- matrix(0, nrow = length(index_p)+length(index_n), ncol = ncol(res@fit@W))
  sig[c(index_p),] <- as.data.frame(res@fit@W) %>% set_colnames(paste0("sig_",1:ncol(.)))

  expo <- as.data.frame(t(res@fit@H))
  colnames(expo)[1:num_sigs] <- paste0("sig_",1:num_sigs)
  
  # Output results, signatures, and exposures
  output <- list(
    sig = as.data.frame(sig),
    expo = expo
  )
  
  save(expo,file=paste0(output,type,'/',type,'_',n_sample,"_expo_",Sys.Date(),".RData"))
  save(sig,file=paste0(output,type,'/',type,'_',n_sample,"_sig_",Sys.Date(),".RData"))
  save(res,file=paste0(output,type,'/',type,'_',n_sample,"_res_",Sys.Date(),".RData"))
  return(output)
}


#' Plotting rank estimate
#' @name rank_estimate_plot
#' @param outputFolder folder stores rank estimate files
#' @param rankfilepath rank file path to output
#' @return save nmf results in output folder and plot signature for all cancer types
#' @export
#' @importFrom cowplot plot_grid
#' @importFrom gridExtra arrangeGrob grid.arrange
#' @import dplyr
#' @import ggplot2
#' @import reshape2
rank_estimate_plot <- function(outputFolder,rankfilepath,format) {
  
  typelist <- unique(unlist(lapply(dir(outputFolder),function(x) strsplit(x,"_")[[1]][[1]])))
  rank_blank <- data.frame(cancertype=typelist,rank=NA,rss_suggested_rank=NA)
  
  i = 1
  for (i in 1:length(typelist)) {
    tryCatch({
      type <- typelist[i]
      if (format=="fraction"){
        estimate <- read.csv(paste0(outputFolder,type,"_ccfFractionMatrix.csv")) %>% mutate(Data='normal')
        estimate_random <- read.csv(paste0(outputFolder,type,"_ccfFractionMatrix.random.csv")) %>% mutate(Data='random')
      }
      
      if (format=="count"){
        estimate <- read.csv(paste0(outputFolder,type,"_ccfCountMatrix.csv")) %>% mutate(Data='normal')
        estimate_random <- read.csv(paste0(outputFolder,type,"_ccfCountMatrix.random.csv")) %>% mutate(Data='random')
      }
      
      estimate_rank <- rbind(estimate,estimate_random) %>% filter(rank>2)
      xx <- reshape2::melt(estimate_rank,id=c("rank","Data")) %>% mutate(Measure=NA)
      xx$variable <- as.character(xx$variable)
      
      xx[which(xx$variable=="dispersion"),]$Measure <- "Best fit"
      xx[which(xx$variable=="cophenetic"),]$Measure <- "Consensus"
      xx[which(xx$variable=="sparseness1"),]$Measure <- "Basis"
      xx[which(xx$variable=="sparseness2"),]$Measure <- "Coefficients"
      xx[c(which(xx$variable=="sparseness2"),which(xx$variable=="sparseness1")),]$variable <- "sparseness"
      xx[which(xx$variable=="rss"),]$Measure <- "Consensus"
      xx <- subset(xx,variable %in% c("cophenetic","dispersion","rss","sparseness"))
      
      min_rank = min(xx$rank)
      
      idx = 1:(nrow(estimate)-1)
      rss_decrease <-  min_rank + which((estimate[order(estimate$rank),]$rss[idx]-estimate[order(estimate$rank),]$rss[idx+1]) - (estimate_random[order(estimate_random$rank),]$rss[idx]-estimate_random[order(estimate_random$rank),]$rss[idx+1])<0)[1]
      rank_blank[which(rank_blank$cancertype == type),]$rss_suggested_rank <- rss_decrease 
      
      write.csv(estimate_rank,paste0(outputFolder,type,'_rank_summary.csv'))
      
      command1 <- paste0('g',i,'<-ggplot(data=xx,aes(x=rank,y=value))+ geom_line(aes(lty=Data,color=Measure),cex=1)+geom_point(aes(shape=Data,color=Measure),cex=2)+facet_wrap(~variable,scales="free",nrow=1)+
          labs(title=paste0("Rank estimate for ",type),subtitle=paste0("- rss suggested rank = ",rss_decrease))+
          scale_x_continuous(breaks = min(xx$rank):max(xx$rank))+theme_grey()+ theme(strip.background = element_rect(fill="orange"),strip.text = element_text(colour ="white",size=14))')
      
      command2 <- paste0("print(g",i,")")
      eval(parse(text=command1))
      eval(parse(text=command2))
      
    },error=function(e) print("error!"))
  }
  if (exists("g1")) {
    # output signature plot for all cancer types
    j <- length(typelist) 
    n_file <-  j %% 8
    
    if (n_file==0) {
      n_pdf <- (j %/% 8) 
    } else {
      n_pdf <- (j %/% 8)+1
    }
    
    for (i in 1:n_pdf) {
      if (i != n_pdf) {
        commands1 <- paste0("pdf(file=paste0(outputFolder",",'rank_estimate_g',(8*i-7),'-',8*i,'_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("p",i,"<- plot_grid(paste0('g',(8*i-7):(8*i),collapse=','),align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))") 
      } else {
        commands1 <- paste0("pdf(file=paste0(outputFolder",",'rank_estimate_g',(8*i-8),'-',8*i-8+n_file,'_',Sys.Date(),'.pdf'),height=25,width=12)")
        commands2 <- paste0("p",i," <- plot_grid(paste0('g',(8*i-8):(8*i-8+n_file),collapse=','),align='V',ncol=1,rel_heights = c(0.1,0.1,0.1,0.1,0.1,0.1,0.1,0.1))")
      }
      commands3 <- paste0("print(p",i,")")
      commands4 <- "dev.off()"
      
      eval(parse(text=commands1))
      eval(parse(text=commands2))
      eval(parse(text=commands3))
      eval(parse(text=commands4))
    }
    write.csv(rank_blank,file=rankfilepath)
  }
}


