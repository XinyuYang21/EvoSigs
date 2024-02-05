

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
combine_sig_nmf <- function(sig_folder,output_folder=NA,cancertype){
  
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
  EvoSig_exposure$samplename = ccfMat[,101])
  EvoSig_exposure$samplename = substr(filename$samplename,1,12)
  EvoSig_exposure$samplename = gsub("[.]","-", EvoSig_exposure$samplename)
    
  if (!is.na(output)) {
    if (!dir.exists(output)) {dir.create(output)} 
    save(EvoDynamics_exposure,file=paste0(output,"evoDynamic_exposure.RData"))
  }
  
  return(EvoSig_exposure)
}

#' Perform Hierarchical clustering number estimate
#' @name hc_cluster_test ccfMat ccf matrix for all samples
#' @param combined_sigs combined signatures for all cancer types
#' @param methods clustering method
#' @param distance distance funciton
#' @param min minimum clustering number
#' @param max maximum clustering number
#' @return exposure
#' @import NbClust
hc_cluster_test <- function(combined_sigs,methods,distance,min = 2,max = 10){
  
  getmode <- function(v) {
    uniqv <- unique(v)
    uniqv[which.max(tabulate(match(v, uniqv)))]
  }
  
  tabla = as.data.frame(matrix(ncol = length(distance), nrow = length(methods)))
  names(tabla) = distance
  
  for (j in 2:length(distance))
    for(i in 1:length(methods)){
      tryCatch({
        nb = NbClust::NbClust(data,distance = distance[j],
                     min.nc = min, max.nc = max, 
                     method = "complete", index =methods[i])
        tabla[i,j] = nb$Best.nc[1]
        tabla[i,1] = methods[i]
      },error=function(e) print("error"))
    } 
  
  tabla <- rbind(tabla,c("Most_common",apply(tabla[,2:5],2,getmode)))
  
  return(tabla)
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

#' @name nmf_sig_plot_type
#' @param type cancer type
#' @param MatType matrix type
#' @param input_folder ccf file path
#' @param output output file path
#' @param rank rank
#' @return save nmf results in output folder and plot signature for this type
#' @import NMF
#' @importFrom magrittr %>% set_colnames
#' @import dplyr
nmf_sig_plot_type <- function(type,MatType="fraction",input_folder,output,rank,nrun){
  
  type_path <- paste0(input_folder,type,"/")
  
  if (MatType=="fraction") {
    file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfFractionMatrix_",dir(type_path))])
  } 
  
  if (MatType=="count"){
    file_path <- paste0(input_folder,type,"/",dir(type_path)[grep("ccfCountMatrix_",dir(type_path))])
  }
  
  if (!dir.exists(paste0(output,type))) {
    dir.create(paste0(output,type)) 
  }
  
  load(file=file_path)
  
  #format rank summary file
  if (exists("ccfFractionMatrix")) ccfMat <- ccfFractionMatrix
  if (exists("ccfCountMatrix")) ccfMat <- ccfCountMatrix
  
  n_sample <- ncol(ccfMat)
  
  ccf <- t(apply(ccfMat[1:100,],1,as.numeric))
  
  #preprocess for rows with all 0
  index_p <- which(rowSums(ccf)>0)
  index_n <- which(!rowSums(ccf)>0)
  ccf<- ccf[which(rowSums(ccf)>0),]
  
  #run NMF
  res <- nmf(ccf,rank,nrun=nrun,.opt='vp4')
  
  sig <- as.data.frame(matrix(0,nrow=length(index_p)+length(index_n),ncol=ncol(res@fit@W)))
  
  sig[c(index_p),] <- as.data.frame(res@fit@W) %>% set_colnames(paste0("sig_",1:ncol(.)))
  
  expo <- as.matrix(res@fit@H)
  
  #output sig and expo
  expo <- as.data.frame(t(expo)) 
  colnames(expo)[1:rank] <- paste0("sig_",1:rank)
  
  save(expo,file=paste0(output,type,'/',type,'_',n_sample,"_expo_",Sys.Date(),".RData"))
  save(sig,file=paste0(output,type,'/',type,'_',n_sample,"_sig_",Sys.Date(),".RData"))
  save(res,file=paste0(output,type,'/',type,'_',n_sample,"_res_",Sys.Date(),".RData"))
}

#' Perform nmf for multiple cancer types
#' @name nmf_sig_all_plot
#' @param input_folder ccf file path
#' @param output output file path
#' @param rank_summary rank summary file
#' @return save nmf results in output folder and plot signature for all cancer types
#' @export
#' @import dplyr
#' @import NMF
nmf_sig_all_plot <- function(input_folder,output,rank_summary,MatType="fraction",nrun=100) {
  
  if (!dir.exists(output)) dir.create(output)
  
  lapply(1:nrow(rank_summary),function(x) nmf_sig_plot_type(rank_summary[x,1],input_folder=input_folder,output=output,rank=rank_summary[x,2],MatType=MatType,nrun=nrun))
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

#' Build count matrix for input samples
#' @name Build_post_summary
#' @param ccf_folder ccf files folder
#' @param output output path 
#' @return a data frame containing the summary for all samples in the ccf files
#' @export
Build_post_summary_TCGA <- function(ccf_folder,output=NA,minsample=30){
  
  sample_list <- dir(ccf_folder)
  n_sample = length(dir(input))
  cat("Start building post summary for",n_sample," CCF files \n")
  
  post_summary_analyse_single_ccf <- function(samplename){
    
    ssm <- load_ccf(samplename,input=input)
    ccf <- unique(ssm$ccube_ccf_mean)
    ccf_mean_order <- sort(ccf,decreasing = T)
    Ncluster <- length(ccf)
    
    post_summary <- data.frame(samplename <- samplename) 
    post_summary$Tumor_Sample_Barcode = unique(ssm$Tumor_Sample_Barcode)
    post_summary$n_mutations = nrow(ssm)
    post_summary$ave_depth <- mean(ssm$total_counts)
    post_summary$ccf_01_percentage = mean(ssm$ccube_ccf<=1)
    post_summary$ccf_02_percentage = mean(ssm$ccube_ccf<=2)
    post_summary$Ncluster = Ncluster
    post_summary$purity = unique(ssm$purity)
    post_summary$ccube_purity = ifelse(exists("ssm$ccube_purity"),unique(ssm$ccube_purity),NA)
    post_summary$ccf_mean_cluster1 = ifelse(Ncluster>=1,ccf_mean_order[1],0)
    post_summary$ccf_mean_cluster2 = ifelse(Ncluster>=2,ccf_mean_order[2],0)
    post_summary$ccf_mean_cluster3 = ifelse(Ncluster>=3,ccf_mean_order[3],0)
    post_summary$ccf_mean_cluster4 = ifelse(Ncluster>=4,ccf_mean_order[4],0)
    post_summary$ccf_mean_cluster5 = ifelse(Ncluster>=5,ccf_mean_order[5],0)
    colnames(post_summary) = c("samplename","Tumor_Sample_Barcode","n_mutations","ave_depth","ccf_0-1_percentage","ccf_0-2_percentage","Ncluster","purity","ccube_purity","ccf_mean_cluster1","ccf_mean_cluster2","ccf_mean_cluster3","ccf_mean_cluster4","ccf_mean_cluster5")
    
    cat(">")
    return(post_summary)
  }
  
  post_summary <- do.call(rbind,lapply(sample_list,post_summary_analyse_single_ccf))
  data("TCGAtypefile") 
  cancertype <- TCGAtypefile
  
  if (colnames(cancertype)[1]=="Tumor_sample_barcode") 
    colnames(cancertype)[1] <- "samplename"
  
  post_summary <- suppressWarnings(left_join(post_summary, cancertype,by="samplename"))
  
  if ("Types" %in% colnames(post_summary)) 
    colnames(post_summary)[which(colnames(post_summary)=="Types")] <- "cancertype"
  
  # Output
  if (!is.na(output)) {
    
    if (!dir.exists(output)) dir.create(output,recursive = T)
    write.csv(post_summary,file=paste0(output,"post_summary_",n_sample,"_",Sys.Date(),".csv"))
    save(post_summary,file=paste0(output,"post_summary_",n_sample,"_",Sys.Date(),".RData"))
    
  }
  
  cat("\n Done \n")
  return(post_summary)
}




#' Correlation calculation function 
#' @name cor.table
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @return A table including correlation r,p,n between variables A and variables B
#' @import dplyr
#' @import reshape2
#' @export
cor.table = function(df,var_idx,sig_idx,delete_0=TRUE,min_n=30){
  
  if (delete_0) df[,c(sig_idx,var_idx)] <- apply(df[,c(sig_idx,var_idx)],2,function(x) replace(x,x==0,NA)) # avoid 0 and zero-deviation
  
  summary = c(0,0,0)
  delete variable with less than `min_n` samples
  for (i in var_idx)
    for (j in sig_idx) {
      summary = rbind(summary,c(i,j,length(which(!apply(is.na(df[,c(i,j)]),1,any)))))
    }
  colnames(summary) = c("var_idx","sig_idx","n")
  
  summary = as.data.frame(summary[-1,]) %>% dcast(var_idx~sig_idx,value.var = "n") %>%
    mutate(min=apply(.[,-1],1,min)) %>%
    filter(min>min_n)
  var_idx= summary$var_idx
  
  res <- psych::corr.test(df[,c(var_idx)],df[,c(sig_idx)],method = "spearman",use = "pairwise",adjust="BH",ci=FALSE)
  if (!is.null(res)) {
    r = melt(res['r']) %>% dplyr::rename(r = value) %>% dplyr::select(-L1)
    p = melt(res['p']) %>% dplyr::rename(adj.p = value) %>% dplyr::select(-L1)
    n = melt(res['n']) %>% dplyr::rename(n = value) %>% dplyr::select(-L1)
    
    cor_table = left_join(r,p,by=c("Var1","Var2"))
    if (ncol(n)>1) {
      cor_table %>%
        left_join(n,by=c("Var1","Var2")) %>%
        mutate(n=as.numeric(n),significant=ifelse(adj.p<0.05,1,0))
    } else {
      cor_table %>%
        mutate(n=as.numeric(n)) %>%
        mutate(significant=ifelse(adj.p<0.05,1,0))
    }
  }
}

#' Plot scatterplot for specific measure with evo sig and subtype - 20200915
#' @name scatter_facet
#' @param df exposure file merged with measures of interest
#' @param sig Column indexs of evo signature
#' @param vb Column indexs of all measures in association heatmap in df
#' @param facet name of facet variable consistent with association heatmap
#' @return Scatterplots
#' @import dplyr
#' @import ggpubr
#' @import scales
#' @import ggplot2
#' @importFrom gridExtra grid.arrange
#' @export
#' @example scatter_facet(df=file,vb=var_index,va=sig_idx,var=x,facet_name="cancertype",x_lab="Signature Exposure (count / MB)",legend.position = "",highlight_color="lightsalmon")),ncol=1))
scatter_facet <- function(var,df,vb,var_idx,facet_name="cancertype",cor_table,label.x=NA,label.y=NA,bin.width.x=NA,bin.width.y=NA,
                          x_lab="",legend.position="right",highlight_color="#1A9993FF",zero_delete=TRUE,subtype_panel=FALSE,varwidth=FALSE,log2_scale=FALSE,
                          shape_subtype=FALSE,compare_subtype=FALSE,side_boxplot=TRUE) {
  sig_var = colnames(df)[vb]
  df=  df[,c(var_idx,vb,which(colnames(df) %in% c("subtype",facet_name)))]
  var_idx = which(colnames(df)==var)
  colnames(df)[which(colnames(df)==facet_name)]="facet"
  
  # if (is.na!(cor_table)) {
  #   cor_table = plyr::ddply(df, .(facet) , .fun =cor.table,va=var_idx,vb=vb) 
  # }
  
  cor_table = cor_table %>%
    mutate(p_label=paste0("r =",round(r,2),"\n p = ",round(adj.p,3),"\n n=",n),significant=as.factor(significant),Var2=as.character(Var2),
           group=as.character(group)) %>%
    filter(Var1==var) 
  colnames(cor_table)[1] ="facet"  
  
  df$Var=df[,var_idx]
  
  df_new <- df %>%
    reshape2::melt(id=colnames(df)[-which(colnames(df) %in% sig_var)]) %>%
    dplyr::rename(Var2=variable) %>%
    drop_na("Var") %>%
    mutate(Var=as.numeric(Var),facet=as.character(facet),Var2=as.character(Var2)) %>%
    left_join(cor_table,by=c("Var2","facet")) 
  #left_join(subset(cor_table,facet=="All") %>% 
  #            dplyr::rename(all_p_label= p_label,all_significant=significant) %>% .[,c(2,7:9)],by="Var1")
  
  if (zero_delete) df_new <- subset(df_new,Var!=0 & value!=0)
  
  if (log2_scale) 
    df_new = df_new %>%
    mutate(value=ifelse(value>0,log2(value),value),
           Var=ifelse(Var>0,log2(Var),Var))
  
  fun_median_y <- function(x){
    return(data.frame(y=median(x),label=round(median(x,na.rm=T),2)))}
  
  xmax = max(df_new$value); xmin=min(df_new$value);xmean=mean(df_new$value)
  ymin = min(df_new$Var); ymax=max(df_new$Var);ymean <- mean(df_new$Var)
  
  if (is.na(label.x)) label.x=xmax*0.75
  if (is.na(label.y)) label.y=ymax
  if (is.na(bin.width.x)) bin.width.x= (ymax-ymin)*0.05
  if (is.na(bin.width.y)) bin.width.y= (xmax-xmin)*0.05
  
  scatter_plot_function = function(df) {
    wilcox =  compare_means(Var ~ subtype,subset(df,subtype %in% c("ESig3","ESig4")),group.by='facet',p.adjust.method = "BH")
    df= left_join(df,wilcox,"facet")
    
    if (shape_subtype) 
      p = ggplot(df)+ geom_point(aes(x=value,y=Var,colour = significant,shape=shape),size=0.8)
    else
      p = ggplot(df)+ geom_point(aes(x=value,y=Var,colour = significant),size=0.8)   # y_axis boxplot
    
    if (side_boxplot) {
      if (compare_subtype) 
        p=p+
          geom_boxplot(data=subset(df,subtype %in% c("ESig3","ESig4")),aes(x=xmin,y=Var,fill=subtype),width=bin.width.y,varwidth=varwidth,outlier.size=0.5)+
          geom_boxplot(data=subset(df),aes(x=value,y=ymin),orientation="y",width=bin.width.x,outlier.size=0.5) # x_axis boxplot
      else
        p=p+
          geom_boxplot(data=df,aes(x=xmin,y=Var,fill=subtype),width=bin.width.y,varwidth=varwidth,outlier.size=0.5)+
          geom_boxplot(data=subset(df),aes(x=value,y=ymin),orientation="y",width=bin.width.x,outlier.size=0.5) # x_axis boxplot
    }
    p +
      scale_colour_manual(values=c("grey",highlight_color))+
      geom_smooth(method = "lm",aes(value,Var,color=significant))+
      #scale_x_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))+
      geom_text(aes(x=xmin,y=ymax,label=p.signif),position = "dodge",size=3,vjust = "inward",size=3,check_overlap = TRUE)+ # Add check_overlap, otherwise produce low print quality
      facet_grid(cols=vars(facet),rows=vars(Var2),scale="free")+
      labs(y=var,x=x_lab,title=var)+
      theme_pubclean()+
      geom_text(aes(y=label.y,x=label.x,label=p_label),position = "dodge",size=3,vjust = "inward", hjust = "inward",check_overlap = TRUE)+
      theme(legend.position = legend.position,plot.margin = unit( c(0.2,0,0,0) , units = "lines" ) )+
      scale_fill_manual(values=c("#eccbae","#abddde"))
  }
  
  if (xmin>=0) df_new = df_new %>% filter(value>0)
  if (ymin>=0) df_new = df_new %>% filter(Var>0)
  
  p1 <-   df_new %>% scatter_plot_function()
  
  if (log2_scale) p1 = p1+labs(y=paste0("log2(",var,")"),x=paste0("log2(",x_lab,")"))
  # if (all(df_new$Var>=0) & log2_scale) p1 <- p1+scale_y_continuous(trans='sqrt',breaks=trans_breaks("sqrt",function(x) x^2),labels=trans_format("sqrt",function(x) x^2))
  
  
  return(p1)
  
}
#' Correlation calculation function by group and plot association heatmap - 20200915
#' @name cor_facet
#' @param df exposure file merged with measures of interest
#' @param va Column indexs of variables A in df
#' @param vb Column indexs of variables B in df
#' @param facet name of group variable
#' @param heatmap specify whether plot association heatmap
#' @param title Title for association heatmap
#' @param empty_row_delete whether to delete rows without any significant correlation value
#' @param flip whether to flip facet and x
#' @param keep_all whether to show all facet
#' @param col_low color for correlation r=-1
#' @param col_high color for correlation r=1
#' @return A table including correlation r,p,n between variables A and variables B during group variable and association heatmap (optional)
#' @import dplyr
#' @importFrom plyr ddply
#' @import reshape2
#' @import ggpubr
#' @import ggplot2
#' @export
#' @Example va:var vb:sig

cor_facet = function(df,var_idx,sig_idx,facet,heatmap=FALSE,title="",empty_row_delete=FALSE,flip=FALSE,keep_all=TRUE,col_low="#4a7b94",col_high="#bb5a39",scatter=FALSE,scatter_var=NA,
                     label.x=NA,label.y=NA,bin.width.x=NA,bin.width.y=NA){
  
  var_va = colnames(df)[var_idx];var_vb = colnames(df)[sig_idx]
  df = df[,c(var_idx,sig_idx,which(colnames(df) %in% c(facet,"subtype","samplename")))]
  colnames(df)[which(colnames(df)==facet)] = "group"
  var_idx = which(colnames(df) %in% var_va)
  sig_idx = which(colnames(df) %in% var_vb)
  
  sig_idx = sig_idx[which(apply(df[,sig_idx],2,function(x) length(unique(x))==1)==FALSE)]
  var_idx = var_idx[which(apply(df[,var_idx],2,function(x) length(unique(x))==1)==FALSE)]
  df[,c(var_idx,sig_idx)] = apply(df[,c(var_idx,sig_idx)],2,as.numeric)
  
  if (length(sig_idx)!=0) {
    
    cor_table_all <- cor.table(df,var_idx=var_idx,sig_idx=sig_idx) %>% mutate(group='All')
    
    cor_table <- plyr::ddply(df, .(group) , .fun =cor.table,var_idx=var_idx,sig_idx=sig_idx) %>% 
      rbind(cor_table_all)  %>%
      mutate(r=ifelse(n<30,0,r)) %>%
      mutate(r=replace(r,is.na(r),0),adj.p=replace(adj.p,is.na(adj.p)|is.nan(adj.p),1)) %>% # set correlation with na or p>0.05 invisible
      mutate(r=round(r,2),Var2=as.character(Var2),label=paste0(group," \n (n=",n,")"))
    
    if (scatter) {
      p_scatter = scatter_facet(df=df,va=var_idx,vb=sig_idx,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
                                label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      
      # p_esig3 = scatter_facet(df= df %>% filter(subtype=="ESig3"),va=va,vb=vb,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
      #                         label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      # p_esig4 = scatter_facet(df=df %>% filter(subtype=="ESig4"),va=va,vb=vb,var=scatter_var,facet_name="group",highlight_color="lightsalmon",
      #                         label.x=label.x,label.y=label.y,bin.width.x=bin.width.x,bin.width.y=bin.width.y,subtype_panel=TRUE)
      # p_subtype = p_esig3/p_esig4
    }
    
    
    if (heatmap) {
      cor_table_plot <- cor_table %>% mutate(r=replace(r,adj.p>0.05,0))  %>% 
        left_join( cor_table %>% group_by(Var1,group) %>% dplyr::summarise(n_show=min(n)),by=c("Var1","group")) 
      
      # delete empty rows
      if (keep_all==FALSE)  cor_table_plot <- subset( cor_table_plot,group!="All")
      
      if (empty_row_delete==TRUE) {
        non_empty_gene <- cor_table_plot %>% group_by(Var2) %>% 
          dplyr::summarise(empty_r=all(r==0)) %>% filter(empty_r==FALSE) %>% as.data.frame() %>% .[,1] %>% as.character()
        
        cor_table_plot <- subset(cor_table_plot,Var2 %in% non_empty_gene) 
      }
      
      if (flip==TRUE) {
        
        n_facet <- length(non_empty_gene)
        
        p_heatmap <- cor_table_plot %>% 
          ggplot(aes(y=Var1,x=group,fill=r)) +
          geom_tile()+ geom_text(aes(y=Var1,x=group,label=r),size=4,col='#ffffff')+
          facet_grid(cols=vars(Var2),switch="y",scale="free",space="free")
        
      } else {
        
        p_heatmap <- cor_table_plot %>%
          ggplot(aes(y=Var1,x=Var2,fill=r)) +
          geom_tile()+ 
          geom_text(aes(y=Var1,x=Var2,label=r),size=4,col='#ffffff')
        
        if (length(unique(cor_table_plot$n))<=5) {
          
          p_heatmap <- p_heatmap + facet_grid(cols=vars(label),switch="y",scale="free",space="free")
          
        } else {
          
          n_var1 <- length(unique(cor_table_plot$Var2))
          p_heatmap <- p_heatmap + 
            geom_text(aes(y=Var1,x=n_var1+1,label=n_show),size=3,col="grey50")+
            facet_grid(cols=vars(group),switch="y",scale="free",space="free")+ 
            coord_cartesian( xlim=c(1,n_var1+0.5),clip = "off")
        }
      }
      
      
      
      p_heatmap <- p_heatmap +
        scale_fill_gradient2(low=col_low,mid="#ffffff",high=col_high,midpoint = 0,name="Spearman Correlation")+
        labs(x="",y="",title=title)+ theme_pubclean()+
        theme(axis.text.x = element_text(angle=90,vjust=0.5),plot.margin = unit(c(1, 7, 1, 1), "lines"),legend.position = "right")
      
    } 
    
    # if (heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,p_scatter,p_subtype,non_empty_gene)),return(list(cor_table,p_heatmap,p_scatter,p_subtype)))
    # if (!heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_scatter,p_subtype,non_empty_gene)),return(list(cor_table,p_scatter,p_subtype)))
    # if (heatmap & !scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,non_empty_gene)),return(list(cor_table,p_heatmap)))
    # if (!heatmap & !scatter) return(cor_table)
    
    if (heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,p_scatter,non_empty_gene)),return(list(cor_table,p_heatmap,p_scatter)))
    if (!heatmap & scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_scatter,non_empty_gene)),return(list(cor_table,p_scatter)))
    if (heatmap & !scatter) ifelse(empty_row_delete==TRUE,return(list(cor_table,p_heatmap,non_empty_gene)),return(list(cor_table,p_heatmap)))
    if (!heatmap & !scatter) return(cor_table)
    
    
    
    
  } else {print("the standard deviation is zero")}
}


#' Load multiple format ccf file
#' @name load_ccf
#' @param samplename cancer type
#' @param input 
#' @export
#' @return ssm
load_ccf <- function(samplename,input){
  Check <- ArgumentCheck::newArgCheck()
  suppressWarnings(rm(ssm,res,ccubeRes))
  
  format <- paste0(input,samplename,"/ccubeRes.RData")
  
  if (file.exists(format )) load_ssm = get(load(format4)) 
  return(load_ssm$ssm) 
}

#' create multiple dir
#' @name multi.dir.create
#' @param list list of directory
#' @return create multiple directory
#' @export
multi_dir_create <- function(dirlist){
  for (i in dirlist) {if (!dir.exists(i)) dir.create(i,recursive = T)}
}

#' Unify format of data frame 
#' @name file_format
#' @param df data frame
#' @param samplenamecol column index of samplename
#' @export
file_format <- function(filename=filename,samplenamecol){
  names <- colnames(filename)
  names[samplenamecol] <- "samplename"
  names -> colnames(filename)
  filename$samplename <- substr(filename$samplename,1,12)
  filename$samplename <- gsub("[.]","-",filename$samplename)
  return(filename)
}

sample_format <- function(df,colname="sample"){
  colnames(filename)[samplenamecol] <- "samplename"
  df$samplename <- substr(filename$samplename,1,12)
  df$samplename <- gsub("[.]","-",filename$samplename)
  return(df)
}

#' Merge ssm and cna files
#' @name ParseSnvCnaPcawgFormat
#' @param ssm ssm
#' @param cna cna
#' @export
#' @import dplyr
ParseSnvCnaPcawgFormat <- function (ssm, cna) {
  
  ssm <- ssm %>%
    mutate(chr = substr(chr,4,length(chr)),
           cn_frac = NA,
           major_cn = NA,
           minor_cn = NA,
           mutation_id = NA)
  
  for (jj in seq_len(nrow(cna)) ) {
    cc = cna[jj,]
    
    idx = which(ssm$chr == cc$chromosome &  (ssm$Start_Position >= cc$start & ssm$End_Position <= cc$end) )
    
    if (length(idx) > 0) {
      ssm[idx,] <- ssm[idx,] %>% 
        mutate( major_cn=cc$major_cn,minor_cn =cc$minor_cn, cn_frac = 1)
    }
  }
  
  ssm$mutation_id = paste0(ssm$chr, "_", ssm$Start_Position )
  
  ssm <- ssm %>%
    select(-chr,-Start_Position,-End_Position,-df.n_alt_count,-n_ref_count) %>%
    rename(var_counts='t_alt_count',ref_counts='t_ref_count') %>%
    mutate(total_counts=var_counts+ref_counts,normal_cn=2) %>%
    filter(!is.na(major_cn) ,!is.na(minor_cn),!is.na(cn_frac),major_cn > 0)
  
  return(ssm)
}

#' Customize top strip color
#' @name strip_color
#' @param p plot
#' @param col customized color
#' @param draw whether to display in plots panal
#' @param direction strip location
#' @import ggplot2 
#' @import grid
#' @export
strip_color <- function(p,col1=signature_col,col2=NULL,draw=FALSE,direction='top'){
  p1 <- ggplot_gtable(ggplot_build(p))
  k <- 1
  
  if (direction=='top') strip_col <- which(grepl('strip-t', p1$layout$name))
  if (direction=='bottom') strip_col <- which(grepl('strip-b', p1$layout$name))
  if (direction=='left') strip_col <- which(grepl('strip-l', p1$layout$name))
  if (direction=='right') strip_col <- which(grepl('strip-r', p1$layout$name))
  if (direction=='both'){
    strip_col_t <- which(grepl('strip-t', p1$layout$name))
    strip_col_l <- which(grepl('strip-l', p1$layout$name))
    
    k <- 1
    for (j in strip_col_t) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col1[k]
      k <- k+1
    }
    
    k <- 1
    for (j in strip_col_l) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col2[k]
      k <- k+1
    }
  } else {
    k <- 1
    for (j in strip_col) {
      j1 <- which(grepl('rect', p1$grobs[[j]]$grobs[[1]]$childrenOrder))
      p1$grobs[[j]]$grobs[[1]]$children[[j1]]$gp$fill <- col1[k]
      k <- k+1
    }
  }
  
  if (draw) grid.draw(p1)
  return(p1)
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




