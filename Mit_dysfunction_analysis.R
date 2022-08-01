################################################################################################################
################################################################################################################
#################################### Mitochondrial dysfunction further analyses ################################ 
################################################################################################################
################################################################################################################
################################################################################################################
#README:
# Done by Mar Muniz. Last update 2022.
# This scrip have a serie of functions that will allow to identify genes linked to mitochondrial dysfunction 
# from scRNASeq analyses using Mitoxplorer database 
#for more info:
## MitoExplorer noted as MX:
# http://mitoxplorer.ibdm.univ-mrs.fr/interactome.php
###########################################################################################
# NOTES:
###########################################################################################
# Be sure that you are working with an up to date version of R and bioconductor environment,
# I noticed that running go.gsets function can give an error if the environment is not up to 
# date specifically go.gsets needs to run in R 3.5.1.
# update R --> https://www.r-bloggers.com/updating-r-on-ubuntu/
# The hidden parts of the code are there for QC purposes or to adapt the script to multiple
# analysis, that part of the code is not setup completely and depends on samples so do not
# expect it to work at the first run fro multiple analyses.
############################################################################################

############################################################################################
##################################### Packages needed ######################################
############################################################################################
library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("xlsx");
library("readxl");library("grid");library("gridBase");library("gridExtra");library("biomaRt");
library("sctransform");library("glmGamPoi");  library("patchwork");library("MAST");library("ggrepel");
library("Matrix");library("igraph");library("pagoda2");library("gtools");library("tidyr");library("ggtree"); 
library("ggdendro");library("RColorBrewer");library("openxlsx");library("grid");library("gridBase");
library("colorspace");library("ggpubr");library("gridExtra");library("RColorBrewer");
library("biomaRt");library("stringr");library("topGO");library("Hmisc");library("gtools");
library("corrplot");library("psych");library("correlation");
###
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()
#wd:
#####################################################

wd <- "C:/Users/Documents/Mitochondrial_dysfunction/"
setwd(wd);


#set output directories
f_Mit <- "input//thresholdDEG_1_minus1/" 

f_results <- "/results/";
f_Rdata <- "RData/";
f_allClusters <-"allClusters/"
f_dea_MAST <- "DEA_MAST/"
f_tables <- "tables/";
f_plots <- "plots/";
f_markers <- "markers/";
f_featurePlot <- "featurePlots/"
f_apoe <- "apoe/"
f_markers_p_Cluster_percond <- "perCluster_percond/"
f_markers_p_Cluster_percond_Sex <- "perCluster_percond_PER_Sex/"
f_markers_p_Cluster <- "perCluster/"
f_heatmap <- "heatmap/"
f_donut <- "donutPlot/"
f_dotplot <- "distribution_dotplot/"

#########################################################  ##  #####################################################
dir.create(file.path(wd, f_results,f_Mit ), showWarnings = F);
dir.create(file.path(wd, f_results,f_Mit,f_tables ), showWarnings = F);
dir.create(file.path(wd, f_results,f_Mit,f_plots ), showWarnings = F);
dir.create(file.path(wd, f_results,f_Mit,f_plots,f_heatmap ), showWarnings = F);
dir.create(file.path(wd, f_results,f_Mit,f_plots,f_donut ), showWarnings = F);
dir.create(file.path(wd, f_results,f_Mit,f_plots,f_dotplot ), showWarnings = F);
dir.create(file.path(wd,f_results,f_corr), showWarnings = F);
dir.create(file.path(wd,f_results,f_corr,f_linearCorr), showWarnings = F);
dir.create(file.path(wd,f_results,f_corr,f_heatmapCorr), showWarnings = F);
#########################################################  ##  #####################################################



################################################################################################################
################################################################################################################
# A) Use all our expressed genes slightly altered or more storngly altered based on log2FC expression levels (DEGs), 
# B) find all the genes in the mitoExplorerDB that are DEGs
# C) Then, together and by MX categories run the dotplot and heatmap in only the DEGs.
################################################################################################################
################################################################################################################

############################################################
##A  Reading mitoXplorer DB files and getting the genes  ####
############################################################
#download the files from Mitoxplorer interactome webpage work with the last version
input_mitoXplorer_hsa<- as.data.frame(read_excel(paste0("C:/Users/Documents/mitoExplorerDB/input/05272022/humanInteractome/human_gene_function.xlsx"),sheet =1, col_names =T))  
input_mitoXplorer_hsa <- unique(input_mitoXplorer_hsa[,c("ENSG_ID","mito_process","gene_name","gene_function","ENSG_name")]);
colnames(input_mitoXplorer_hsa)[c(1:5)]<- c("ensembl_gene_id","Mitochondrial_process","full_gene_name","gene_function","Gene")
length(unique(input_mitoXplorer_hsa$full_gene_name)) # Nb unique genes associated with Mit act: 1229 form v2020 to version 2022 no change

input_mitoXplorer_hsa$Group <- input_mitoXplorer_hsa$Mitochondrial_process
input_mitoXplorer_hsa <- left_join(input_mitoXplorer_hsa,colDonut_Mit ,by="Group")
write.xlsx(input_mitoXplorer_hsa, file = paste0(wd,f_results,f_Mit,f_tables, "allMitochondria_genes_mitoXplorer_colour.xlsx"));



#all MY=TD genes are unique for each category
MTDcategories_perGene <- unique(input_mitoXplorer_hsa[,c("Gene","Mitochondrial_process","gene_function")])
  MTDcategories_perGene$Id <- paste0(MTDcategories_perGene$Gene,"___" ,MTDcategories_perGene$gene_function);
  MTDcategories_perGene_aggregateDf <- unique(aggregate( Mitochondrial_process ~ Id  , data = MTDcategories_perGene, paste, collapse = "// "));
  MTDcategories_perGene_aggregateDf<- separate(MTDcategories_perGene_aggregateDf, Id, c("Gene" ,"gene_function"  ),sep="___");
  MTDcategories_perGene_aggregateDf$nbMTDs <- str_count(MTDcategories_perGene_aggregateDf$Mitochondrial_process, "//") +1;

#DEGs per cluster per analysis
############################################################

# More deregulated genes based in adj pval<0.05,-1 < log2Fc>1
load( file = paste0(wd,f_Mit,"DEGs_mit_theshold1_minus1.RData"));  #DEGS_percond_fem, DEGS_percond_male,DEGS_percond

# More permissive analysis, getting also genes slighly altered that can be contributing to the deregulation
#  genes based in adj pval<0.05,-0.32 < log2Fc>0.25
load( file = paste0(wd,f_Mit,"DEGs_logsFC025_minus032_vF.RData"));  #DEGS_percond_fem_fc12, DEGS_percond_male_fc12,DEGS_percond_fc12

############################################################
## Selecting MTD genes that are DEGs  ########
############################################################

###################################
mit_percond <- unique( left_join(input_mitoXplorer_hsa, DEGs_per_cond ,by=c("Gene")));
mit_percond <- mit_percond[!is.na(mit_percond$cluster ),]
mit_percond_FC12 <- unique( left_join(input_mitoXplorer_hsa, DEGS_percond_fc12 ,by=c("Gene")));
mit_percond_FC12 <- mit_percond_FC12[!is.na(mit_percond_FC12$cluster ),]

##################################################################################################################################
##################################################################################################################################

# 2) for sex specific analysis
###################################
mit_Fem_percond <- unique( left_join(input_mitoXplorer_hsa, DEGS_percond_fem ,by=c("Gene")));
mit_Fem_percond <- mit_Fem_percond[!is.na(mit_Fem_percond$cluster ),]
mit_Male_percond <- unique( left_join(input_mitoXplorer_hsa, DEGS_percond_male ,by=c("Gene")));
mit_Male_percond <- mit_Male_percond[!is.na(mit_Male_percond$cluster ),]

# 2) for sex specific analysis
###################################
mit_Fem_percond_FC12 <- unique( left_join(input_mitoXplorer_hsa, DEGS_percond_fem_fc12 ,by=c("Gene")));
mit_Fem_percond_FC12 <- mit_Fem_percond_FC12[!is.na(mit_Fem_percond_FC12$cluster ),]
mit_Male_percond_FC12 <- unique( left_join(input_mitoXplorer_hsa, DEGS_percond_male_fc12 ,by=c("Gene")));
mit_Male_percond_FC12 <- mit_Male_percond_FC12[!is.na(mit_Male_percond_FC12$cluster ),]

res<- list(LA=mit_percond, Fem=mit_Fem_percond, Male=mit_Male_percond)
write.xlsx(res, file = paste0(wd,f_results,f_Mit,f_tables, "allMitochondria_DEGs_found_thresholds_1_minus1.xlsx"));

save(mit_percond,mit_Fem_percond,mit_Male_percond, file = paste0(wd,f_results,f_Rdata,f_allClusters,"DEGs_mit_theshold1_minus1.RData"));
#############################

############################################################
##  3) percentages of genes per Mit processes per model in 
##############  all  DEGs  #################
############################################################
#input files
####################
categories <- unique(input_mitoXplorer_hsa$`Mitochondrial_process`)

#stats
MTcategories_stats.function <- function(inputDf,nameOutput){
  inputDf <- inputDf[which(inputDf$Group=="Significant & log2FC"),];
  DEGs_percond_per_category_MitoXplorer_gage_genes <- as.data.frame(unique(inputDf[,c('Gene','Mitochondrial_process','cellType_res05')]) %>% dplyr::group_by(Mitochondrial_process,cellType_res05 )  %>% dplyr::count());
  DEGs_percond_per_sense_category_MitoXplorer_gage_genes <- as.data.frame(unique(inputDf[,c('Gene','Mitochondrial_process','cellType_res05','Regulation')]) %>% dplyr::group_by(Mitochondrial_process,cellType_res05,Regulation )  %>% dplyr::count());

  DEGs_percond_per_cluster_category_MitoXplorer_gage_genes <- as.data.frame(unique(inputDf[,c('Gene','Mitochondrial_process','cellType_res05')]) %>% dplyr::group_by(cellType_res05, Mitochondrial_process)  %>% dplyr::count());
  DEGs_percond_sense_per_cluster_category_MitoXplorer_gage_genes <- as.data.frame(unique(inputDf[,c('Gene','Mitochondrial_process','cellType_res05','Regulation')]) %>% dplyr::group_by(cellType_res05, Mitochondrial_process,Regulation)  %>% dplyr::count());

  #stats wuth genes aggregated
  inputDf_annot <- unique(inputDf[,c('Gene','Mitochondrial_process','Regulation','cellType_res05')])
  inputDf_annot <-  as.data.frame(inputDf_annot %>% dplyr::group_by(Mitochondrial_process,Regulation))

  inputDf_annot$Id <- paste0(inputDf_annot$cellType_res05,"___" ,inputDf_annot$Mitochondrial_process,"___",inputDf_annot$Regulation);
  inputDf_annot$Id2 <- paste0(inputDf_annot$cellType_res05,"___" ,inputDf_annot$Gene,"___",inputDf_annot$Regulation);
  inputDf_annot$Id3 <- paste0(inputDf_annot$cellType_res05,"___" ,inputDf_annot$Gene);

  aggregateDf <- unique(aggregate( Gene ~ Id  , data = inputDf_annot, paste, collapse = ", "));
  aggregateDf2 <- unique(aggregate( Mitochondrial_process ~ Id2  , data = inputDf_annot, paste, collapse = ", "));
  aggregateDf3 <- unique(aggregate( Mitochondrial_process ~ Id3  , data = inputDf_annot, paste, collapse = ", "));

  removeDupsSpecialBarChar.function <- function(input,column,nameOutputCollapse){
    resultsList <- list();
    for (j in 1:length(input[,column])){
      resultsList[[j]] <- paste(unique(trimws(unlist(strsplit(input[j,column],split="// ",fixed=F,perl=T)))),collapse=", ");
      j=j+1
    };
    input[,column]<-as.data.frame(unlist(resultsList))[,1];
    assign(paste0(nameOutputCollapse),input, envir=parent.frame());
  };

  removeDupsSpecialBarChar.function(aggregateDf,'Gene',"aggregateDf");
  removeDupsSpecialBarChar.function(aggregateDf2,'Mitochondrial_process',"aggregateDf2");
  removeDupsSpecialBarChar.function(aggregateDf3,'Mitochondrial_process',"aggregateDf3");

  aggregateDf<- separate(aggregateDf, Id, c("cellType_res05" ,"Mitochondrial_process","Regulation"  ),sep="___");
  aggregateDf2<- separate(aggregateDf2, Id2, c("cellType_res05" ,"Gene","Regulation"  ),sep="___");
  aggregateDf3<- separate(aggregateDf3, Id3, c("cellType_res05" ,"Gene"  ),sep="___");
  aggregateDf$nbGenes <- str_count(aggregateDf$Gene, ", ") +1;
  aggregateDf2$nbMTcategories <- str_count(aggregateDf2$Mitochondrial_process, ", ") +1;
  aggregateDf3$nbMTcategories <- str_count(aggregateDf3$Mitochondrial_process, ", ") +1;

  results_MitoActivity <- list(input=inputDf,DEGsName_perSenseMitCat=aggregateDf,sharedDEGs_perCatSense=aggregateDf2,sharedDEGs_perCat=aggregateDf3, DEGs_perMitCat=DEGs_percond_per_category_MitoXplorer_gage_genes,DEGs_perMitCat_perREG=DEGs_percond_per_category_MitoXplorer_gage_genes ,
  DEGs_perCluster_MitCat=DEGs_percond_per_cluster_category_MitoXplorer_gage_genes,DEGs_perCluster_MitCatReg=DEGs_percond_sense_per_cluster_category_MitoXplorer_gage_genes)

  write.xlsx(results_MitoActivity, file = paste0(wd,f_results,f_Mit,f_tables,nameOutput, "_Mt_Results_stats.xlsx"));

};
MTcategories_stats.function(mit_percond,"percond")
MTcategories_stats.function(mit_Fem_percond,"Fem_percond")
MTcategories_stats.function(mit_Male_percond,"Male_percond")


MTcategories_stats.function(mit_percond_FC12,"cond_FC12")
MTcategories_stats.function(mit_Fem_percond_FC12,"Fem_percond_FC12")
MTcategories_stats.function(mit_Male_percond_FC12,"Male_percond_FC12")

############################################################
##  4) donut PLOT 
############################################################

#creating input file for donut plot with percentages and colors included
#1. Creating the col palette
 colDonut_Mit <- data.frame(Group=as.character(categories), stringsAsFactors=FALSE)
 preferred.order <- data.frame(Group=as.character(categories), stringsAsFactors=FALSE)
 colDonut_Mit$colDonut <-c("#00C5CD","#008B8B","#AFEEEE","#1E90FF", "#7FFFD4",
      "#FFC0CB","#FF1493","#b30000","#145214", "#FF6685","#ffc400","#080404",
      "#FFFF00","#ff6600","#996600","#913a11","#4B0082","#ecb3ff","#BA55D3",
      "#FF0000","#FF6A6A","#800000","#96784d","#8B008B", "#9aff91","#6d72b0",
      "#66676e","#b300b3","#A9A9A9","#ffaf00","#ffaf8a","#d2afff","#adafff","#59afff","#1ac48f","#086f1b","#823661","#fa59f1")

donut.inputFile.function <- function(inputDf,colDonut,nameOutput){
  outputDfSumup <- list();
  inputDf <- inputDf[which(inputDf$Group=="Significant & log2FC"),];

  for (i in 1:length(unique(inputDf$cellType_res05))){
    print(i)
    df.tmp <- as.data.frame(inputDf[which(inputDf$cellType_res05==unique(inputDf$cellType_res05)[i]),]);
    df <- unique(df.tmp[,c('Gene','Mitochondrial_process','cellType_res05')])
    
    dfDonut1 <-  as.data.frame(df %>% dplyr::group_by(Mitochondrial_process) %>% dplyr::count())
    dfDonut1$total <- sum(length(unique(df[,1])));
    dfDonut1$perc <-round((dfDonut1$n/dfDonut1$total)*100,digits=2)
    dfDonut1$Group <- dfDonut1$Mitochondrial_process
    dfDonut1f <- left_join(dfDonut1, colDonut, by="Group");
    dfDonut1f<- dfDonut1f[order(match(dfDonut1f[,1],colDonut[,1])),];
    outputDfSumup[[i]] <-dfDonut1f;
    names(outputDfSumup)[i] <- unique(inputDf$cellType_res05)[i];
    i=i+1;
  };
  assign(paste0("sumupGroupDonutSumup_",nameOutput), outputDfSumup,.GlobalEnv);

};


donut.inputFile.function(mit_percond,colDonut_Mit, "percond"); #ignore the warnings
donut.inputFile.function(mit_Fem_percond,colDonut_Mit, "Fem_percond"); #ignore the warnings
donut.inputFile.function(mit_Male_percond,colDonut_Mit, "Male_percond"); #ignore the warnings


donut.inputFile.function(mit_percond_FC12 ,colDonut_Mit, "cond_FC12"); #ignore the warnings
donut.inputFile.function(mit_Fem_percond_FC12 ,colDonut_Mit, "Fem_percond_FC12"); #ignore the warnings
donut.inputFile.function(mit_Male_percond_FC12,colDonut_Mit, "Male_percond_FC12"); #ignore the warnings



donut.inputFile_ordered_byPercentage.function <- function(inputDf,colDonut,nameOutput){
  outputDfSumup_ordered <- list();
  inputDf <- inputDf[which(inputDf$Group=="Significant & log2FC"),];

  for (i in 1:length(unique(inputDf$cellType_res05))){
    print(i)
    df.tmp <- as.data.frame(inputDf[which(inputDf$cellType_res05==unique(inputDf$cellType_res05)[i]),]);
    df <- unique(df.tmp[,c('Gene','Mitochondrial_process','cellType_res05')])
    #name <- as.character(unique(df[,"cellType_res05"]))
    dfDonut1 <-  as.data.frame(df %>% dplyr::group_by(Mitochondrial_process) %>% dplyr::count())
    dfDonut1$total <- sum(length(unique(df[,1])));
    dfDonut1$perc <-round((dfDonut1$n/dfDonut1$total)*100,digits=2)
    dfDonut1$Group <- dfDonut1$Mitochondrial_process
    dfDonut1f <- left_join(dfDonut1, colDonut, by="Group");
    dfDonut1f<- dfDonut1f[order(match(dfDonut1f[,1],colDonut[,1])),];
        dfDonutff<- dfDonut1f[mixedorder(dfDonut1f$perc, decreasing=TRUE),];
        print(head(dfDonutff))
        outputDfSumup_ordered[[i]] <-dfDonutff;
        names(outputDfSumup_ordered)[i] <-  as.character(unique(df$cellType_res05));
  };
      assign(paste0("sumupGroupDonutSumup_ordered_",nameOutput), outputDfSumup_ordered,.GlobalEnv);

};
donut.inputFile_ordered_byPercentage.function(mit_percond,colDonut_Mit, "percond"); #ignore the warnings
donut.inputFile_ordered_byPercentage.function(mit_Fem_percond,colDonut_Mit, "Fem_percond"); #ignore the warnings
donut.inputFile_ordered_byPercentage.function(mit_Male_percond,colDonut_Mit, "Male_percond"); #ignore the warnings

donut.inputFile_ordered_byPercentage.function(mit_percond_FC12,colDonut_Mit, "cond_FC12"); #ignore the warnings
donut.inputFile_ordered_byPercentage.function(mit_Fem_percond_FC12,colDonut_Mit, "Fem_percond_FC12"); #ignore the warnings
donut.inputFile_ordered_byPercentage.function(mit_Male_percond_FC12,colDonut_Mit, "Male_percond_FC12"); #ignore the warnings



donuts_plot <- function(
                        panel = runif(3), # counts
                        pctr = c(.5,.2,.9), # percentage in count
                        legend.condbel='',
                        cols = c('chartreuse', 'chocolate','deepskyblue'), # colors
                        outradius = 1, # outter radius
                        radius = .7,   # 1-width of the donus 
                        add = F,
                        innerradius = .5, # innerradius, if innerradius==innerradius then no suggest line
                        legend = F,
                        pilabels=F,
                        legend_offset=.25, # non-negative number, legend right position control
                        borderlit=c(F,F,T,T)
                        ){
    par(new=add,mar=c(14,14,14,16))
    if(sum(legend.condbel=='')>=1) legend.condbel=paste("Series",1:length(pctr))
    if(pilabels){
        pie(panel, col=cols,border = borderlit[1],labels = legend.condbel,radius = outradius)
    }
    panel = panel/sum(panel)

    pctr2= panel*(1 - pctr)
    pctr3 = c(pctr,pctr)
    pctr_indx=2*(1:length(pctr))
    pctr3[pctr_indx]=pctr2
    pctr3[-pctr_indx]=panel*pctr
    cols_fill = c(cols,cols)
    cols_fill[pctr_indx]='white'
    cols_fill[-pctr_indx]=cols
    par(new=TRUE,mar=c(14,14,14,16))
    pie(pctr3, col=cols_fill,border = borderlit[2],labels = '',radius = outradius)
    par(new=TRUE,mar=c(14,14,14,16))
    pie(panel, col='white',border = borderlit[3],labels = '',radius = radius)
    par(new=TRUE,mar=c(14,14,14,16))
    pie(1, col='white',border = borderlit[4],labels = '',radius = innerradius)
    if(legend){
        par(mar=c(6,13,32,6), xpd=TRUE)
        legend("topright",inset=c(-legend_offset,0),pt.cex=1.2,legend=legend.condbel, cex=0.8,pch=rep(15,'.',length(pctr)), 
               col=cols,bty='n')
    }
    par(new=FALSE)
}
## col- > subcor(change hue/alpha)
subcolors <- function(.dta,main,mainCol){
    tmp_dta = cbind(.dta,1,'col')
    tmp1 = unique(.dta[[main]])
    for (i in 1:length(tmp1)){
        tmp_dta$"col"[.dta[[main]] == tmp1[i]] = mainCol[i]
    }
    u <- unlist(by(tmp_dta$"1",tmp_dta[[main]],cumsum))
    n <- dim(.dta)[1]
    subcol=rep(rgb(0,0,0),n);
    for(i in 1:n){
        t1 = col2rgb(tmp_dta$col[i])/256
        subcol[i]=rgb(t1[1],t1[2],t1[3],1/(1+u[i]))
    }
    return(subcol);
};

#modify pt.cex for the symbols size, the mar parameters for plot and legend
######################
#1.following my chosen order in the groups of pathways  organization:
f_donut<- "DonutPlot/"
dateDonut<- "010122"
plot_donut.function <- function(inputDfSumup,elem,date,nameOutput){
  i=elem
  df <- inputDfSumup[[i]];
  nameDf <- names(inputDfSumup)[i];
  totalPaths <- as.character(unique(df$total));
  ngroups <- nrow(df);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut), showWarnings = F);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut), showWarnings = F);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut,"/",nameOutput), showWarnings = F);
  pdf(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut,"/",nameOutput,"/",nameOutput,"_", nameDf,"_",totalPaths,".pdf"),width = 8.27, height=12.69);  #8.27 × 11.69)
  donuts_plot(df[,4],rep(1,ngroups),paste0(df$Group," ",df$perc," %"),
        cols=df$colDonut,pilabels=F,legend=T,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,ngroups) )  
  dev.off();
  };

######################
#2. ordered groups by higher % of pathways
plot_donut_ordered.function <- function(inputDfSumup_ordered,elem,date,nameOutput){
  i=elem
  df <- inputDfSumup_ordered[[i]];
  nameDf <- names(inputDfSumup_ordered)[i];
  #nameDf <- gsub(".*_","",names(inputDfSumup_ordered)[i]);
  totalPaths <- as.character(unique(df$total));
  ngroups <- nrow(df);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut), showWarnings = F);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut), showWarnings = F);
  dir.create(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut,"/",nameOutput), showWarnings = F);
  pdf(paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut,"/",nameOutput,"/",nameOutput,"_", nameDf,"_",totalPaths,"_ord.pdf"),width = 8.27, height=12.69);  #8.27 × 11.69)
  donuts_plot(df[,4],rep(1,ngroups),paste0(df$Group," ",df$perc," %"),
        cols=df$colDonut,pilabels=F,legend=T,legend_offset=1.3,
        outradius = 0.9,radius = .5,innerradius=.2,add=T,
        borderlit = rep(F,ngroups) )  
dev.off();
};


#############@@@@
#A) plot # the ordered doesnt print all the donuts
#############@@@@


##

for (i in 1: length(names(sumupGroupDonutSumup_ordered_percond))){
       plot_donut_ordered.function(sumupGroupDonutSumup_ordered_percond,i,date,"percond");

      i = i+1
};
##

#
for (i in 1: length(names(sumupGroupDonutSumup_ordered_Male_percond))){
      plot_donut_ordered.function(sumupGroupDonutSumup_ordered_Male_percond,i,date,"Male_percond");

      i = i+1
};


for (i in 1: length(names(sumupGroupDonutSumup_ordered_Fem_percond))){
      plot_donut_ordered.function(sumupGroupDonutSumup_ordered_Fem_percond,i,date,"Fem_percond");

      i = i+1
};


#############@@@@
# B) creating table input for dotplot comparison of the % of pathways one each groups
#############@@@@
create_table_DonutPlot.function <- function(allPathsList_sumupGroupDonutSumup_ordered, nameOutput){
  for (i in 1: length(names(allPathsList_sumupGroupDonutSumup_ordered))){
        df <-allPathsList_sumupGroupDonutSumup_ordered[[i]]
        df <- df[,c('Group','perc','colDonut')]
        colnames(df)[2]<- names(allPathsList_sumupGroupDonutSumup_ordered)[i]
        if (i==1){
           tableDotplot_percPaths.input <-df   
        } else {
           tableDotplot_percPaths.input <-full_join(tableDotplot_percPaths.input,df, by=c("Group","colDonut"))
           tableDotplot_percPaths.input[is.na(tableDotplot_percPaths.input) ] <- 0 
        };
        i = i+1
  };
  
  write.xlsx(tableDotplot_percPaths.input, file = paste0(wd,f_results,f_Mit,f_plots,f_donut,dateDonut,"/",nameOutput,"/",nameOutput,".xlsx"));
  assign(paste0("donutPlot_table_",nameOutput), tableDotplot_percPaths.input,.GlobalEnv);

};
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_percond, "percond");
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_Fem_percond, "Fem_percond");
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_Male_percond, "Male_percond");


create_table_DonutPlot.function(sumupGroupDonutSumup_percond_FC12, "DEGs_percond_FC12");
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_percond_FC12, "DEGs_percond_ordered_FC12");
create_table_DonutPlot.function(sumupGroupDonutSumup_Fem_percond_FC12, "DEGs_Fem_percond_FC12");
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_Fem_percond_FC12, "DEGs_Fem_percond_ordered_FC12");
create_table_DonutPlot.function(sumupGroupDonutSumup_Male_percond_FC12, "DEGs_Male_percond_FC12");
create_table_DonutPlot.function(sumupGroupDonutSumup_ordered_Male_percond_FC12, "DEGs_Male_percond_ordered_FC12");

############################################################################################################

log2FC_matrixInput_4heatmap.function <-  function(data,nameOutput ) {
  #####################################
  #####################################
  resList <- list();
  data <- data[which(data$Group=="Significant & log2FC"),]
  for (i in 1:length(unique(data$Mitochondrial_process))){
    process <- unique(data$Mitochondrial_process)[i];
    df <- unique(data[which(data$Mitochondrial_process==unique(data$Mitochondrial_process)[i]),]);
    df<- unique(df[,c(grep("Gene" ,colnames(df)),grep("avg_log2FC" ,colnames(df)),grep("cellType_res05" ,colnames(df)))]);
    df <- df %>% spread(cellType_res05,avg_log2FC );
    df[is.na(df)] = 0
    names4Rows <- df[,1];
    log2Matrix <- as.matrix(df[,-1]);     
    colNamesDf <- gsub("avg_log2FC","",gsub("_"," ",colnames(df)));
    colNamesDf <- colNamesDf[-1];
    colnames(log2Matrix) <- colNamesDf;
    rownames(log2Matrix) <- names4Rows;
    resList[[i]] <- log2Matrix;
    names(resList)[i] <- process;
    i <- i+1;

  };
  assign(paste0("resList_",nameOutput), resList,.GlobalEnv);
};


log2FC_matrixInput_4heatmap.function(mit_percond,"percond" );
log2FC_matrixInput_4heatmap.function(mit_Fem_percond,"Fem_percond" );
log2FC_matrixInput_4heatmap.function(mit_Male_percond,"Male_percond" );

log2FC_matrixInput_4heatmap.function(mit_percond_FC12,"DEGs_percond_FC12" );
log2FC_matrixInput_4heatmap.function(mit_Fem_percond_FC12,"DEGs_Fem_percond_FC12" );
log2FC_matrixInput_4heatmap.function(mit_Male_percond_FC12,"DEGs_Male_percond_FC12" );


#ok for the inputs lets go for the heatmaps

# A) per individual datatset
###################### colour palette  ###########################
col1 <- colorRampPalette(rev(brewer.pal(9, "RdPu")))(100)  # We choose the colours for the heatmap.
pink  <-   colorRampPalette(brewer.pal(9, "PiYG"))(9)[1:4]
pur  <-   colorRampPalette(brewer.pal(9, "PRGn"))(12)[1:4]
colGroup.tmp1 <- append("#8B0A50",pink)
colGroup.tmp2 <-append(rev(pur),"#35133a")
colGroup <- append(colGroup.tmp1,"#f9f2fa")
colGroup <- append(colGroup, "#f9f2fa")
colGroup <- append(colGroup,colGroup.tmp2)
colGroup1 <- append(colGroup.tmp1,"#f9f2fa")
colGroup1 <- append(colGroup1,colGroup.tmp2)
#################################################################
space<-"                                 "

f_dateHeatmap<- "70122/"
f_heatmap<- "heatmap/"
dir.create(file.path(wd,f_results,f_Mit,f_plots,f_heatmap), showWarnings = F);
dir.create(file.path(wd,f_results,f_Mit,f_plots,f_heatmap,f_dateHeatmap), showWarnings = F);


heatmapPLot.function<- function(inputList,nameSave){
pdf(paste0(wd,f_results,f_Mit,f_plots,f_heatmap,f_dateHeatmap,nameSave,"_Mit_heatmap.pdf") , width =10, height =12 );
  par(cex.main=0.7)
  found <- lapply(inputList,function(x) length(x));
  found <- names(found[which(found>1)]); 
  inputListF<- inputList[found];
  
  for (i in 1: length(names(inputListF))){
    m <- inputListF[[i]];
    nameMit <- names(inputListF)[i];
    if (nrow(m)<2){
      m <- rbind(m,0) 
      rownames(m)[2] <- ""
      gplots::heatmap.2(m, trace="none", Colv = F, dendrogram="row", scale="none",  density.info="none", lhei = c(1.5,6), keysize =1.5,col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);
      gplots::heatmap.2(m, trace="none",Colv = T,  dendrogram="colum", scale="none",density.info="none", lhei = c(1.5,6), keysize =1.5,  col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);
      gplots::heatmap.2(m, trace="none",  dendrogram="both", scale="none",density.info="none", lhei = c(1.5,6), keysize =1.5,  col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);

    } else {
      gplots::heatmap.2(m, trace="none", Colv = F, dendrogram="row", scale="none", density.info="none", lhei = c(1.5,6), keysize =1.5, col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);
      gplots::heatmap.2(m, trace="none",Colv = T,  dendrogram="colum", scale="none",  density.info="none", lhei = c(1.5,6), keysize =1.5,col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);
      gplots::heatmap.2(m, trace="none",  dendrogram="both", scale="none", density.info="none", lhei = c(1.5,6), keysize =1.5, col = colGroup1, main= paste0("DEGs involded in \n", nameMit), margin=c(26,8), breaks=c(-6,-2,-1.2,-0.9,-0.5,-0.1,+0.1,+0.5,+0.9,+1.2,+2,6), cexCol=1,cexRow=0.9);
    };
    i = i+1;

  };
  dev.off();

};


#1 run for the DEGs strict log2fc values
nameExp <- "DEGs: FDR<0.05 -1 <log2FC> 1";
nameSave <- "Male_percond";
inputList<- resList_Male_percond
heatmapPLot.function(resList_Male_percond,"Male_percond");

nameSave <- "Fem_percond";
inputList<- resList_Fem_percond
heatmapPLot.function(resList_Fem_percond,"Fem_percond" );


nameSave <- "DEGs_perCond";
inputList<- resList_percond
heatmapPLot.function(resList_percond,"DEGs_perCond");


##2 run fro the DEGs more permissive thresholds of log2fc

nameExp <- "DEGs: FDR<0.05 -0.32 <log2FC>0.25";
nameSave <- "DEGs_Male_percond_FC12";
inputList<- resList_DEGs_Male_percond_FC12
heatmapPLot.function(resList_DEGs_Male_percond_FC12,"DEGs_Male_percond_FC12");

nameSave <- "DEGs_Fem_percond_FC12";
inputList<- resList_DEGs_Fem_percond_FC12
heatmapPLot.function(resList_DEGs_Fem_percond_FC12,"DEGs_Fem_percond_FC12" );


nameSave <- "DEGs_percond_FC12";
inputList<- resList_DEGs_percond_FC12
heatmapPLot.function(resList_DEGs_percond_FC12,"DEGs_percond_FC12");




##1) stats tables #DEA_thresholds adj p-value< 0.05 -1<log2(FC)> 1
#####################################################################

mit_percond_Cutoofs <-mit_percond[which(mit_percond$p_val_adjPlot<0.05),]
mit_Fem_percond_Cutoofs <-mit_Fem_percond[which(mit_Fem_percond$p_val_adjPlot<0.05),]
mit_Male_percond_Cutoofs <-mit_Male_percond[which(mit_Male_percond$p_val_adjPlot<0.05),]

mit_percond_Cutoofs <-bind_rows(mit_percond_Cutoofs[which(mit_percond_Cutoofs$avg_log2FC>1),],mit_percond_Cutoofs[which(mit_percond_Cutoofs$avg_log2FC < -1),])
mit_Fem_percond_Cutoofs <-bind_rows(mit_Fem_percond_Cutoofs[which(mit_Fem_percond_Cutoofs$avg_log2FC>1),],mit_Fem_percond_Cutoofs[which(mit_Fem_percond_Cutoofs$avg_log2FC < -1),])
mit_Male_percond_Cutoofs <-bind_rows(mit_Male_percond_Cutoofs[which(mit_Male_percond_Cutoofs$avg_log2FC>1),],mit_Male_percond_Cutoofs[which(mit_Male_percond_Cutoofs$avg_log2FC < -1),])

length(unique(mit_percond_Cutoofs$Gene)) 
length(unique(mit_Fem_percond_Cutoofs$Gene)) 
length(unique(mit_Male_percond_Cutoofs$Gene)) 
DEGs_stats <- as.data.frame(mit_percond_Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 
DEGs_statsfem <- as.data.frame(mit_Fem_percond_Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 
DEGs_statsmale <- as.data.frame(mit_Male_percond_Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 
resDegs <- list( LA=DEGs_stats, fem=DEGs_statsfem,male=DEGs_statsmale)
write.xlsx(resDegs, file = paste0(wd,f_results,f_Mit,f_tables,"Mit_DEGs_Stats2.xlsx"));



#2) stats tables #DEGs: DEA_thresholds adj p-value< 0.05 -0.32<log2(FC)> 0.25 
#####################################################################
mit_percond_FC12Cutoofs <-mit_percond[which(mit_percond$p_val_adjPlot<0.05),]
mit_Fem_percond_FC12Cutoofs <-mit_Fem_percond[which(mit_Fem_percond$p_val_adjPlot<0.05),]
mit_Male_percond_FC12Cutoofs <-mit_Male_percond[which(mit_Male_percond$p_val_adjPlot<0.05),]

mit_percond_FC12Cutoofs <-bind_rows(mit_percond_FC12Cutoofs[which(mit_percond_FC12Cutoofs$avg_log2FC>0.25),],mit_percond_FC12Cutoofs[which(mit_percond_FC12Cutoofs$avg_log2FC < -0.32),])
mit_Fem_percond_FC12Cutoofs <-bind_rows(mit_Fem_percond_FC12Cutoofs[which(mit_Fem_percond_FC12Cutoofs$avg_log2FC>0.25),],mit_Fem_percond_FC12Cutoofs[which(mit_Fem_percond_FC12Cutoofs$avg_log2FC < -0.32),])
mit_Male_percond_FC12Cutoofs <-bind_rows(mit_Male_percond_FC12Cutoofs[which(mit_Male_percond_FC12Cutoofs$avg_log2FC>0.25),],mit_Male_percond_FC12Cutoofs[which(mit_Male_percond_FC12Cutoofs$avg_log2FC < -0.32),])

dim(mit_percond_FC12Cutoofs)
dim(mit_Fem_percond_FC12Cutoofs)
dim(mit_Male_percond_FC12Cutoofs)

length(unique(mit_percond_FC12Cutoofs$Gene)) 
length(unique(mit_Fem_percond_FC12Cutoofs$Gene)) 
length(unique(mit_Male_percond_FC12Cutoofs$Gene)) 
DEGs_permissive_stats <- as.data.frame(mit_percond_FC12Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 
DEGs_permissive_statsfem <- as.data.frame(mit_Fem_percond_FC12Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 
DEGs_permissive_statsmale <- as.data.frame(mit_Male_percond_FC12Cutoofs %>% dplyr::group_by(cellType_res05) %>% dplyr::summarise(n())) 

resDegs_permissive <- list( LA=DEGs_permissive_stats, fem=DEGs_permissive_statsfem,male=DEGs_permissive_statsmale)
  write.xlsx(resDegs_permissive, file = paste0(wd,f_results,f_Mit,f_tables,"Mit_DEGs_permissive_Stats2.xlsx"));


##########
load(file = paste0(wd,f_results,f_Rdata,f_allClusters,"scRNAseq_Seurat_object.RData"));#scRNASeq_seurat

PercentageFeatureSet_addFeatures.function <- function(scRNASeq_seurat,features,namePattern4col,nameOutput){
    nameColPattern <- paste0("percent.",namePattern4col);
    scRNASeq_seurat[[nameColPattern]] <- PercentageFeatureSet(scRNASeq_seurat, features = features);
  assign(paste0("scRNASeq_seurat_",nameOutput),scRNASeq_seurat,.GlobalEnv);
};
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat,unique(mit_percond_Cutoofs$Gene),"cond_MtDysfunction","mit")
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat_mit,unique(mit_Fem_percond_Cutoofs$Gene),"cond_Fem_MtDysfunction","mit")
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat_mit,unique(mit_Male_percond_Cutoofs$Gene),"cond_Male_MtDysfunction","mit")




scRNASeq_seurat_mit_permissive <-scRNASeq_seurat_mit
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat_mit_permissive,unique(mit_percond_FC12Cutoofs$Gene),"permissivecond_MT","mit_permissive")
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat_mit_permissive,unique(mit_Fem_percond_FC12Cutoofs$Gene),"permissivecond_Fem_MT","mit_permissive")
PercentageFeatureSet_addFeatures.function(scRNASeq_seurat_mit_permissive,unique(mit_Male_percond_FC12Cutoofs$Gene),"permissivecond_Male_MT","mit_permissive")



# calculate the percentage of all the counts belonging to a subset of the possible features for each cell.
#  This is useful when trying to compute the percentage of transcripts that map to mitochondrial genes for example. 
FeaturePlot_FinalClustering.function <- function(findClustersList,width,height,nameOutput,date, metadataCondition,featureName,folderName){
  library(RColorBrewer);

  dir.create(file.path(wd, f_results,f_foundClusters,f_plots,f_featurePlot,folderName), showWarnings = F);
  DefaultAssay(findClustersList) <- "RNA";
  #
  if (metadataCondition=="Sample") {
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition) & scale_color_viridis_c() 
  
  pdf(file=paste0(wd,f_results,f_foundClusters,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".pdf"),onefile=TRUE,width = width, height = height)
  print(plot);
  dev.off();


  } else if (metadataCondition=="NULL"){
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE) & scale_color_viridis_c() & 
      theme(legend.text=element_text(size=15),legend.spacing.y = unit(2, 'cm'),legend.key.height = unit(2.3, "cm"),
        axis.title = element_text(size=22,face="bold"),axis.text.x = element_text(size=18, color="black"),axis.text.y = element_text(size=18, color="black"), 
        axis.line = element_line(size=1.5))


  png(file=paste0(wd,f_results,f_foundClusters,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".png"), height = height, width = width, units = "in", res = 1200);
  print(plot);
  dev.off();


  } else{
    plot <- FeaturePlot(findClustersList, features = featureName, raster = FALSE,split.by = metadataCondition) && scale_color_viridis_c() & 
    theme(legend.text=element_text(size=15),legend.spacing.y = unit(2, 'cm'),legend.key.height = unit(2.3, "cm"),
            axis.title = element_text(size=22,face="bold"),axis.text.x = element_text(size=18, color="black"),axis.text.y = element_text(size=18, color="black"), 
            axis.line = element_line(size=1.5))
;
  png(file=paste0(wd,f_results,f_foundClusters,f_plots,folderName,featureName,"_umap_",nameOutput,
    "_",metadataCondition,"_",date,".png"), height = height, width = width, units = "in", res = 1200);
  print(plot);
  dev.off();

  };
  #
};

f_annotCluster <- "MTD_annotated_Cluster_cellType/"
date <-"2022";
  

#DEGs: DEA_thresholds adj p-value< 0.05 -0.32<log2(FC)> 0.25
FeaturePlot_FinalClustering.function(scRNASeq_seurat_mit_permissive,7,6,"permissivecond_Mt",date, "NULL","percent.permissivecond_MT",f_annotCluster);
FeaturePlot_FinalClustering.function(scRNASeq_seurat_mit_permissive,7,6,"permissivecond_Fem_Mt",date, "NULL","percent.permissivecond_Fem_MT",f_annotCluster);
FeaturePlot_FinalClustering.function(scRNASeq_seurat_mit_permissive,7,6,"permissivecond_Male_Mt",date, "NULL","percent.permissivecond_Male_MT",f_annotCluster);


######## table with mean percent of genes linked to MT dysfunction altered per cluster per condition  or cond and sex  #############

#DEGs: DEA_thresholds adj p-value< 0.05 -0.32<log2(FC)> 0.25 
Idents(scRNASeq_seurat_mit_permissive) <- "cluster.percond_sex"
identsDf <- data.frame(cell=names(Idents(scRNASeq_seurat_mit_permissive)), cellType_res05=as.character(Idents(scRNASeq_seurat_mit_permissive)))
mitDisfLA<-data.frame(cell=names(scRNASeq_seurat_mit_permissive$percent.permissivecond_MT), mit=as.numeric(scRNASeq_seurat_mit_permissive$percent.permissivecond_MT))
mitDisfLA<- left_join(mitDisfcond,identsDf,by="cell")
mitDisfLA<-separate(mitDisfcond,cellType_res05,c("cellType_res05" ,"percond","Sex" ),sep="__");

permissiveDEGs_condition_Group_MT_perc <- as.data.frame(mitDisfLA  %>% dplyr::group_by(cond,cellType_res05 ) %>% dplyr::summarise(quantile095_premissive_MT = quantile(mit, 0.95)))
permissiveDEGs_condition_Group_MT_Sex_perc <- as.data.frame(mitDisfLA  %>% dplyr::group_by(cond,Sex,cellType_res05 ) %>% dplyr::summarise(quantile095_premissive_MT_sex = quantile(mit, 0.95)))
permissiveDEGs_condition_Group_MT_Sex_perc$cond_Sex <- paste0(permissiveDEGs_condition_Group_MT_Sex_perc$cond, "__", permissiveDEGs_condition_Group_MT_Sex_perc$Sex)
permissiveDEGs_condition_Group_MT_perc <- spread(permissiveDEGs_condition_Group_MT_perc, cond, quantile095_premissive_MT)
permissiveDEGs_condition_Group_MT_Sex_perc<- permissiveDEGs_condition_Group_MT_Sex_perc[,c(5,3,4)]
permissiveDEGs_condition_Group_MT_Sex_perc <- spread(permissiveDEGs_condition_Group_MT_Sex_perc, cond_Sex, quantile095_premissive_MT_sex)
scRNASeq_seurat_mit_permissive <- unique(scRNASeq_seurat_mit_permissive)
#

Idents(scRNASeq_seurat_mit) <- "cluster.percond_sex"
identsDf2 <- data.frame(cell=names(Idents(scRNASeq_seurat_mit)), cellType_res05=as.character(Idents(scRNASeq_seurat_mit)))
mitDisfLA2<-data.frame(cell=names(scRNASeq_seurat_mit$percent.cond_MtDysfunction), mit=as.numeric(scRNASeq_seurat_mit$percent.cond_MtDysfunction))
mitDisfLA2<- left_join(mitDisfLA2,identsDf2,by="cell")
mitDisfLA2<-separate(mitDisfLA2,cellType_res05,c("cellType_res05" ,"percond","Sex" ),sep="__");


Idents(scRNASeq_seurat_mit) <- "cluster.percond_sex"
identsDf3 <- data.frame(cell=names(scRNASeq_seurat_mit$cluster.percond_sex), cellType_res05=as.character(scRNASeq_seurat_mit$cluster.percond_sex))
mitDisfLA3<-data.frame(cell=names(scRNASeq_seurat_mit$percent.cond_MtDysfunction), mit=as.numeric(scRNASeq_seurat_mit$percent.cond_MtDysfunction))
mitDisfLA3<- left_join(mitDisfLA3,identsDf3,by="cell")
mitDisfLA3<-separate(mitDisfLA3,cellType_res05,c("cellType_res05" ,"percond","Sex" ),sep="__");
scRNASeq_seurat_mit <- unique(scRNASeq_seurat_mit)


