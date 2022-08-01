library("pathfindR");library("readxl");library("Seurat");library("dplyr");library("cowplot");library("ggplot2");library("MAST");library(topGO);
library("openxlsx");library("grid");library("gridBase");library("gridExtra");library("patchwork");library("org.Hs.eg.db");
library("stringr");library("tidyr");library("stringr"); library(readxl);
options(future.globals.maxSize = 100000 * 1024^2) #Default is 500 * 1024 ^ 2 = 500 Mb) using the following code:


Sys.setenv(RSTUDIO_PANDOC="home/software/anaconda3/envs/r4-base/bin/pandoc")
# Sys.getenv("RSTUDIO_PANDOC")
#Before begining:
#rm(list =ls()) ## erasing all the enviroment variables
set.seed(22); # need to set it to get always the same random results and plots
#sessionInfo()

#wd:
#####################################################
wd <- "home/omicData/scRNASeq/"
setwd(wd);

f_results <- "/results/";
f_Rdata <- "RData/";
f_tables <- "tables/";
f_plots <- "plots/";
f_dfa_topGo <- "dfa_topGO/"
f_list <- "res_list_dfs/"
f_pathfindR <-"pathfindR/"
f_customPaths <-"enrichment_MTDpaths/"
#########################################################  ##  #####################################################
dir.create(file.path(wd, f_results), showWarnings = F);
dir.create(file.path(wd, f_results,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results,f_pathfindR), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_Rdata,f_list), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_tables), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_plots), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_dfa_topGo,f_Rdata,f_list), showWarnings = F);

dir.create(file.path(wd, f_results, f_pathfindR,f_customPaths), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_customPaths,f_Rdata), showWarnings = F);
dir.create(file.path(wd, f_results, f_pathfindR,f_customPaths,f_tables), showWarnings = F);

#########################################################  ##  #####################################################

#MX database
input_mitoXplorer_hsa<- as.data.frame(read_excel(paste0("home/mitoExplorerDB/input/05272022/humanInteractome/human_gene_function.xlsx"),sheet =1, col_names =T))  
input_mitoXplorer_hsa <- unique(input_mitoXplorer_hsa[,c("ENSG_ID","mito_process","gene_name","gene_function","ENSG_name")]);
colnames(input_mitoXplorer_hsa)[c(1:5)]<- c("ensembl_gene_id","Mitochondrial_process","full_gene_name","gene_function","Gene")
length(unique(input_mitoXplorer_hsa$full_gene_name)) # Nb unique genes associated with Mit act: 1229 form v2020 to version 2022 no change
input_mitoXplorer_hsa$Group <- input_mitoXplorer_hsa$Mitochondrial_process;

#DEGs input
dataInput<- as.data.frame(read_excel(paste0("home/omicData/scRNASeq/rData/input/DEGs/MAST_EGs.xlsx"),sheet =1, col_names =T)); 
dataName <-"cond1_vs_cond2";

########################################################################################################
########################################################################################################
# Function to calculate significance of custom pathways we need each pathway to be in a list 
# and the genes involved in a vector without header.

########################################################################################################
########################################################################################################

inputcustomPath_formatting_pathFindR.function <- function(customPathsDf,nameOutput,date,input_testDf,nameAnalysis){
	# Function to calculate significance of custom pathways we need each pathway to be in a list 
	# and the genes involved in a vector without header.	customPath.list<- list();
	############################################################################################
	dir.create(file.path(wd,f_results, f_customPaths,"MTDpaths/"), showWarnings = F);
	
	customPath.list<- list();
	res_clusters.List <- list();


	for (i in 1:length(unique(customPathsDf$Group))){
		genes_in_path <- unique(customPathsDf[which(customPathsDf$Group==unique(customPathsDf$Group)[i]),"Gene"]);

		customPath.list[[i]] <- genes_in_path;
        names(customPath.list)[i] <- paste0("Path",i);
        i <- i+1;
	};
	customPath_descriptions <- setNames(unique(customPathsDf$Group), names(customPath.list))
	save(customPath.list, customPath_descriptions, file=paste0(wd,f_results, f_pathfindR,f_customPaths,f_Rdata, nameOutput,"_",date,"_custom_pathwayFiles_4pathFindR.RData"))
	

	for ( j in 1:length(unique(dataInput$cluster))){
		print(paste0("cluster being analysed: c",unique(dataInput$cluster)[j] ));
		input_testDf.tmp <- dataInput[which(dataInput$cluster==unique(dataInput$cluster)[j]),]
		input_testDf.tmp <- input_testDf.tmp[!is.na(input_testDf.tmp$Gene),]
		input_testDf <- unique(input_testDf.tmp[,c("Gene","avg_log2FC","p_val_adj")]);
		cluster <- unique(input_testDf.tmp$cluster);
		cellType_res05 <- unique(input_testDf.tmp[which(input_testDf.tmp$cluster==cluster),"cellType_res05"]);
		dir.create(file.path(wd,f_results,f_pathfindR, f_customPaths,"MTDpaths/",cluster,"/"), showWarnings = F);
	
		custom_pathFindR_res <- run_pathfindR(input_testDf,
	                               gene_sets = "Custom",
	                               custom_genes = customPath.list,
	                               custom_descriptions = customPath_descriptions,plot_enrichment_chart =FALSE,visualize_enriched_terms=FALSE,
	                               max_gset_size = Inf, # DO NOT LIMIT GENE SET SIZE
	                               output_dir = paste0(wd,f_results, f_pathfindR,f_customPaths,"MTDpaths/",cluster,"/",nameAnalysis));

		custom_pathFindR_res$cluster <- rep(cluster, times=length(custom_pathFindR_res$ID))
		custom_pathFindR_res$cellType_res05 <- rep(cellType_res05, times=length(custom_pathFindR_res$ID));
		
		res_clusters.List[[j]] <- custom_pathFindR_res;
		names(res_clusters.List)[j] <- cluster;
		j<- j+1;
	};


	for (ii in 1:(length(res_clusters.List))){
		if(ii==1){
		  df <- res_clusters.List[[ii]];
		  sign_pathsDf=df;
		} else{
		  df.tmp <- res_clusters.List[[ii]];
		  sign_pathsDf<- bind_rows(sign_pathsDf,df.tmp);

		};
		ii <- ii+1;
	};

	openxlsx::write.xlsx(sign_pathsDf, file = paste0(wd,f_results,f_pathfindR,f_customPaths,f_tables,nameOutput,"_",nameAnalysis, "_MTDpaths_significance.xlsx"),rowNames =FALSE, colNames =TRUE);
	save(res_clusters.List,sign_pathsDf, file=paste0(wd,f_results,f_pathfindR,f_customPaths,f_tables,nameOutput,"_",nameAnalysis,"_",date,"_MTDpaths_significance.RData"));
};



#DEGs input
degs_input <-dataInput[which(dataInput$p_val_adjPlot<0.05),]
degs_input <-bind_rows(degs_input[which(degs_input$avg_log2FC>1),],degs_input[which(degs_input$avg_log2FC < -1),]);
inputcustomPath_formatting_pathFindR.function(input_mitoXplorer_hsa,dataName,"072122",degs_input,"DEGs");

#EGs imput by log2FC 0.25 and -0.32
degs_input2 <-dataInput[which(dataInput$p_val_adjPlot<0.05),]
degs_input2 <-bind_rows(degs_input2[which(degs_input2$avg_log2FC>0.25),],degs_input2[which(degs_input2$avg_log2FC < -0.32),]);
inputcustomPath_formatting_pathFindR.function(input_mitoXplorer_hsa,dataName,"072122",degs_input2,"permissiveDEGs");
