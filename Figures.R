library(dplyr)
library(workflowsets)
library(tidymodels)
library(data.table)
library(ggplot2)
library(finetune)
library(parsnip)
library(baguette)
library(DALEX)
library(fastshap)
library(vip)
library(NeuralNetTools)
library(stringr)
library(patchwork)

options(bitmapType='cairo')

IncludeFigHere <- "/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/Oct_20231108_1935/"

#Oct_20231031_1118/"
#Oct_20231030_1425
setwd(IncludeFigHere)


load("Prediction.RData")
load("stacks.RData")
load("ModelingWorkflows.RData")
load("Modeling2data3mods.RData")


###################
################### Fig 1

mae_stats <- collect_metrics(race_results, summarize = FALSE) %>%
  filter(.metric == "mae") %>%
  separate(wflow_id, into = c("Recipe", "Model_Type"), sep = "_", remove = F, extra = "merge") 


PerfMAE <- list()

for (i in unique(mae_stats$Model_Type)){
  
  full_model <- filter(mae_stats, 
                       Recipe == "AllRec2J", Model_Type == i) %>%
    dplyr::select(id, full_model= .estimate)
   
  replacements <- c("Imaging" = "x IDPs ", "Education" = "x CR ", "Clinical" = "x AgeSex ", "Genes" = "x SNPs ")
  
  
  rest <- filter(mae_stats, 
                 Recipe != "AllRec2J", Model_Type == i) %>%
    dplyr::select(id, .estimate, Recipe) %>% 
    mutate(NewRecipe = str_replace_all(Recipe, replacements)) %>%
    mutate(NewRecipe = sub("^x ", "", NewRecipe))

  
 rest2 <- left_join(rest, full_model, by = 'id') %>%
    mutate(loss = .estimate - full_model) %>%
    mutate(term =  NewRecipe) %>%
    nest_by(term) %>%
    mutate(stats = list(t.test(data$loss ~ 1, data = data))) %>%
    summarise(tidy(stats), .groups = "drop") 
    
   #label_colors <- ifelse(rest2$term %in% c(2, 3), "blue", "red") - impossible
    
    PerfMAE[[i]] <- rest2 %>%
   ggplot(aes(x= estimate, y = reorder(term , estimate))) +
    geom_point() + 
    geom_vline(xintercept = 0, lty = 2, col = "red") + 
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = .25) + 
    labs(x = "loss of MAE", y = NULL) + theme_bw() + ggtitle(paste0("Algorithm: ", i))
  
  #PerfMAE[[i]]  + theme(axis.text.x = element_text(color = label_colors))
  
  
}

#PerfMAE[["xgboost"]] <- NULL

pdf(paste0(IncludeFigHere, "/PerfromanceMAEFinalFig.pdf"), 10, 7)
print(wrap_plots(PerfMAE))
dev.off()

library(patchwork)
pdf(paste0(IncludeFigHere, "/PerfromanceMAEFinalFig.pdf"), 13, 7)
(PerfMAE[[1]] + PerfMAE[[2]] + PerfMAE[[3]])/(plot_spacer() + PerfMAE[[4]] + PerfMAE[[5]] + plot_spacer()  + plot_layout(widths = c(0.3, 2, 2, 1)))
dev.off()

###################
################### Fig 2
#############





prettify_engine_col <- function(data) {
  res <- data %>% mutate_at(vars(Model), ~sprintf('{%s}', Model))
}

factor_src <- function(x) {
  ordered(x, levels = c('vip_model', 'vip_shap', 'vip_permute', 'dalex'))
}

VarImp <- race_results$wflow_id %>% 
  setdiff(c("Education_nnet", "Education_glm", "Education_ranger", "Education_xgboost", "Education_lasso")) ####add all education related results here, why was lasso not here?

VarImp <- VarImp[!grepl("glm", VarImp)]


names(result) <- VarImp

vi_rnks2 <- result %>%
  bind_rows() %>%
  as.data.frame()

Threshold <- 10

#IncludeFigHere <- "/rds/projects/g/gkoutosg-variant-prediction/Laura/Magda/FinalPipeline_20211024_2033"

library(biomaRt)
ensembl_snp1=useMart("ENSEMBL_MART_SNP", dataset="hsapiens_snp")
ensembl_snp2=useMart("ENSEMBL_MART_ENSEMBL", dataset="hsapiens_gene_ensembl")

vi_rnks2 <- vi_rnks2 %>%
  mutate(refsnp_id = sub("_.*", "", var))


filtered_vector <- grep("^rs", vi_rnks2$refsnp_id, value = TRUE)  # Filter elements that start with "rs"

Map <- getBM(attributes=c("refsnp_source",'refsnp_id','chr_name','chrom_start','chrom_end',"consequence_type_tv", "clinical_significance","ensembl_gene_stable_id"), 
             filters = 'snp_filter', 
             values = filtered_vector,
             mart = ensembl_snp1)

#Add gene names!

Map2 <- getBM(attributes=c("ensembl_gene_id", "hgnc_symbol" ), 
              filters = 'ensembl_gene_id', 
              values = Map$ensembl_gene_stable_id,
              mart = ensembl_snp2) %>%
  dplyr::rename("ensembl_gene_stable_id" = "ensembl_gene_id")

Table <- left_join(Map, Map2) %>%
  dplyr::select(refsnp_id,hgnc_symbol ) 

modified_table <- Table %>%
  group_by(refsnp_id) %>%
  summarise(V3 = paste(unique(hgnc_symbol), collapse = "-"))

save(modified_table, file = "modified_table.RData")

load("modified_table.RData")

replacement_rules <- c(
  "WM_ICVF_ForcepsMinor" = "ForcepsMinor (ICVF)",
  "WM_OD_SLF_Right" = "SLF_R (ODI)",
  "WM_OD_SLF_Left" = "SLF_L (ODI)",
  "WM_ICVF_SLF_Left" = "SLF_L (ICVF)",
  "WM_ICVF_SLF_Right" = "SLF_R (ICVF)",
  "WM_OD_ForcepsMinor" = "ForcepsMinor (ODI)",
  "WM_OD_IFOF_Left" = "IFOF_L (ODI)",
  "WM_OD_IFOF_Right" = "IFOF_R (ODI)",
  "WM_ICVF_IFOF_Right" = "IFOF_R (ICVF)",
  "WM_ICVF_IFOF_Left" = "IFOF_L (ICVF)",
  "c7" = "rFPN-dDMN",
  "c35" = "vDMN-rDMN",
  "c33" = "rFPN-rDMN",
  "C21" = "lFPN-vDMV",
  "C16" = "dDMN-vDMN",
  "C15" = "lFPN-rFPN",
  "C111" = "lFPN-ECN",
  "C106" = "dDMN-ECN", 
  "c11" =  "dDMN-lFPN", 
  "c20" =   "rFPN-vDMN",
  "c29" = "dDMN-rDMN",
  "c34" = "lFPN-rDMN",
  "c110" = "rFPN-ECN",
  "c112" = "vDMN-ECN",
  "c114"= "rDMN-ECN",
  "e4e41" = "e4e4_1(APOE)", 
  "e3e31" = "e3e3_1(APOE)", 
  "ClinicalAge" = "Age", 
  "Clinical_Sex1" = "Sex",
  "ClinicalEducation" = "CR", 
  "e2e30" = "e2e3_0(APOE)",
  "Clinical_Sex0" = "Sex", 
  "e4e4" = "e4e4(APOE)", 
  "Clinical_EducationMeasurement" = "CR (Education)", 
  "Clinical_Sex" = "Sex", 
  "GeneralCog" = "GCAS", 
  "e2e3" = "e2e3(APOE)",
  "e3e3" = "e3e3(APOE)",
  "e3e4" = "e3e4(APOE)"
)

vi_rnks2b <- vi_rnks2 %>%
  left_join(., modified_table) %>%
  mutate(Final = ifelse(!is.na(V3), paste0(var,"(", V3,")"), var)) %>%
  mutate(Final = gsub("_X", "", Final)) %>%
  mutate(Final2 = str_replace_all(Final, regex(replacement_rules, ignore_case = TRUE))) 


for (j in unique(vi_rnks2b$Dataset)){
  
  viz <-
    vi_rnks2b %>% 
    filter(Dataset == j) %>%
    group_by(var) %>% 
    mutate(rnk_mean = rnk %>% mean(na.rm = TRUE)) %>% 
    ungroup() %>% 
    #mutate_at(vars(var), ~forcats::fct_reorder(., -rnk_mean)) %>% 
    #ungroup() %>% 
    #prettify_engine_col() %>% 
    #mutate_at(vars(src), ~ordered(., levels = c('vip_model', 'vip_shap', 'vip_permute', 'dalex'))) %>% 
    mutate(lab = sprintf('%2d (%s)', rnk, scales::percent(imp_abs_norm, accuracy = 1, width = 2, justify = 'right'))) %>% 
    filter(rnk < Threshold)  %>%
    ggplot() +
    aes(x = src, y = Final2) +
    geom_tile(aes(fill = rnk), alpha = 0.5, show.legend = F) +
    geom_text(aes(label = lab)) +
    scale_fill_viridis_c(direction = -1, option = "D", na.value = 'white') +
    theme_minimal(base_family = '') +
    facet_wrap(~Model) +
    theme(
      strip.text = element_text(size=15),
      axis.text.x=element_text(size=12),
      plot.title.position = 'plot',
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = 'bold'),
      plot.subtitle = ggtext::element_markdown(),
    ) +
    labs(x = NULL, y = NULL)
  
  pdf(paste0(IncludeFigHere,"/NewHeatmaps_",j,".pdf"), 15, 12)
  
  print(viz)
  
  dev.off()
  
  ggsave(paste0(IncludeFigHere,"/NewHeatmapSave_",j,".pdf"), plot = viz, device = "pdf")
  #ggsave(paste0(IncludeFigHere,"/HeatmapSave_",j,".png"), plot = viz, device = "png")
  
}


###################
################### Fig 3

SelectFeatures <- result %>%
  bind_rows() %>%
  filter(rnk < 3) 

SelectFeatures <- unique(SelectFeatures$var) %>% 
  setdiff(c("Clinical_EducationMeasurement", "Clinical_Sex_X0", "e2e3_X0","rs17046359_G_X2"))

#https://towardsdatascience.com/doing-and-reporting-your-first-mediation-analysis-in-r-2fe423b92171


WithCogAndSNPs4ModeAPOE <- read.csv("/rds/projects/c/chechlmy-brain-ageing-ukbiobank/Laura/WithCogAndSNPs4ModeAPOE.csv")

data <- WithCogAndSNPs4ModeAPOE %>% 
  dplyr::select(-c(FID, Age2)) %>%
  dplyr::mutate_at(vars(Clinical_Sex, rs80136977_A:e2e3), as.factor)

col_y <- "GeneralCog" #outcome
plotType <- "D" #plot format
y <- data[[col_y]]

.seed <- 11
set.seed(.seed)

splits <- initial_split(data, strata = col_y, prop = 5/6)
data_trn <- splits %>% rsample::training()
data_tst <- splits %>% rsample::testing()


AllRec <- recipe(GeneralCog ~ ., data = data) #not data_trn
PreparedPreProc <- AllRec %>% prep()
AllRecJ <- juice(PreparedPreProc) #training data - basis of the rest

AllRec2 <- AllRec %>%
  step_normalize(all_numeric()) %>% #had skip = TRUE : took it out
  step_dummy(all_nominal(),  one_hot = FALSE)  %>% #change all from true to false
  step_zv(all_predictors())
PreparedPreProc <- AllRec2 %>% prep()
AllRec2J <- juice(PreparedPreProc)  # All variables dummified


# rename(dv = GeneralCog,
#       iv = Clinical_EducationMeasurement)



#install.packages("psych")
library(psych)

# has a cool comorbidity function? 
#HAVE TO ZERO CENTRE: glbwarmc <-data.frame(scale(globalWarm,scale=FALSE))

#registerDoSEQ() 

mod.glb <- list()
Params <- list()

SelectFeatures <- setdiff(SelectFeatures, c("rs118186402_A_X2", "rs77433227_T_X2", "rs141892450_A_X2", "rs73353563_T_X2")) #inside calculation - Error in dimnames(x) <- dn : 
#length of 'dimnames' [2] not equal to array extent

ToMediate <- AllRec2J %>%
  dplyr::select(c(SelectFeatures), GeneralCog, Clinical_EducationMeasurement) #%>%


ToMediate2 <- data.frame(Features = names(ToMediate)) %>%
  mutate(refsnp_id = ifelse(grepl("^rs", Features),
                            sub("^(rs[^_]+)_.*$", "\\1", Features),
                            Features)) %>%
  left_join(., modified_table) %>%
  mutate(Final = ifelse(!is.na(V3), paste0(refsnp_id,"(", V3,")"), refsnp_id)) %>%
  mutate(Final = gsub("_X", "", Final)) %>%
  mutate(Final2 = str_replace_all(Final, regex(replacement_rules, ignore_case = TRUE))) #%>%
  #mutate(Selected = ifelse(Features %in% SelectFeatures, TRUE,FALSE ))
  
  
names(ToMediate) <- ToMediate2$Final2

#dplyr::filter(ToMediate2,Selected == TRUE)$Final2

m <-setdiff(names(ToMediate), c("GCAS","CR (Education)"))

for (i in setdiff(names(ToMediate), c("GCAS","CR (Education)", "rs798507(GNA12)"))) {
  
  print(i)
  
  
  ToMediate3 <- ToMediate %>% 
    #dplyr::select(GeneralCog, Clinical_EducationMeasurement, i)
   dplyr::select(GCAS, 'CR (Education)', i)
  
  #debugonce(psych::mediate)
  mod.glb[[i]] <-  mediate(y = "GCAS", x = "CR (Education)", m= eval(i), data = ToMediate3, plot = TRUE)
  
  Params[[i]] <- data.frame(Name = i, 
                            a = paste0(round(as.numeric(mod.glb[[i]][["a"]]),4)," (", round(mod.glb[[i]][["a.reg"]][["se"]][2],4), ")"),
                            b = paste0(round(as.numeric(mod.glb[[i]][["b"]]),4), " (", round(mod.glb[[i]][["b.reg"]][["se"]][2],4), ")"),
                            c = paste0(round(as.numeric(mod.glb[[i]][["c"]]),4), " (", round(mod.glb[[i]][["cprime.reg"]][["se"]][2],4), ")"),
                            #ab = as.numeric(round(mod.glb[[i]][["ab"]],5)), 
                            'ab boost' = paste0(round(as.numeric(mod.glb[[i]][["boot"]]$mean[1]),4), " (", round(mod.glb[[i]][["boot"]]$sd[1],4), ")"),
                            'CI 2.5' = paste0(round(as.numeric(mod.glb[[i]][["boot"]]$ci.ab[1]),4)),
                            'CI 97.5' = paste0(round(as.numeric(mod.glb[[i]][["boot"]]$ci.ab[2]),4))
                            # = paste0(as.numeric(mod.glb[[i]][["boot"]]$mean[1]), "(", mod.glb[[i]][["boot"]]$sd[1], ")"),
                            #  ab = as.numeric(mod.glb[[i]][["ab"]]), 
  )
    
    
    #data.frame(Mean = mod.glb[[i]][["boot"]]$mean, Sd = mod.glb[[i]][["boot"]]$sd, ci.low = mod.glb[[i]][["boot"]]$ci[1,], 
     #                       ci.high = mod.glb[[i]][["boot"]]$ci[2,])
  
  
}

save(Params, file = "Params.RData")

save(mod.glb, file = "mod.glb.RData") # Table with values better!

#load("mod.glb.RData")
#load("Params.RData")


s <- do.call(rbind, Params)

fwrite(s, "MediationDiagrams2.csv")



print("Mediation Figures")

pdf(paste0(IncludeFigHere, "/MediateDiagramAll.pdf"), 10, 5)

for (i in setdiff(names(ToMediate), c("GCAS","CR (Education)", "rs798507(GNA12)"))) {
  
  
print(mediate.diagram(mod.glb[[i]]))
  
}

dev.off()

#print("Mediation Figures Parameters")
#
#ParamsPlot <-   Params %>%
#  bind_rows(.,.id = 'Features') %>%
#  mutate(Type = case_when(startsWith(Features, "WM") | startsWith(Features, "c") ~ "IDP FPN", 
#                          startsWith(Features, "rs") | startsWith(Features, "e") ~ "SNPs", 
#                          Features %in% c("ClinicalAge","Clinical_Sex_X1")  ~ "Demographics", 
#                          Features == "Clinical_EducationMeasurement" ~ "CR",
#                          TRUE ~ "GCAS")) %>% 
#  mutate(refsnp_id = ifelse(grepl("^rs", Features),
#                                   sub("^(rs[^_]+)_.*$", "\\1", Features),
#                                   Features)) %>%
#  left_join(., modified_table) %>%
#  mutate(Final = ifelse(!is.na(V3), paste0(refsnp_id,"(", V3,")"), refsnp_id)) %>%
#  mutate(Final = gsub("_X", "", Final)) %>%
#  mutate(Final2 = str_replace_all(Final, regex(replacement_rules, ignore_case = TRUE)))
#  
#  
#ParamsPlot2 <- ParamsPlot %>% 
#  ggplot(aes(x=Final2, y = `Mean.1`,color = Type)) + #Before Mean
#  geom_hline(yintercept = 0, alpha = 0.2) +
#  geom_point() +
#  geom_errorbar(aes(ymin = `Mean.1` - `Sd.1`, ymax = `Mean.1` + `Sd.1`)) + theme_bw() + coord_flip()  +
#  labs(x = "Weight", y = "Features")
#
#pdf(paste0(IncludeFigHere, "/MediationPerf.pdf"), 10, 5)
#print(ParamsPlot2)
#dev.off()
#

################## Supp!

#Fig 4


Model1aPivot <- WithCogAndSNPs4ModeAPOE %>%
  dplyr::select(-Age2) %>%
  pivot_longer(-c(FID), names_to = "Features", values_to = "Value") %>%
  mutate(Type = case_when(startsWith(Features, "WM") | startsWith(Features, "c") ~ "FPC IDPs", 
                          startsWith(Features, "rs") | startsWith(Features, "e") ~ "SNPs", 
                          Features %in% c("ClinicalAge","Clinical_Sex")  ~ "Demographics", 
                          Features == "Clinical_EducationMeasurement" ~ "CR (Education)",
                          TRUE ~ "Cognitive Tests"))

library(RColorBrewer)
myColors <- brewer.pal(5,"Set1")
names(myColors) <- levels(as.factor(Model1aPivot$Type))
colScale <- scale_fill_manual(name = "Type",values = myColors)

Networks <- Model1aPivot %>%
  filter(Type == "FPC IDPs" | Features %in% c("ClinicalAge", "GeneralCog", "Clinical_EducationMeasurement") ) # %>% 
# mutate(Type2 = str_replace(Type, "Networks", "Brain Imaging"))
# mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%
# mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%
#   mutate(Type = ifelse(Type == "Networks", "Brain Imaging", .)) %>%  

#  General cognitive ability score (GCAS) 


NetworksPlot <- Networks %>%
  mutate(Features = str_replace_all(Features, regex(replacement_rules, ignore_case = TRUE))) %>%
  ggplot(aes(Value)) +
  facet_wrap(~Features, scales = "free") +
  geom_histogram(aes(fill = Type)) + 
  theme_bw() +  colScale


pdf(paste0(IncludeFigHere, "/ContinousPLot.pdf"), 16,  9)
print(NetworksPlot)
dev.off()

Genes <- Model1aPivot %>%
  filter( Features %in% c( "Clinical_Sex") | Type == "SNPs") %>%
  mutate( Value = as.character(Value)) %>%
  count( Features, Value ,sort = TRUE) %>%
  group_by(Features) %>%
  mutate(Sum = sum(n), 
         Freq = round(n/Sum, 2 )) %>%
  mutate(Type = case_when(Features %in% c("ClinicalAge","Clinical_Sex") ~ "Demographics", TRUE ~ "SNPs")) %>%
  ungroup() %>% 
  mutate(refsnp_id = ifelse(grepl("^rs", Features),
                            sub("^(rs[^_]+)_.*$", "\\1", Features),
                            Features)) %>%
  left_join(., modified_table) %>%
  mutate(Final = ifelse(!is.na(V3), paste0(refsnp_id,"(", V3,")"), refsnp_id)) %>%
  mutate(Final = gsub("_X", "", Final)) %>%
  mutate(Final2 = str_replace_all(Final, regex(replacement_rules, ignore_case = TRUE))) %>%
  ggplot(aes(Freq, Value )) +
  geom_col(aes(fill = Type)) +
  facet_wrap(~Final2, scales = "free_y") + 
  theme_bw() +   colScale


pdf(paste0(IncludeFigHere, "/AllPlots.pdf"), 16,  15) #13
NetworksPlot/Genes
dev.off()

############# MISSING VALUES

options(bitmapType = "cairo")

load("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/Magda/SelectedParticipants3.RData")

library(naniar)
library(finalfit)

ll <- data.frame(is.na(SelectedParticipants3))  

NewNames <- str_replace_all(names(ll), regex(replacement_rules, ignore_case = TRUE)) %>% 
  gsub("Cog_", "Cognitive Test: ", .)  %>% 
  gsub("Clinical_EducationMeasurement", "CR (Education)", .) %>% 
  gsub("Clinical_", "", .) %>% 
  gsub("Clinical", "", .) %>% 
  gsub("DateBirth", "Date of Birth", .) %>% 
  gsub("Date_Assessment", "Date of Assessment", .)

names(ll) <- NewNames 
  
cols <- sapply(ll, is.logical)
ll[, cols] <- lapply(ll[, cols], as.numeric)

Miss1 <- UpSetR::upset(ll,
                       nsets = 40, number.angles = 20, point.size = 3.5, line.size = 2,
                       mainbar.y.label = "Missing Values", sets.x.label = "Total Number Missing Values",
                       text.scale = c(2.3, 2.3, 2, 2, 2, 1.75), order.by = "freq", sets.bar.color = "red3"
)

pdf("SFig1_MissingVal1.pdf", 25, 30)
print(Miss1)
dev.off()



Data <- SelectedParticipants3 %>%
  drop_na(Clinical_EducationMeasurement) %>%
  dplyr::select(-c("Clinical_Date_Assessment", "ClinicalDateBirth"))

ll <- data.frame(is.na(Data)) # %>%
cols <- sapply(ll, is.logical)
ll[, cols] <- lapply(ll[, cols], as.numeric)


NewNames <- str_replace_all(names(ll), regex(replacement_rules, ignore_case = TRUE)) %>% 
  gsub("Cog_", "Cognitive Test: ", .)  %>% 
  gsub("Clinical_EducationMeasurement", "CR (Education)", .) %>% 
  gsub("Clinical_", "", .) %>% 
  gsub("Clinical", "", .) %>% 
  gsub("DateBirth", "Date of Birth", .) %>% 
  gsub("Date_Assessment", "Date of Assessment", .)

names(ll) <- NewNames 



Miss1 <- UpSetR::upset(ll,
                       nsets = 40, number.angles = 20, point.size = 3.5, line.size = 2,
                       mainbar.y.label = "Missing Values", sets.x.label = "Total Number Missing Values",
                       text.scale = c(2.3, 2.3, 2, 2, 2, 1.75), order.by = "freq", sets.bar.color = "red3"
)

pdf("SFig2_MissingVal.pdf", 34, 23)
print(Miss1)
dev.off()


pdf("Missing2.pdf", 20, 10)
print(vis_miss(Data, warn_large_data = FALSE) + theme(axis.text.x =  element_text(angle = 90)))
dev.off()


#pdf("MissingPlotPatterns.pdf",25,25)
# Missing data patterns
  
  
  NewNames <- str_replace_all(names(SelectedParticipants3), regex(replacement_rules, ignore_case = TRUE)) %>% 
  gsub("Cog_", "Cognitive Test: ", .)  %>% 
  gsub("Clinical_EducationMeasurement", "CR (Education)", .) %>% 
  gsub("Clinical_", "", .) %>% 
  gsub("Clinical", "", .) %>% 
  gsub("DateBirth", "Date of Birth", .) %>% 
  gsub("Date_Assessment", "Date of Assessment", .)

names(SelectedParticipants3) <- NewNames 

library(mice)
  
PlotMissingMagda <- function(x){
  
  #x <- SelectedParticipants3
  
  R <- is.na(x)
  nmis <- colSums(R)
  # sort columnwise
  R <- matrix(R[, order(nmis)], dim(x))
  pat <- apply(R, 1, function(x) paste(as.numeric(x), collapse = ""))
  # sort rowwise
  sortR <- matrix(R[order(pat), ], dim(x))
  if (nrow(x) == 1) {
    mpat <- is.na(x)
  } else {
    mpat <- sortR[!duplicated(sortR), ]
  }
  
  # update row and column margins
  if (all(!is.na(x))) {
    cat(" /\\     /\\\n{  `---'  }\n{  O   O  }\n==>  V <==")
    cat("  No need for mice. This data set is completely observed.\n")
    cat(" \\  \\|/  /\n  `-----'\n\n")
    mpat <- t(as.matrix(mpat, byrow = TRUE))
    rownames(mpat) <- table(pat)
  } else {
    if (is.null(dim(mpat))) {
      mpat <- t(as.matrix(mpat))
    }
    rownames(mpat) <- table(pat)
  }
  
  r <- cbind(abs(mpat - 1), rowSums(mpat))
  r <- rbind(r, c(nmis[order(nmis)], sum(nmis)))
  
  
  op <- par(mar = rep(0, 4))
  on.exit(par(op))
  plot.new()
  
  R <- t(as.matrix(r[1:nrow(r) - 1, 1:ncol(r) - 1]))
  R <- r[1:nrow(r) - 1, 1:ncol(r) - 1]
  
  adj <- c(0, 0.5)
  srt <- 90
  length_of_longest_colname <- max(nchar(colnames(r))) / 2.6
  
  plot.window(
    xlim = c(-1, ncol(R) + 1), #c(-1, 2),
    ylim = c(-1, nrow(R) + length_of_longest_colname), #c(-1, 100), #,
    asp = 1
  )
  
  TextSize <- 0.3
  
  M <- cbind(c(row(R)), c(col(R))) - 1
  shade <- ifelse(R[nrow(R):1, ], mdc(1), mdc(2))
  rect(M[, 2], M[, 1], M[, 2] + 1, M[, 1] + 1, col = shade)
  for (i in 1:ncol(R)) {
    text(i - .5, nrow(R) + .3, colnames(r)[i], adj = adj, srt = srt, cex = TextSize )
    text(i - .5, -2.5, nmis[order(nmis)][i], adj = c(0.5,0.0), srt = 90,  cex = TextSize) #.3
  }
  for (i in 1:nrow(R)) {
    text(ncol(R) + .3, i - .5, r[(nrow(r) - 1):1, ncol(r)][i], adj = 0,  cex = TextSize)
    text(-.3, i - .5, rownames(r)[(nrow(r) - 1):1][i], adj = 1,  cex = TextSize)
  }
  text(ncol(R) + .3, -.3, r[nrow(r), ncol(r)],  cex = TextSize)
  
  return(r)
  
}



pdf("PlotMissing.pdf", 4,11)

PlotMissingMagda(SelectedParticipants3)

dev.off()



