##### LM Iceland Crop Samples Statistical Analysis ##########

library(ChemoSpecUtils)
library(ChemoSpec)
library(R.utils)
library(lattice)
library(ggplot2)
library(gridExtra)
library(tidyverse)
library(doParallel)
library(randomForest)
library(caret)
library(prospectr)
library(vip)
library(dplyr)
library(e1071)
library(mosaic)


rm(list=ls())

setwd("") # Assign working directory


##### ChemoSpec PCA ------------------------------------------------------------

#For color-blind safe palette see  ggplot of scores (line 63 ff) 

data_species <- matrix2SpectraObject(gr.crit =c("E.nigrum","Vaccinium", "B.nana_C","B.nana_IFR","B.pubescens","S.phylicifolia","S.herbacea","D.octopetala","ind"), gr.cols = c("darkolivegreen3","red","black","cornflowerblue","gold", "darkblue","chocolate","chocolate4","purple"),
                                     freq.unit = "wavenumber",
                                     int.unit = "absorbance",
                                     descrip = "LM Iceland", in.file = "LM_Iceland_crop_samples_FTIR_data_ChemoSpec.csv",
                                     out.file = "data", chk = TRUE, sep=";" , dec= ",")


data_clean <- removeFreq(data_species, rem.freq = data_species$freq > 1800 & data_species$freq < 2500)

data_clean <- removeFreq(data_species, rem.freq = data_species$freq > 1800)

pca <- r_pcaSpectra(data_clean, choice= "noscale")

plotScores(data_clean, pca, pcs = c(1,2))
plotLoadings(data, pca, loads = c(1,2))


#####Extract PCA scores to csv 

write.csv(pca$x, "scores.csv")
write.csv(data_clean$names, "ID.csv")


#Extract PCA loadings to csv 
write.csv(pca$rotation, "loadings_LM.csv")
write.csv(data_clean$freq, "wavenumbers.csv")



##### PCA Scores plot for plant taxa and plant parts ---------------------------

# Prepared PCA scores data (merged scores, ID and sample infromatione) is saved in "scores_LM_full.csv"
#Prepared PCA loadings data (merged loadings and wavenumbers for PC1 and PC2) is saved in loadings_PC1.csv" and loadings_PC2.csv"


LM <- read.csv("scores_LM_full.csv", header=TRUE, sep=";", dec=".")
LM$species <- factor(LM$species, levels=c('Vaccinium', 'E.nigrum', "B.nana", "B.pubescens","D.octopetala", "S.herbacea", "S.phylicifolia"))
LM$part <- factor(LM$organ, levels=c('B', 'C', "IFR", "L", "Sb"))

p1 <- ggplot(LM, aes(x = PC1, y = PC2, colour = species, shape = organ)) +
  geom_point(size = 4) +
  labs(x = "\nPC1 (70 %)", y = "PC2 (13 %)\n") +
  scale_colour_viridis_d(name = NULL, option = "H") +  
  scale_shape_manual(
    name = NULL,
    labels = c('Berries', 'Catkins', "Infructescence","Leaves", "Stems with buds"),
    values = c(1, 0, 4, 3, 2)
  ) +
  guides(
    colour = guide_legend(label.theme = element_text(face = "italic", size = 23)),
    shape = guide_legend(label.theme = element_text(face = "plain", size = 23))
  ) +
  theme(
    legend.position = "right",
    legend.text = element_text(size= 23),
    axis.line = element_line(color = "black"),
    panel.background = element_rect(fill = "white"),
    axis.text = element_text(size = 20),
    axis.title = element_text(size = 23),
    plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches")
  )

p1

ggsave(
  "LM_scores_plot_full.tiff",
  plot = p1,
  device = "tiff",
  path = "",                  # Assign path
  scale = 1,
  width = 4500,
  height = 4000,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE)



##### PCA Loadings plots #######

#Loadings Plot - Fingerprint

PC1 <- read.csv("loadings_LM_fingerprint_PC1.csv", header=TRUE, sep=";", dec=".")
PC2 <- read.csv("loadings_LM_fingerprint_PC2.csv", header=TRUE, sep=";", dec=".")

P1 <- ggplot(PC1, aes(x = wavenumber, y = PC1)) +
  geom_line()+
  labs(title= NULL,x=" corresponding wavenumber (1/cm)", y="PC1 (70 %)")+
  scale_x_continuous(trans = scales::reverse_trans(),
                     breaks = c(3800,3500, 3200,  2900,  2600,   2300,  2000, 1700, 1400,1100,800,500),
                     limits=c(3800, 500))+
  theme(panel.background= element_rect("white"),
        axis.line = element_line(color="black", size = 0.3),
        axis.title.y=element_text(size=17), axis.text = element_text(size=17),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))


P1

P2 <- ggplot(PC2, aes(x = wavenumber, y = PC2)) +
  geom_line()+
  labs(title= NULL,x=" corresponding wavenumber (1/cm)", y="PC2 (13 %)")+
  scale_x_continuous(trans = scales::reverse_trans(),
                     breaks = c(3800,3500, 3200,  2900,  2600,   2300,  2000, 1700, 1400,1100,800,500),
                     limits=c(3800, 500))+
  theme(panel.background= element_rect("white"),
        axis.line = element_line(color="black", size = 0.3),
        axis.title.y=element_text(size=17), axis.text = element_text(size=17),
        plot.margin = unit(c(0.2, 0.2, 0.2, 0.2), "inches"))


P2

PC <- grid.arrange(P1, P2, nrow=2, ncol=1)

ggsave("loadings_fingerprint_PC1.tiff",
       plot = P1,
       device = "tiff",
       path = "",                  # Assign path       
       scale = 1,
       width = 3500,
       height = 2500,
       units = c("px"),
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,
       create.dir = FALSE)


ggsave("loadings_fingerprint_PC2.tiff",
       plot = P2,
       device = "tiff",
       path = "",                  # Assign path       
       scale = 1,
       width = 3500,
       height = 2500,
       units = c("px"),
       dpi = 300,
       limitsize = TRUE,
       bg = NULL,
       create.dir = FALSE)






##### Band aggregation function ------------------------------------------------
# If RandomForest on all (raw) wavenumbers is applied, skip to line 517
 
# Numerische Hilfsfunktionen

# Numerische Hilfsfunktionen
trapz <- function(x, y) {
  nx <- length(x)
  if (nx < 2) return(NA_real_)
  sum((x[-1] - x[-nx]) * (y[-1] + y[-nx])) / 2
}

lin_baseline <- function(x, y) {
  y0 <- y[1]; y1 <- y[length(y)]
  x0 <- x[1]; x1 <- x[length(x)]
  if (!is.finite(y0) || !is.finite(y1) || x1 == x0) return(y)
  m  <- (y1 - y0) / (x1 - x0)
  b  <- y0 - m * x0
  y - (m * x + b)
}

centroid <- function(x, y_pos) {
  S <- sum(y_pos)
  if (!is.finite(S) || S <= 0) return(NA_real_)
  sum(x * y_pos) / S
}

fwhm <- function(x, y) {
  ymax <- max(y, na.rm = TRUE)
  if (!is.finite(ymax) || ymax <= 0) return(NA_real_)
  half <- ymax / 2
  above <- which(y >= half)
  if (length(above) < 2) return(NA_real_)
  x[above[length(above)]] - x[above[1]]
}

# Band-aggregierte Merkmale für EIN Spektrum
band_features_one <- function(x_wn, y_abs, bands, calc_fwhm = FALSE,
                              min_points_area = 2,
                              min_points_height = 1,
                              adapt_expand = TRUE,
                              expand_step = 2,
                              max_expands = 3) {
  o <- order(x_wn)
  x <- x_wn[o]; y <- y_abs[o]
  out <- list()
  for (b in bands) {
    lo <- min(b$lower, b$upper); hi <- max(b$lower, b$upper)
    sel <- which(x >= lo & x <= hi)
    if (adapt_expand && length(sel) < min_points_height) {
      ex <- 0
      while (length(sel) < min_points_height && ex < max_expands) {
        ex <- ex + 1
        lo <- lo - expand_step
        hi <- hi + expand_step
        sel <- which(x >= lo & x <= hi)
      }
    }
    nm <- b$name
    
    A <- NA_real_; H <- NA_real_; Cen <- NA_real_; F <- NA_real_
    if (length(sel) >= min_points_height) {
      xx <- x[sel]; yy <- y[sel]
      yy_bc  <- lin_baseline(xx, yy)
      H <- max(yy_bc, na.rm = TRUE)
      
      if (length(sel) >= min_points_area) {
        yy_pos <- pmax(yy_bc, 0)
        A <- trapz(xx, yy_pos)
        S <- sum(yy_pos)
        if (is.finite(S) && S > 0) Cen <- centroid(xx, yy_pos)
        if (calc_fwhm) F <- fwhm(xx, yy_bc)
      } else {
        A <- 0
        Cen <- NA_real_
        if (calc_fwhm) F <- NA_real_
      }
    }
    vals <- c(A, H, Cen)
    names(vals) <- c(paste0(nm, "_area"),
                     paste0(nm, "_height"),
                     paste0(nm, "_centroid"))
    if (calc_fwhm) vals <- c(vals, setNames(F, paste0(nm, "_fwhm")))
    out[[nm]] <- vals
  }
  unlist(out)
}

# Definition der Bandfenster
make_default_bands <- function() {
  list(
    # OH
    list(name = "OH_3600_3500", lower = 3500, upper = 3700),
    list(name = "OH_3490_3100", lower = 3100, upper = 3490),
    # cis-Alkenes
    list(name = "CisAlkene_3020_3000", lower = 3000, upper = 3020),
    # CH2
    list(name = "CH2_asym_2930_2910", lower = 2910, upper = 2930),
    list(name = "CH2_sym_2860_2840",  lower = 2840, upper = 2860),
    # Carbonyl
    list(name = "TriGly_1745_1738",   lower = 1738, upper = 1745),
    list(name = "Cutin_1733_1730",    lower = 1730, upper = 1733),
    list(name = "Cutin_1711_1708",    lower = 1708, upper = 1711),
    
    # Proteine - Aromatics
    list(name = "AmideI_1670_1625", lower = 1625, upper = 1670),
    list(name = "AmideII_1580_1520",lower = 1500, upper = 1580),
    
    list(name = "Phenolics_1610_1600",lower = 1600, upper = 1610),
    list(name = "Lignin_1515_1500",   lower = 1500, upper = 1515),
    
    # Carbohydrates
    list(name = "Carb_1165_1156", lower = 1156, upper = 1165),
    list(name = "Carb_1155_1140", lower = 1140, upper = 1155),
    
    list(name = "Carb_1058_1046", lower = 1046, upper = 1058),
    list(name = "Carb_1035_1010", lower = 1010, upper = 1035),
    list(name = "Cellu_1000_990", lower = 990, upper = 1000),
    list(name = "Pectin_975_960", lower = 960, upper = 975),
    
    # Fructose
    list(name = "Fructose_820_815", lower = 815, upper = 820),
    list(name = "Fructose_779_773", lower = 773, upper = 779),
    
    # Oxalat/Tannine
    list(name = "Oxalate_782_780",   lower = 780, upper = 782),
    list(name = "Oxalate_1320_1310", lower = 1310, upper = 1320),
    list(name = "Oxalate_1620_1618", lower = 1618, upper = 1620),
    
    list(name = "Tannin_1288_1282",  lower = 1282, upper = 1288),
    list(name = "Tannin_763_758",    lower = 758,  upper = 763),
    list(name = "Tannin_1722_1716",  lower = 1716, upper = 1722)
  )
}

# Hauptfunktion NUR für aggregierte Band-Features
aggregate_bands_df <- function(df,
                               meta_cols = 1:9,
                               remove_region = c(2500, 1800),
                               bands = make_default_bands(),
                               calc_fwhm = FALSE) {
  stopifnot(all(meta_cols %in% seq_len(ncol(df))))
  spec_cols <- setdiff(seq_len(ncol(df)), meta_cols)
  if (length(spec_cols) < 3) stop("Not enough columns found")
  
  # Wellenzahlen aus Spaltennamen extrahieren
  coln <- colnames(df)[spec_cols]
  wn_vals <- suppressWarnings(as.numeric(coln))
  if (anyNA(wn_vals)) {
    wn_vals <- suppressWarnings(as.numeric(sub(".*?(\\d+\\.?\\d*).*", "\\1", coln)))
  }
  if (anyNA(wn_vals)) stop("Colnames must include wavenumbers ('3304' or 'wn_3304').")
  
  # Artefakt-Region entfernen
  keep <- rep(TRUE, length(wn_vals))
  if (!is.null(remove_region) && length(remove_region) == 2) {
    lo <- min(remove_region); hi <- max(remove_region)
    keep <- !(wn_vals <= hi & wn_vals >= lo)
  }
  wn_keep <- wn_vals[keep]
  if (length(wn_keep) < 3) stop("Not enough wavenumbers left after artefact regions were removed")
  
  # Spektralmatrix
  mat <- as.matrix(df[, spec_cols, drop = FALSE])[, keep, drop = FALSE]
  if (!is.numeric(mat)) stop("Spectral data is not numeric.")
  
  # Band-Features berechnen
  feat_list <- lapply(seq_len(nrow(mat)), function(i) {
    band_features_one(wn_keep, mat[i, ], bands, calc_fwhm)
  })
  lens <- vapply(feat_list, length, integer(1))
  if (length(unique(lens)) != 1) {
    stop("Band feature vectors have different lengths. Check band windows.")
  }
  feat_mat <- do.call(rbind, feat_list)
  feat_df  <- as.data.frame(feat_mat, check.names = TRUE)
  rownames(feat_df) <- rownames(df)
  
  # Ausgabe: Meta + Band-Features
  cbind(df[, meta_cols, drop = FALSE], feat_df)
}



##### Apply band aggregation on data -------------------------------------------
###### Data for plant parts ################

setwd("...")

data <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",")


# Apply band aggregation function on data
result_agg <- aggregate_bands_df(data,
                                 meta_cols    = 1:9,
                                 remove_region= c(2500, 1800),
                                 bands        = make_default_bands(),
                                 calc_fwhm    = FALSE
)

result <- result_agg

meta_cols    <- 1:9
target_col   <- "organ"  

cand_cols <- setdiff(seq_along(result), c(meta_cols, which(names(result) == target_col)))

#Decide which metric/s to be used in RF grepl("(_area|_height)$" 

feature_mask <- grepl("(_height)$", names(result))

feature_cols <- intersect(which(feature_mask), cand_cols)

#Remove constants and NA
is_constant <- vapply(result[ , feature_cols, drop=FALSE], function(v) {
  v_allna <- all(is.na(v))
  v_allna || (length(unique(v[!is.na(v)])) <= 1)
}, logical(1))
feature_cols <- feature_cols[!is_constant]

X_all <- as.matrix(result[ , feature_cols, drop = FALSE])
storage.mode(X_all) <- "numeric"

organ <- as.factor(result[[target_col]])
data_organ <- cbind.data.frame(organ = organ , as.data.frame(X_all), stringsAsFactors = FALSE)




###### Data for plant species and part ################

data <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",")

# Apply band aggregation function on data
result_agg <- aggregate_bands_df(data,
                                 meta_cols    = 1:9,
                                 remove_region= c(2500, 1800),
                                 bands        = make_default_bands(),
                                 calc_fwhm    = FALSE
)

result <- result_agg

meta_cols    <- 1:9
target_col   <- "species"  

cand_cols <- setdiff(seq_along(result), c(meta_cols, which(names(result) == target_col)))

#Decide which metric/s to be used in RF grepl("(_area|_height|_centroid)$" 

feature_mask <- grepl("(_height)$", names(result))

feature_cols <- intersect(which(feature_mask), cand_cols)

#Remove constants and NA
is_constant <- vapply(result[ , feature_cols, drop=FALSE], function(v) {
  v_allna <- all(is.na(v))
  v_allna || (length(unique(v[!is.na(v)])) <= 1)
}, logical(1))
feature_cols <- feature_cols[!is_constant]

X_all <- as.matrix(result[ , feature_cols, drop = FALSE])
storage.mode(X_all) <- "numeric"

species <- as.factor(result[[target_col]])
data_species <- cbind.data.frame(species = species , as.data.frame(X_all), stringsAsFactors = FALSE)


###### Random Forests -----------------------------------------------------

setwd("...")

# Using raw data for RandomForest on all wavenumbers ---------------------------
data <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",")

data <- data %>% select(-sample)
data <- data %>% select(-id)
data <- data %>% select(-species)
#data <- data %>% select(-organ)   # --> plant part
data <- data %>% select(-sex_age)
data <- data %>% select(-elevation)
data <- data %>% select(-year)
data <- data %>% select(-dry_mass)
data <- data %>% select(-biomass)


# Using Data after band aggregation ---------------------------------------------
data <- data_organ
#data <- data_species --> Choose variable


# Save classification variable as factor
data$organ <- as.factor(data$organ)

#Split data into training and test data
set.seed(7)
train_set <- data %>%
  group_by(organ) %>%
  slice_sample(prop = 0.8)  

test_set <- anti_join(data, train_set)  

table(train_set$organ)
table(test_set$organ)


set.seed(7)
cv_rf_grid <- 
  cv_rf_grid <- function(data, target_var,
                         mtry_vals = c(5,10, 20),
                         nodesize_vals = c(1, 5,10),
                         maxnodes_vals = c(50, 100),
                         ntree_vals = c(100, 300, 500),
                         cv_folds = 10)
{
  
  set.seed(7)
  
  results <- data.frame()
  
  
  folds <- createFolds(data[[target_var]], k = cv_folds, list = TRUE, returnTrain = FALSE)
  
  for (mtry in mtry_vals) {
    for (nodesize in nodesize_vals) {
      for (maxnodes in maxnodes_vals) {
        for (nt in ntree_vals) {
          
          accs <- c()  
          
          for (i in 1:cv_folds) {
            test_idx <- folds[[i]]
            test_data <- data[test_idx, ]
            train_data <- data[-test_idx, ]
            
            rf <- randomForest(
              formula = as.formula(paste(target_var, "~ .")),
              data = train_data,
              mtry = mtry,
              nodesize = nodesize,
              maxnodes = maxnodes,
              ntree = nt
            )
            
            preds <- predict(rf, test_data)
            acc <- mean(preds == test_data[[target_var]])
            accs <- c(accs, acc)
          }
          
          mean_acc <- mean(accs)
          
          results <- rbind(results, data.frame(
            mtry = mtry,
            nodesize = nodesize,
            maxnodes = maxnodes,
            ntree = nt,
            mean_accuracy = round(mean_acc, 4)
          ))
          
          cat("✔ mtry:", mtry, "| nodesize:", nodesize, "| maxnodes:", maxnodes, "|🌲 ntree:", nt, "| mean CV-accuracy:", round(mean_acc, 4), "\n")
        }
      }
    }
  }
  
  
  results <- arrange(results, desc(mean_accuracy))
  return(results)
}



# Save target variable as factor
train_set$organ <- as.factor(train_set$organ)
test_set$organ <- as.factor(test_set$organ)


# Hyperparameter tuning through cross-validation
cv_results <- cv_rf_grid(
  data = train_set,
  target_var = "organ",
  mtry_vals = c(5,10, 20),
  nodesize_vals = c(1, 5,10),
  maxnodes_vals = c(50, 100),
  ntree_vals = c(100, 300, 500),
  cv_folds = 10
)

head(cv_results)
stopCluster(cl)


###### RandomForest for organ (plant part) -------------------------------------
train_set$organ <- as.factor(train_set$organ)
test_set$organ <- as.factor(test_set$organ)

#Hyperparameters are set for RF on aggregated bands -> For hyperparameters on raw wavenumbers see Supplemental or run hyperparameter optimization (line 540 ff)
set.seed(7)
rfp_model <- randomForest(organ ~ ., 
                          data = train_set, ntree = 500, mtry = 5, maxnodes=50, nodesize=1,
                          importance = TRUE)

print(rfp_model)

###### RandomForest for species ------------------------------------------------
train_set$species <- as.factor(train_set$species)
test_set$species <- as.factor(test_set$species)

#Hyperparameters are set for RF on aggregated bands -> For hyperparameters on raw wavenumbers see Supplemental or run hyperparameter optimization (line 540 ff)
set.seed(7)
rftp_model <- randomForest(species ~ ., 
                           data = train_set, ntree = 100, mtry = 5, maxnodes=100, nodesize=1, importance =TRUE)

print(rftp_model)



##### RandomFroest Prediction --------------------------------------------------
rf_pred <- predict(rfp_model, test_set)
rf_cm <- confusionMatrix(rf_pred, test_set$organ)
print(rf_cm)


#Print accuracy
cat("RF Accuracy:", rf_cm$overall["Accuracy"], "\n")




#####Random Forest Variable Importance -----------------------------------------
n_runs <- 30
importance_results <- list()

for (i in 1:n_runs) {
  set.seed(7 + i)
  
  rfp_model <- randomForest(organ ~ ., 
                            data = train_set, ntree = 500, mtry = 5, maxnodes=50, nodesize=1,
                            importance = TRUE)
  
  
  imp <- importance(rfp_model)
  
  imp_df <- data.frame(
    Feature = rownames(imp),
    MeanDecreaseAccuracy = imp[, "MeanDecreaseAccuracy"],
    MeanDecreaseGini = imp[, "MeanDecreaseGini"],
    Run = i
  )
  
  importance_results[[i]] <- imp_df
}

# Combine results
importance_all <- bind_rows(importance_results)

# Reduce to the first 26 feature
top_features <- importance_all %>%
  group_by(Feature) %>%
  summarise(MeanAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE)) %>%
  arrange(desc(MeanAccuracy)) %>%
  slice_head(n = 26) %>%
  pull(Feature)

importance_filtered <- importance_all %>%
  filter(Feature %in% top_features)

# Export to csv
write.csv(importance_filtered, "species_variable_importance_runs.csv", row.names = FALSE)

# Print for control
print(head(importance_filtered, 10))






##### Food composition in crop contents ----------------------------------------
rm(list=ls())

setwd("...") # Assign working directory

df <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",", stringsAsFactors = FALSE)

##### Add Zeros for unselected plant parts in individual ptarmigan an years
all_combos <- expand_grid(
  id  = unique(df$id),
  species = unique(df$species)
)

df_complete <- all_combos %>%
  left_join(df, by = c("id", "species")) %>%
  mutate(biomass = if_else(is.na(biomass), 0, biomass))
# Years have to be adpated manually for new created cells of zero, the file LM_Iceland_food_comp_zeros.csv already holds the adapted list.

write.csv(df_complete, "LM_Iceland_zeros.csv", row.names = FALSE)


#Data + zeros
data <- read.csv("LM_Iceland_food_comp_zeros.csv", header=TRUE, sep=";", dec=",")

#Data without zeros
data <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",", stringsAsFactors = FALSE)



##### Subset data to plant parts
C <- data[data$sample %in% c("B.nana_C","B.pubescens_C"), ]
B <- data[data$sample %in% c("Vaccinium_B", "E.nigrum_B"), ]
Sb <- data[data$sample %in% c("S.phylicifolia_Sb","S.herbacea_Sb"), ]
IFR <- data[data$sample %in% c("B.nana_IFR"), ]
L <- data[data$sample %in% c("D.octopetala_L"), ]




##### Violin plots #####
catkins <- ggplot(C, aes(x= as.factor(year), y= biomass, fill=species))+
  geom_violin(position = position_dodge(width = 1),width=0.7,  drop=FALSE) +
  labs(title= "Catkins",y = "Biomass (%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("skyblue","cornflowerblue"), label =c("B. nana", "B. pubescens"))+ 
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 


catkins


berries <- ggplot(B, aes(x= as.factor(year), y=biomass, fill=species))+
  geom_violin(position = position_dodge(width = 1), width=0.7) +
  labs(title= "Berries",y = "Biomass (%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18,),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("grey", "red2"))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 

berries


sb <- ggplot(Sb, aes(x= as.factor(year), y=biomass, fill=species))+
  geom_violin(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs( title= "Stems and buds",y = "Biomass(%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("brown","chocolate1"))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 

sb


leaves <- ggplot(L, aes(x=as.factor(year), y=biomass, group=year,fill=species))+
  geom_violin(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs(title= "Leaves", x = "Year",
       y = "Biomass (%)")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  scale_fill_manual(name=NULL, values=c("green3"))+
  stat_summary(fun="median", geom="point") 

leaves


infructescence <- ggplot(IFR, aes(x=as.factor(year), y=biomass, group=year, fill=species))+
  geom_violin(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs(title= "Infructescence", x = "Year",
       y = "Biomass (%)")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  scale_fill_manual(name=NULL, values=c("orange"),label=c("B. nana"))+
  stat_summary(fun="median", geom="point") 

infructescence

p1 <- grid.arrange( leaves, infructescence, catkins, sb, berries, nrow=3, ncol=2)

ggsave(
  "boxplots_biomass_violin.tiff",
  plot = p1,
  device = "tiff",
  path = "C:/Users/schae/Desktop",
  scale = "" ,     # assign path  
  width = 4000,
  height = 4000,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE)



##### Boxplots ####

catkins <- ggplot(C, aes(x= as.factor(year), y= biomass, fill=species))+
  geom_boxplot(position = position_dodge(width = 1),width=0.7) +
  labs(title= "Catkins",y = "Biomass (%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("skyblue","cornflowerblue"), label =c("B. nana", "B. pubescens"))+ 
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 


catkins


berries <- ggplot(B, aes(x= as.factor(year), y=biomass, fill=species))+
  geom_boxplot(position = position_dodge(width = 1), width=0.7) +
  labs(title= "Berries",y = "Biomass (%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18,),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("grey", "red2"))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 

berries


sb <- ggplot(Sb, aes(x= as.factor(year), y=biomass, fill=species))+
  geom_boxplot(position = position_dodge(width = 1),width=0.7) +
  labs( title= "Stems and buds",y = "Biomass(%)", x= "Year")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_fill_manual(name=NULL, values=c("brown","chocolate1"))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  stat_summary(fun="median", geom="point") 

sb


leaves <- ggplot(L, aes(x=as.factor(year), y=biomass, group=year,fill=species))+
  geom_boxplot(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs(title= "Leaves", x = "Year",
       y = "Biomass (%)")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  scale_fill_manual(name=NULL, values=c("green3"))+
  stat_summary(fun="median", geom="point") 

leaves


infructescence <- ggplot(IFR, aes(x=as.factor(year), y=biomass, group=year, fill=species))+
  geom_boxplot(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs(title= "Infructescence", x = "Year",
       y = "Biomass (%)")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  scale_fill_manual(name=NULL, values=c("orange"),label=c("B. nana"))+
  stat_summary(fun="median", geom="point") 

infructescence

p1 <- grid.arrange( leaves, infructescence, catkins, sb, berries, nrow=3, ncol=2)


ggsave(
  "boxplots_biomass_zeros_separate.tiff",
  plot = p1,
  device = "tiff",
  path = "",  # assign path
  width = 4500,
  height = 4500,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE)
infructescence <- ggplot(IFR, aes(x=as.factor(year), y=biomass, group=year, fill=species))+
  geom_boxplot(position = position_dodge(width = 1),width=0.7, drop=FALSE) +
  labs(title= "Infructescence", x = "Year",
       y = "Biomass (%)")+
  theme(legend.text = element_text(size=18, face ="italic"), legend.position ="bottom", plot.title = element_text(size=18),axis.title = element_text(size=18), axis.text = element_text(size=16),
        plot.margin = unit(c(0.25, 0.25, 0.25, 0.25), "inches"),
        panel.background = element_rect(fill ="white"),
        panel.grid = element_line(colour= "grey", linewidth = 0.2))+
  scale_x_discrete(limits = as.character(2006:2014)) +
  scale_fill_manual(name=NULL, values=c("orange"),label=c("B. nana"))+
  stat_summary(fun="median", geom="point") 

infructescence

p1 <- grid.arrange( leaves, infructescence, catkins, sb, berries, nrow=3, ncol=2)


ggsave(
  "boxplots_biomass_zeros_separate.tiff",
  plot = p1,
  device = "tiff",
  path = "",  # assign path
  width = 4500,
  height = 4500,
  units = c("px"),
  dpi = 300,
  limitsize = TRUE,
  bg = NULL,
  create.dir = FALSE)



