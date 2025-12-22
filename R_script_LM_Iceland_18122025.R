
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


##### ChemoSpec PCA #########################################################################################
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

# Prepared PCA scores data (merged scores, ID and sample infromatione) is saved in "scores_LM_full.csv"
#Preapred PCA loadings data (merged loadings and wavenumbers for PC1 and PC2) is saved in loadings_PC1.csv" and loadings_PC2.csv"

##### PCA Scores plot for plant taxa and plant parts #######

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








###### Random Forest #######################################################################################
rm(list=ls())

setwd("") # Assign working directory

data <- read.csv("LM_Iceland_crop_samples_FTIR_data.csv", header=TRUE, sep=";", dec=",")

sample	id	species	organ	sex_age	elevation	year	dry_mass	biomass



# Remove all columns except for the classification variable
data <- data %>% select(-sample)
data <- data %>% select(-id)
data <- data %>% select(-species)
#data <- data %>% select(-organ)   # --> plant part
data <- data %>% select(-sex_age)
data <- data %>% select(-elevation)
data <- data %>% select(-year)
data <- data %>% select(-dry_mass)
data <- data %>% select(-biomass)


# Save classification variable as factor
data$organ <- as.factor(data$organ)

#Split data into training and test data
set.seed(123)
train_set <- data %>%
  group_by(organ) %>%
  slice_sample(prop = 0.8)  

test_set <- anti_join(data, train_set)  

table(train_set$organ)
table(test_set$organ)


# Function for hyperparameter tuning
num_cores <- detectCores()
cl <- makeCluster(num_cores - 1)
registerDoParallel(cl)

set.seed(123)
cv_rf_grid <- function(data, target_var,
                       mtry_vals = c(10, 20),
                       nodesize_vals = c(1, 5),
                       maxnodes_vals = c(50, 100),
                       ntree_vals = c(50, 100, 200, 300, 500),
                       cv_folds = 10)
{
  
  set.seed(123)
  
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
  mtry_vals = c(10, 20, 30),
  nodesize_vals = c(1, 5, 10),
  maxnodes_vals = c(50, 100),
  ntree_vals  = c(50, 100, 200, 300, 500),
  cv_folds = 10
)

head(cv_results)
stopCluster(cl)



###### RandomForest for organ (plant part) ###########
set.seed(123)
rfp_model <- randomForest(organ ~ ., 
                         data = train_set, ntree = 100, mtry = 20, maxnodes=100, nodesize=1,
                         importance = TRUE)

print(rfp_model)

###### RandomForest for species ###########
set.seed(123)
rftp_model <- randomForest(species ~ ., 
                         data = train_set, ntree = 300, mtry = 30, maxnodes=100, nodesize=1, importance =TRUE)

print(rftp_model)


##### RandomFroest Predictions #############
rf_pred <- predict(rf_model, test_set)
rf_cm <- confusionMatrix(rf_pred, test_set$organ)
print(rf_cm)


#Print accuracy
cat("RF Accuracy:", rf_cm$overall["Accuracy"], "\n")




#####Random Forest Variable Importance for RFtp  ####
n_runs <- 30
importance_results <- list()

for (i in 1:n_runs) {
  set.seed(123 + i)
  
  rftp_model <- randomForest(species ~ ., 
                             data = train_set, ntree = 300, mtry = 30, maxnodes=100, nodesize=1, importance =TRUE)
  
  imp <- importance(model)
  
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

# Reduce to the first 30 feature
top_features <- importance_all %>%
  group_by(Feature) %>%
  summarise(MeanAccuracy = mean(MeanDecreaseAccuracy, na.rm = TRUE)) %>%
  arrange(desc(MeanAccuracy)) %>%
  slice_head(n = 30) %>%
  pull(Feature)

importance_filtered <- importance_all %>%
  filter(Feature %in% top_features)

# Export to csv
write.csv(importance_filtered, "species_variable_importance_runs.csv", row.names = FALSE)

# Print for control
print(head(importance_filtered, 10))






##### Food composition in crop contents #####
rm(list=ls())

setwd("C:/Users/schae/Desktop/Finale Dokumente Publikation Island/Data") # Assign working directory

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




