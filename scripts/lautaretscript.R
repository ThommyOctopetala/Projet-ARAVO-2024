#library####
library(googledrive)
library(factoextra)
library(FactoMineR)
drive<-drive_find()
View(drive)
RELEVES<-drive_download(file = "2022_RELEVES.txt",overwrite = TRUE)
aravo<-read.csv(aravo$local_path)


###

source("veget_tools_share.r")

# 2022_RELEVES.txt contains the sites x species table
# 2022_SITES.txt contains metadata on sites - correspondence with new site labelling
# 2022_SPECIES.txt contains species full names

ARAP   <- read.table(RELEVES$local_path,quote="",header=TRUE,sep="\t",row.names=1)
ARAP   <- ARAP[,!is.na(apply(ARAP,2,sum))]
dim(ARAP)
ARAP   <- ARAP[apply(ARAP,1,sum)!=0,]
dim(ARAP)
REL_TOT   <- t(ARAP)
SPL_TOT   <- colnames(REL_TOT)

FLO(REL_TOT,1)  # the function generates REL_TOT.xx 
# S in % is the threshold to discard rare species

(dim(REL_TOT.db)) # remains 97 rel x  79 species with S = 5%
(dim(REL_TOT.db)) # remains 97 rel x 129 species with S = 1%

SPL_TOTf  <- SPL_TOT[-REL_TOT.i[[3]]]  
REL_TOTf  <- REL_TOT[-REL_TOT.i[[4]],]
REL_TOT.b <- REL_TOT.b[-REL_TOT.i[[4]],]

# Partitionning usinf function pam
nclu  <- 10             # number of tested clusters (from 2 to 10)
TRY   <- REL_TOT.da

PARTa <- part(TRY,nclu,typedist="nonbinary",method="jaccard")
CLUST <- PARTa[[1]]

# Writing output : beware  " " becomes .
CURDAY <- gsub("-","",Sys.Date())
write.csv(CLUST, paste0(CURDAY,"_Dom_Flo_Clust.csv"))

for (i in 1:nclu){
  # save the dominant species per cluster
  write.csv(PARTa[[3]][[i]], paste0(CURDAY,"_Dom_Flo_Part_",i,".csv")) 
  
  # save the list of indicator species
  write.csv(PARTa[[4]][[i]], paste0(CURDAY,"_Indic_Spec_",i,".csv"))
  
  # save the assignation to cluster
  write.csv(PARTa[[5]][[i]], paste0(CURDAY,"_SignAssoc_Spec_",i,".csv")) 
}

# Extract the number of releves per partition / cluster
mypart <- list()
for (i in 1:nclu){
  mypart[[i]] <- table(PARTa[[1]][,i])
}
(mypart)

# Search for specific dominant species
SEARCHpart <- function(Dom_Flo_Part,mynclu,species){
  RES <- rep(NA,mynclu);names(RES)<-1:mynclu
  for (i in 1:mynclu){
    # when species are written genus.epithet use 
    # tmp <- unlist(lapply(strsplit(Dom_Flo_Part[[mynclu]][,i],"\\."),function(x) paste(x[1],x[2],sep=" ")))
    # otherwise use
    tmp <- Dom_Flo_Part[[mynclu]][,i]
    pos <- which(tmp==species)
    if (length(pos)>0) RES[i]<-pos
  }
  return(RES)
}

tryclu <- 8    # example with nclu = 8
mypart[[tryclu]]
SEARCHpart(PARTa[[3]],tryclu,"Alch pent")
SEARCHpart(PARTa[[3]],tryclu,"Care foet")
SEARCHpart(PARTa[[3]],tryclu,"Gnap supi")
SEARCHpart(PARTa[[3]],tryclu,"Sali herb")
SEARCHpart(PARTa[[3]],tryclu,"Sali retu")
SEARCHpart(PARTa[[3]],tryclu,"Care myos")
SEARCHpart(PARTa[[3]],tryclu,"Alch colo")  

# Species richness
RIC <- apply(REL_TOT.db,1,sum)
COMM <- c("CM",	"MS",	"CF",	"AR",	"FV",	"SH",	"PS",	"CS")
tmp <- as.factor(CLUST[,8])
levels(tmp) <- COMM

ORD  <- c(1,5,6,2,3,7,8,4)
graphics.off(); windows(10,10); boxplot(RIC~CLUST[,8],names=COMM[ORD], xlab="Cluster")

# Afetr a bit of reordering
data      <- data.frame(richness=RIC, comm=tmp)
new_order <- with(data, reorder(comm , RIC, na.rm=T))
graphics.off(); windows(10,10); boxplot(RIC~comm,data=data,xlab="Cluster")
graphics.off(); windows(10,10); boxplot(RIC~new_order,data=data,xlab="Cluster")
####PCA####

library(readxl)
aravo2001 <- read_excel("~/Desktop/M2/Lautaret/ARAVO_RESURVEYxlsx.xlsx", 
                        sheet = "2001_GRID.SPEC", col_types = c("text", 
                                                                "text", "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "numeric", "numeric", "numeric", 
                                                                "skip"))
aravo2022 <- read_excel("~/Desktop/M2/Lautaret/ARAVO_RESURVEYxlsx.xlsx", 
                                 sheet = "2022_GRID.SPEC", col_types = c("numeric", 
                                                                         "numeric", "skip", "skip", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "text", "numeric", "text", "text", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric", "text", "numeric", "numeric", 
                                                                         "text", "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "skip"))

aravo2022 <- aravo2022 %>%
  mutate(
    lig = as.numeric(lig),
    col = as.numeric(col),
    # Vérifier les grands nombres et les convertir en valeurs décimales
    lig = ifelse(lig > 40000, as.numeric(as.Date("1847-08-30") + lig+3.5), lig),
    col = ifelse(col > 40000, as.numeric(as.Date("1847-08-30") + col+3.5), col)
  )

aravo2001 <- aravo2001 %>%
  mutate(ID = paste0(...1, "_", col)) %>%
  select(-`...1`, -col)  # Supprimez les colonnes d'origine
aravo2001<- aravo2001 %>%
  select(ID, everything())
aravo2022 <- aravo2022 %>%
  mutate(ID = paste0(lig, "_", col)) %>%
  select(-`lig`, -col)  # Supprimez les colonnes d'origine
aravo2022<- aravo2022 %>%
  select(ID, everything())

aravo2022 <- aravo2022 %>%
  rename_with(~ gsub(" ", ".", .x))
aravo2001 <- aravo2001 %>%
  rename_with(~ gsub(" ", ".", .x))
aravo2001 <- aravo2001 %>% select(-ID,-Alch.colo)

listesp2022<-colnames(aravo2022)
listesp2001<-colnames(aravo2001)
aravo2022 <- aravo2022 %>% select(-ID,-Alch.colo,-Nardus.stri,-Pilosella.officinarum)


aravo_traits_filtered2001 <- aravo$traits %>%
  filter(rownames(aravo$traits) %in% listesp2001)
length(afcL.aravo2001$cw)  
nrow(aravo_traits_filtered2001)   
afcL.aravo2001 <- dudi.coa(aravo2001, scannf = FALSE)
scatter(afcL.aravo2001, xax = 1, yax = 2, 
        clab.row = 0.8, clab.col = 0.8, 
        main = "Analyse des correspondances pour aravo2001")
acpQ.aravo2001 <- dudi.pca(aravo_traits_filtered2001, row.w = afcL.aravo2001$cw,
                       scannf = FALSE)
scatter(acpQ.aravo2001, xax = 1, yax = 2, 
        clab.row = 0.8, clab.col = 0.8, 
        main = "ACP des traits avec les poids de l'AFC")
#2022###
aravo2022 <- aravo2022 %>%
  mutate(across(everything(), ~ as.numeric(.)))

# Remplacer les valeurs manquantes (NA) par zéro
aravo2022[is.na(aravo2022)] <- 0
listesp2022<-colnames(aravo2022)
aravo_traits_filtered2022 <- aravo$traits %>%
  filter(rownames(aravo$traits) %in% listesp2022)
afcL.aravo <- dudi.coa(aravo2022, scannf = FALSE)
length(afcL.aravo$cw)  
nrow(aravo_traits_filtered2022)   
afcL.aravo <- dudi.coa(aravo2022, scannf = FALSE)
scatter(afcL.aravo, xax = 1, yax = 2, 
        clab.row = 0.8, clab.col = 0.8, 
        main = "Analyse des correspondances pour aravo2001")
acpQ.aravo <- dudi.pca(aravo_traits_filtered2022, row.w = afcL.aravo$cw,
                       scannf = FALSE)
scatter(acpQ.aravo, xax = 1, yax = 2, 
        clab.row = 0.8, clab.col = 0.8, 
        main = "ACP des traits avec les poids de l'AFC")
####
data(aravo)
View(aravo$traits)
View(aravo$spe)
View(aravo$spe.names)
aravo$env
afcL.aravo <- dudi.coa(aravo$spe, scannf = FALSE)
acpR.aravo <- dudi.hillsmith(aravo$env, row.w = afcL.aravo$lw,
                             scannf = FALSE)
acpQ.aravo <- dudi.pca(aravo$traits, row.w = afcL.aravo$cw,
                       scannf = FALSE)
rlq.aravo <- rlq(acpR.aravo, afcL.aravo, acpQ.aravo,
                 scannf = FALSE)

plot(rlq.aravo)
par(mfrow = c(1, 3))
s.arrow(rlq.aravo$l1)
s.arrow(rlq.aravo$c1)
s.label(rlq.aravo$lQ, boxes = FALSE)
nrepet <- 49999
four.comb.aravo <- fourthcorner(aravo$env, aravo$spe,
                                aravo$traits, modeltype = 6, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)
####


 # Nombre de lignes dans aravo$traits
acpQ.aravo2001 <- PCA(aravo_traits_filtered2001, row.w = afcL.aravo2001$cw, graph = FALSE)

# Visualiser les individus (espèces) en coloriant chaque point
fviz_pca_biplot(acpQ.aravo2001, 
                geom.ind = "point",             # Représenter les individus par des points
                pointshape = 21, 
                pointsize = 3,
                col.ind = "cos2",               # Colorer les points en fonction de cos² (qualité de représentation)
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Dégradé de couleurs
                col.var = "contrib",            # Colorer les flèches des variables par contribution
                repel = TRUE,                   # Évite le chevauchement des étiquettes
                arrows = TRUE,                  # Ajoute des flèches pour les variables
                title = "Biplot ACP des traits et espèces")





acpQ.aravo2022 <- PCA(aravo_traits_filtered2022, row.w = afcL.aravo$cw, graph = FALSE)

# Visualiser les individus (espèces) en coloriant chaque point
fviz_pca_biplot(acpQ.aravo2022, 
                geom.ind = "point",             # Représenter les individus par des points
                pointshape = 21, 
                pointsize = 3,
                col.ind = "cos2",               # Colorer les points en fonction de cos² (qualité de représentation)
                gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"), # Dégradé de couleurs
                col.var = "contrib",            # Colorer les flèches des variables par contribution
                repel = TRUE,                   # Évite le chevauchement des étiquettes
                arrows = TRUE,                  # Ajoute des flèches pour les variables
                title = "Biplot ACP des traits et espèces")

#####
library(dplyr)
community_data2001 <- data.frame(
  species = c("Agro.rupe", "Alch.pent", "Care.foet", "Kobr.myos", "Care.rosa", 
              "Fest.quad", "Fest.viol", "Geum.mont", "Omal.supi", "Plan.alpi","Sali.herb","Sali.retu","Trif.alpi","Poa.supi"),
  community = c("AR", "CF", "CF", "CM", "AR", "CM","FV", "MS","PS","SH","SH", "CM", "MS","PS")  # Associez les communautés respectives ici
)

library(dplyr)
library(tibble)  # Assure que le package tibble est chargé pour utiliser rownames_to_column()
library(tidyr)
# Ajouter les noms des lignes comme une colonne "species" et joindre les données de communauté
aravo2001_long <- aravo2001 %>%
  pivot_longer(cols = everything(), names_to = "species", values_to = "presence")

# Joindre avec les données de communauté
aravo2001_long <- aravo2001_long %>%
  left_join(community_data2001, by = "species")

aravo_traits_filtered <- aravo_traits_filtered %>%
  rownames_to_column(var = "species")

# Joindre les données de traits avec aravo2001_long en fonction de la colonne "species"
aravo_combined <- aravo2001_long %>%
  left_join(aravo_traits_filtered, by = "species")

# Afficher le résultat
print(aravo_combined)



traits_data <- aravo_combined %>%
  select(Height, Spread, Angle, Area, Thick, SLA, N_mass, Seed)

# Réaliser le PCA avec FactoMineR
pca_result <- PCA(traits_data, scale.unit = TRUE, graph = FALSE)

# Créer le biplot avec les individus et les variables
aravo_combined$community <- as.factor(aravo_combined$community)

# Créer le biplot sans les ellipses
fviz_pca_ind(pca_result,
             geom.ind = "point",
             col.ind = aravo_combined$community, # Couleur selon la communauté
             palette = "jco",
             addEllipses = FALSE,               # Désactiver les ellipses
             legend.title = "Communities",
             repel = TRUE)                      # Évite le chevauchement des étiquettes

# Visualiser les variables (traits) en tant que flèches
fviz_pca_var(pca_result,
             col.var = "contrib",               # Colorer les flèches par contribution
             gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)    




library(ggrepel)

library(ggsci) # pour scale_color_jco si besoin

# Extraire les données numériques pour le PCA
traits_data <- aravo_combined %>%
  select(Height, Spread, Angle, Area, Thick, SLA, N_mass, Seed)

# Réaliser le PCA avec FactoMineR
pca_result <- PCA(traits_data, scale.unit = TRUE, graph = FALSE)

# Extraire les scores des individus et des variables pour ggplot2
individuals <- as.data.frame(pca_result$ind$coord) %>%
  mutate(community = aravo_combined$community) %>%
  mutate(species = aravo_combined$species)

variables <- as.data.frame(pca_result$var$coord) %>%
  rownames_to_column(var = "trait")

# Créer le biplot avec ggplot2
ggplot() +
  # Couche pour les individus (espèces)
  geom_point(data = individuals, aes(x = Dim.1, y = Dim.2, color = community), size = 2) +
  scale_color_jco() + # Utiliser la palette jco de ggsci pour avoir assez de couleurs
  labs(color = "Community") +
  # Ajouter les étiquettes des individus
  
  # Couche pour les variables (traits) avec des flèches
  geom_segment(data = variables, aes(x = 0, y = 0, xend = Dim.1, yend = Dim.2),
               arrow = arrow(length = unit(0.2, "cm")), color = "black") +
  # Ajouter les étiquettes des traits
  geom_text(data = variables, aes(x = Dim.1, y = Dim.2, label = trait), color = "red", size = 4) +
  # Personnalisation des axes et du titre
  labs(title = "Biplot ACP des traits et communautés",
       x = "PC1", y = "PC2") +
  theme_minimal()


#Kriegage Humidité####
library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Desktop/M2/Lautaret/Spatialaravo/Merged_DataFrame.csv")
data <- subset(data, Horizon == "H1")
# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$hum.pond <- as.numeric(gsub(",", ".", data$hum.pond))
data$x_l93 <- as.numeric(data$x_l93)
data$y_l93 <- as.numeric(data$y_l93)

# Supprimer les lignes avec des valeurs manquantes
data <- na.omit(data)
# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points")
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(hum.pond ~ 1, data)
plot(variogram_model, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(250, "Exp", 40, 100))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for humidity")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(hum.pond ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of Soil Moisture", fill = "Humidity") +
  theme_minimal()

library(raster)
khumidit<-krige_result
krige_raster <- raster(krige_result)

writeRaster(krige_raster, filename = "Humidité.tif", format = "GTiff", overwrite = TRUE)
#Kriegage MO####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Desktop/M2/Lautaret/Spatialaravo/Merged_DataFrame.csv")
data <- subset(data, Horizon == "H1")
# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$mo <- as.numeric(gsub(",", ".", data$mo))
data$x_l93 <- as.numeric(data$x_l93)
data$y_l93 <- as.numeric(data$y_l93)

# Supprimer les lignes avec des valeurs manquantes
data <- na.omit(data)
# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points", pch = 20, col = "blue")

# Ajouter les coordonnées x et y sur le graphique
plot(coords[,1], coords[,2], 
     xlab = "X Coordinate", 
     ylab = "Y Coordinate", 
     main = "Sampled Data Points", 
     pch = 5, col = "blue", asp = 1)
data_df <- as.data.frame(data)
coords <- coordinates(data)
data_df$x <- coords[,1]  # Ajouter la colonne des coordonnées x
data_df$y <- coords[,2]  # Ajouter la colonne des coordonnées y

# Créer un graphique avec ggplot
ggplot(data_df, aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +  # Tracer les points
  labs(title = "Sampled Data Points", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +  # Thème pour un rendu propre
  coord_fixed() 
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(mo ~ 1, data)

plot(variogram_model, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(120, "Sph", 20, 100))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Organic matter ")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(mo ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of  mo", fill = "mo") +
  theme_minimal()

library(raster)
kmo<-krige_result
krige_raster <- raster(krige_result)

writeRaster(krige_raster, filename = "krige_mo.tif", format = "GTiff", overwrite = TRUE)

















#Kriegage Caca####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Desktop/M2/Lautaret/Spatialaravo/Donn_es_Fusionn_es.csv")
data <- data %>% 
  filter(data$Dry.mass.poop..g.m.. != 0)
data <- data %>%
  group_by(ID) %>%
  summarize(
    Field1 = first(Field1),
    Field2 = first(Field2),
    Field3 = first(Field3),
    Field4 = first(Field4),
    ID = first(ID),
    secteur = first(secteur),
    corresp = first(corresp),
    Date = first(Date),
    avg_Dry_mass = mean(Dry.mass.poop..g.m.., na.rm = TRUE)
  )

# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$Dry.mass.poop..g.m.. <- as.numeric(gsub(",", ".", data$avg_Dry_mass))
data$x_l93 <- as.numeric(data$Field2)
data$y_l93 <- as.numeric(data$Field3)
data$z_l93 <- as.numeric(data$Field4)
# Supprimer les lignes avec des valeurs manquantes
data <- na.omit(data)
# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points", pch = 20, col = "blue")

# Ajouter les coordonnées x et y sur le graphique
plot(coords[,1], coords[,2], 
     xlab = "X Coordinate", 
     ylab = "Y Coordinate", 
     main = "Sampled Data Points", 
     pch = 5, col = "blue", asp = 1)
data_df <- as.data.frame(data)
coords <- coordinates(data)
data_df$x <- coords[,1]  # Ajouter la colonne des coordonnées x
data_df$y <- coords[,2]  # Ajouter la colonne des coordonnées y

# Créer un graphique avec ggplot
ggplot(data_df, aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +  # Tracer les points
  labs(title = "Sampled Data Points", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +  # Thème pour un rendu propre
  coord_fixed() 
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(Dry.mass.poop..g.m.. ~ 1, data)

plot(variogram_model, main = "Empirical Variogram")
variogram_model<-variogram_model[-1, ]
variogram_model<-variogram_model[-2, ]
# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(1.095, "Gau", 230, 0.1))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Caca ")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(Dry.mass.poop..g.m.. ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of  caca de mouton", fill = "caca") +
  theme_minimal()

library(raster)
kdry<-krige_result
krige_raster <- raster(krige_result)
setwd("~/Desktop/M2/Lautaret/Spatialaravo")
writeRaster(krige_raster, filename = "krige_poop.tif", format = "GTiff", overwrite = TRUE)





#Kriegage Caca with Z ####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Desktop/M2/Lautaret/Spatialaravo/Donn_es_Fusionn_es.csv")
data <- data %>% 
  filter(data$Dry.mass.poop..g.m.. != 0)
data <- data %>%
  group_by(ID) %>%
  summarize(
    Field1 = first(Field1),
    Field2 = first(Field2),
    Field3 = first(Field3),
    Field4 = first(Field4),
    ID = first(ID),
    secteur = first(secteur),
    corresp = first(corresp),
    Date = first(Date),
    avg_Dry_mass = mean(Dry.mass.poop..g.m.., na.rm = TRUE)
  )

# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$Dry.mass.poop..g.m.. <- as.numeric(gsub(",", ".", data$avg_Dry_mass))
data$x_l93 <- as.numeric(data$Field2)
data$y_l93 <- as.numeric(data$Field3)
data$z_l93 <- as.numeric(data$Field4)
# Supprimer les lignes avec des valeurs manquantes
data <- na.omit(data)
# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points", pch = 20, col = "blue")

# Ajouter les coordonnées x et y sur le graphique
plot(coords[,1], coords[,2], 
     xlab = "X Coordinate", 
     ylab = "Y Coordinate", 
     main = "Sampled Data Points", 
     pch = 5, col = "blue", asp = 1)
data_df <- as.data.frame(data)
coords <- coordinates(data)
data_df$x <- coords[,1]  # Ajouter la colonne des coordonnées x
data_df$y <- coords[,2]  # Ajouter la colonne des coordonnées y

# Créer un graphique avec ggplot
ggplot(data_df, aes(x = x, y = y)) +
  geom_point(color = "blue", size = 3) +  # Tracer les points
  labs(title = "Sampled Data Points", x = "X Coordinate", y = "Y Coordinate") +
  theme_minimal() +  # Thème pour un rendu propre
  coord_fixed() 
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(Dry.mass.poop..g.m.. ~ 1+y_l93, data)

plot(variogram_model, main = "Empirical Variogram")
variogram_model<-variogram_model[-1, ]
variogram_model<-variogram_model[-2, ]
# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(19.095, "Sph", 25, 0.1))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Caca ")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(Dry.mass.poop..g.m.. ~ data$y_l93, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of  caca de mouton", fill = "caca") +
  theme_minimal()

library(raster)

krige_raster <- raster(krige_result)
setwd("~/Desktop/M2/Lautaret/Spatialaravo")
writeRaster(krige_raster, filename = "krige_poop.tif", format = "GTiff", overwrite = TRUE)



data %>%
  filter(Field4 >= 2020L & Field4 <= 2680L) %>%
  ggplot() +
  aes(x = avg_Dry_mass, y = Field4) +
  geom_point(colour = "#112446") +
  theme_minimal()
data$mo<-as.numeric(data$mo)
data$paf<-as.numeric(data$paf)
summary(lm(data=data,mo~paf))


#Kriegage matiere sol MnO SURFACE de la merde il converge pas ce fils de putes ####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/Sol_roche/output data/fichier avec coord/XRF_grille_placette.csv")
data <- data[-43, ]
data <- subset(data, Type.de.sol == "Surface")
head(data$Type.de.sol)

# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$hum.pond <- as.numeric(gsub(",", ".", data$MnO))
data <- data[data$MnO < quantile(data$MnO, 0.95), ]
data$x_l93 <- as.numeric(data$x_l93)
data$y_l93 <- as.numeric(data$y_l93)

# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points")
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(MnO ~ 1, data)
plot(variogram_model, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(psill = 0.0005, model = "Sph", range = 550, nugget = 0.0009))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Oxyde")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(MnO ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of MnO", fill = "MnO") +
  theme_minimal()

library(raster)

krige_raster <- raster(krige_result)

writeRaster(krige_raster, filename = "krige_output.tif", format = "GTiff", overwrite = TRUE)

data <- data[data$MnO < quantile(data$MnO, 0.95), ]

dd<-krige_result@data["var1.pred"]
min(dd)
max(dd)


#Kriegage matiere sol MnO SURFACE de la merde il converge pas ce fils de putes ####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/Sol_roche/output data/fichier avec coord/XRF_grille_placette.csv")
data <- data[-43, ]
data <- subset(data, Type.de.sol == "Surface")
head(data$Type.de.sol)

# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$hum.pond <- as.numeric(gsub(",", ".", data$MnO))
data <- data[data$MnO < quantile(data$MnO, 0.95), ]
data$x_l93 <- as.numeric(data$x_l93)
data$y_l93 <- as.numeric(data$y_l93)

# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points")
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(hum.pond ~ 1, data)
plot(variogram_model, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(psill = 0.0005, model = "Sph", range = 550, nugget = 0.0009))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Oxyde")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(hum.pond ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of MnO", fill = "MnO") +
  theme_minimal()

library(raster)

krige_raster <- raster(krige_result)
kmnO<-krige_result
writeRaster(krige_raster, filename = "krige_output.tif", format = "GTiff", overwrite = TRUE)

data <- data[data$MnO < quantile(data$MnO, 0.95), ]

dd<-krige_result@data["var1.pred"]
min(dd)
max(dd)
#transformation data 
library(magick)
library(parallel)

# Charger les bibliothèques nécessaires
library(magick)
library(parallel)

# Chemin vers les images et où sauvegarder les images annotées
image_dir <- "/Volumes/Untitled1/ARAVOPHTO"
annotated_dir <- sprintf("%s/annotated", image_dir)

# Créer le répertoire annoté s'il n'existe pas
if (!dir.exists(annotated_dir)) {
  dir.create(annotated_dir, recursive = TRUE)
}

# Lire les chemins des images
image_paths <- list.files(path = image_dir, full.names = TRUE)

# Vérification pour éviter de continuer sans images
if (length(image_paths) == 0) {
  stop("Aucune image trouvée dans le répertoire spécifié")
}

# Définir le nombre de cœurs à utiliser
num_cores <- detectCores() - 1  # Garder un cœur libre

# Fonction pour convertir et annoter les images
convert_and_annotate <- function(image_path) {
  # Vérifier si le fichier est au format .heic
  if (tolower(tools::file_ext(image_path)) == "heic") {
    # Lire et convertir l'image HEIC en JPG
    img <- magick::image_read(image_path)
    img <- magick::image_convert(img, format = "jpg")
    # Remplacer l'extension pour le fichier de sortie
    output_path <- sprintf("%s/%s.jpg", annotated_dir, tools::file_path_sans_ext(basename(image_path)))
  } else {
    # Lire l'image si ce n'est pas un fichier HEIC
    img <- magick::image_read(image_path)
    output_path <- sprintf("%s/%s", annotated_dir, basename(image_path))
  }
  
  # Ajouter une annotation à l'image
  site_name <- substr(basename(image_path), 1, nchar(basename(image_path)) - 4)  # Nom du site sans l'extension
  annotated_img <- magick::image_annotate(img, text = site_name, size = 100, color = "red", location = "+30+30")
  
  # Sauvegarder l'image annotée
  magick::image_write(annotated_img, path = output_path)
  return(invisible())
}

# Charger et annoter les images en parallèle
images_annotated <- mclapply(image_paths, convert_and_annotate, mc.cores = num_cores)

cat("Traitement terminé pour toutes les images.\n")




image_dir <- "/Volumes/Untitled1/ARAVOPHTO"
annotated_dir <- sprintf("%s/annotated", image_dir)

# Créer le répertoire annoté s'il n'existe pas
if (!dir.exists(annotated_dir)) {
  dir.create(annotated_dir, recursive = TRUE)
}

# Lire les chemins des images dans les répertoires source et annoté
image_paths <- list.files(path = image_dir, full.names = TRUE)
annotated_paths <- list.files(path = annotated_dir, full.names = TRUE)

# Exclure les images déjà traitées
annotated_basenames <- basename(annotated_paths)
to_process <- image_paths[!basename(image_paths) %in% annotated_basenames]

# Vérification pour éviter de continuer sans images à traiter
if (length(to_process) == 0) {
  stop("Toutes les images ont déjà été traitées.")
}

# Définir le nombre de cœurs à utiliser
num_cores <- detectCores() - 1  # Garder un cœur libre

# Fonction pour convertir et annoter les images
convert_and_annotate <- function(image_path) {
  # Vérifier si le fichier est au format .heic
  if (tolower(tools::file_ext(image_path)) == "heic") {
    # Lire et convertir l'image HEIC en JPG
    img <- magick::image_read(image_path)
    img <- magick::image_convert(img, format = "jpg")
    # Remplacer l'extension pour le fichier de sortie
    output_path <- sprintf("%s/%s.jpg", annotated_dir, tools::file_path_sans_ext(basename(image_path)))
  } else {
    # Lire l'image si ce n'est pas un fichier HEIC
    img <- magick::image_read(image_path)
    output_path <- sprintf("%s/%s", annotated_dir, basename(image_path))
  }
  
  # Ajouter une annotation à l'image
  site_name <- substr(basename(image_path), 1, nchar(basename(image_path)) - 4)  # Nom du site sans l'extension
  annotated_img <- magick::image_annotate(img, text = site_name, size = 40, color = "red", location = "+30+30")
  
  # Sauvegarder l'image annotée
  magick::image_write(annotated_img, path = output_path)
  return(invisible())
}

# Charger et annoter les images en parallèle
images_annotated <- mclapply(to_process, convert_and_annotate, mc.cores = num_cores)

cat("Traitement terminé pour toutes les images restantes.\n")

#### Install using the devtools package

API_URL <- "https://my-api.plantnet.org/v2/identify"

key <- "2b10fnqc31WV7nagO70MqEfe" # Your API key here
project <- "weurope" # try "weurope" or "canada"

lang <- "fr"
includeRelatedImages <- FALSE # try TRUE

URL <- paste0(API_URL,
              "/", project, "?",
              "lang=", lang,
              "&include-related-images=", includeRelatedImages,
              "&api-key=", key)

image_1 <- "~/Google Drive/Mon Drive/PHOTO/L12C00_CC.jpg"


data <- list(
  "images" = httr::upload_file(image_1),
  "organs" = "flower"
  
)

response <- httr::POST(URL, body=data, encode="multipart")

status <- response$status_code
message(status)

result <- httr::content(response, as = 'parsed')

message(result)
library(httr)
library(jsonlite)
library(dplyr)
results <- result$results

# Créer un dataframe propre à partir des résultats
plant_data <- lapply(results, function(plant) {
  
  # Vérification si les données existent pour chaque plante
  scientific_name <- ifelse(!is.null(plant$species$scientificName), plant$species$scientificName, NA)
  genus <- ifelse(!is.null(plant$species$genus$scientificName), plant$species$genus$scientificName, NA)
  family <- ifelse(!is.null(plant$species$family$scientificName), plant$species$family$scientificName, NA)
  common_names <- ifelse(!is.null(plant$species$commonNames), paste(plant$species$commonNames, collapse = ", "), NA)
  gbif_id <- ifelse(!is.null(plant$gbif$id), plant$gbif$id, NA)
  powo_id <- ifelse(!is.null(plant$powo$id), plant$powo$id, NA)
  iucn_category <- ifelse(!is.null(plant$iucn), plant$iucn$category, NA)
  
  # Créer un dataframe pour chaque plante avec les données vérifiées
  data.frame(
    score = plant$score,
    scientific_name = scientific_name,
    genus = genus,
    family = family,
    common_names = common_names,
    gbif_id = gbif_id,
    powo_id = powo_id,
    iucn_category = iucn_category,
    stringsAsFactors = FALSE
  )
})

# Convertir la liste en un seul dataframe
df <- bind_rows(plant_data)


#BOUCLE COUPER 2####
library(magick)
library(dplyr)

# Dossier contenant les images (y compris les HEIC convertis)
image_folder <- "/Volumes/Untitled1/ARAVOPHTO"  # Remplacez par le chemin vers votre dossier d'images
setwd(image_folder)  # Définir le répertoire de travail vers le dossier des images

# Récupérer la liste des fichiers d'image (y compris jpg, png, et heic)
image_files <- list.files(image_folder, pattern = "\\.jpg$|\\.png$|\\.JPG$|\\.HEIC$", full.names = TRUE)

# Créer un sous-dossier "couper" si il n'existe pas
output_folder <- file.path(image_folder, "couper")
if (!dir.exists(output_folder)) {
  dir.create(output_folder)  # Créer le dossier "couper"
}

# Fonction pour découper chaque image en 4 parties égales et sauvegarder les parties
process_image <- function(image_path, image_index, total_images) {
  # Afficher l'avancement
  message(sprintf("Traitement de l'image %d sur %d : %s", image_index, total_images, basename(image_path)))
  
  # Si l'image est en HEIC, la convertir en JPEG avant de continuer
  if (grepl("\\.heic$", image_path, ignore.case = TRUE)) {
    img <- image_read(image_path)  # Charger directement avec magick si libheif est installé
  } else {
    img <- image_read(image_path)  # Charger les autres formats comme JPG ou PNG
  }
  
  # Obtenir la taille de l'image
  img_width <- image_info(img)$width
  img_height <- image_info(img)$height
  
  # Nom du fichier sans extension
  file_name <- tools::file_path_sans_ext(basename(image_path))
  
  # Découper l'image en 4 parties égales
  
  # Partie haute gauche
  img_top_left <- image_crop(img, geometry = paste0(img_width / 2, "x", img_height / 2, "+0+0"))
  image_write(img_top_left, file.path(output_folder, paste0(file_name, "_top_left.jpg")))
  
  # Partie haute droite
  img_top_right <- image_crop(img, geometry = paste0(img_width / 2, "x", img_height / 2, "+", img_width / 2, "+0"))
  image_write(img_top_right, file.path(output_folder, paste0(file_name, "_top_right.jpg")))
  
  # Partie basse gauche
  img_bottom_left <- image_crop(img, geometry = paste0(img_width / 2, "x", img_height / 2, "+0+", img_height / 2))
  image_write(img_bottom_left, file.path(output_folder, paste0(file_name, "_bottom_left.jpg")))
  
  # Partie basse droite
  img_bottom_right <- image_crop(img, geometry = paste0(img_width / 2, "x", img_height / 2, "+", img_width / 2, "+", img_height / 2))
  image_write(img_bottom_right, file.path(output_folder, paste0(file_name, "_bottom_right.jpg")))
}

# Nombre total d'images à traiter
total_images <- length(image_files)

# Appliquer la fonction à chaque image (y compris les HEIC convertis) avec suivi
lapply(seq_along(image_files), function(i) {
  process_image(image_files[i], i, total_images)
})

lapply(seq(366, total_images), function(i) {
  process_image(image_files[i], i, total_images)
})
message("Toutes les images ont été découpées et sauvegardées dans le dossier 'couper'.")





#BOucle plante V1####
library(httr)
library(jsonlite)
library(dplyr)
library(plantnet)
library(curl)
# URL de l'API
API_URL <- "https://my-api.plantnet.org/v2/identify"

# Clé API
key <- "2b10fnqc31WV7nagO70MqEfe"  # Remplacez par votre clé API

# Projet (choisissez "weurope" ou "canada")
project <- "weurope"

# Langue (français)
lang <- "fr"

# Inclure ou non des images liées
includeRelatedImages <- FALSE

# URL de l'API complète
URL <- paste0(API_URL, "/", project, "?", 
              "lang=", lang,
              "&include-related-images=", includeRelatedImages, 
              "&api-key=", key)

# Dossier contenant les images
image_folder <- "/Volumes/Untitled1/ARAVOPHTO/couper"

# Récupérer la liste des fichiers d'image (par exemple, fichiers jpg et png)
image_files <- list.files(image_folder, pattern = "\\.jpg$|\\.png$", full.names = TRUE)

# Créer un dataframe pour stocker tous les résultats
all_results <- list()

# Parcourir chaque image du dossier
for (image_1 in image_files) {
  # Préparer les données pour la requête
  data <- list(
    "images" = httr::upload_file(image_1),
    "organs" = c("leaf")  # Vous pouvez ajouter plusieurs organes si nécessaire
  )
  Sys.sleep(0.2)
  # Faire la requête POST à l'API
  response <- httr::POST(URL, body = data, encode = "multipart", config = httr::config(http_version = 1))
  
  # Vérifier le code de statut
  status <- response$status_code
  message("Status for image ", image_1, ": ", status)
  
  # Si la requête a réussi (code 200)
  if (status == 200) {
    # Analyser la réponse JSON
    result <- httr::content(response, as = 'parsed')
    
    # Extraire les résultats
    results <- result$results
    
    # Créer un dataframe propre à partir des résultats
    plant_data <- lapply(results, function(plant) {
      scientific_name <- ifelse(!is.null(plant$species$scientificName), plant$species$scientificName, NA)
      genus <- ifelse(!is.null(plant$species$genus$scientificName), plant$species$genus$scientificName, NA)
      family <- ifelse(!is.null(plant$species$family$scientificName), plant$species$family$scientificName, NA)
      common_names <- ifelse(!is.null(plant$species$commonNames), paste(plant$species$commonNames, collapse = ", "), NA)
      gbif_id <- ifelse(!is.null(plant$gbif$id), plant$gbif$id, NA)
      powo_id <- ifelse(!is.null(plant$powo$id), plant$powo$id, NA)
      iucn_category <- ifelse(!is.null(plant$iucn), plant$iucn$category, NA)
      
      # Créer un dataframe pour chaque plante avec les données vérifiées
      data.frame(
        image = image_1,  # Ajouter le nom du fichier image
        score = plant$score,
        scientific_name = scientific_name,
        genus = genus,
        family = family,
        common_names = common_names,
        gbif_id = gbif_id,
        powo_id = powo_id,
        iucn_category = iucn_category,
        stringsAsFactors = FALSE
      )
    })
    
    # Ajouter les résultats de cette image à la liste principale
    all_results <- c(all_results, plant_data)
  } else {
    message("Failed to process image: ", image_1)
  }
}

# Convertir la liste en un seul dataframe
df <- bind_rows(all_results)

# Afficher le dataframe final
print(df)

df$scientific_name


#filter####
library(dplyr)

# Liste des espèces à garder
species_to_keep <- c(
  "Sol nu", "Caillou", "Andro", "Anten", "Alch", "Alch",
  "Alopec", "Carex", 
  "Fest", "Geum ", "Gnap ", "Minuart", "Nardus", "Pilosella",
  "Plan", "Poa", "Potent", "Potent", "Salix ",
  "Senecio", "Sibbal", "Trif alpina","Jacobea"
)

# Filtrer le dataframe en gardant seulement les lignes où scientific_name contient une des espèces de la liste
df_filtered <- df %>%
  filter(sapply(scientific_name, function(x) {
    any(sapply(species_to_keep, function(species) grepl(species, x, ignore.case = TRUE)))
  }))

# Afficher le dataframe filtré
print(df_filtered)


df





#boucle PlanteV2#####
library(httr)
library(jsonlite)
library(dplyr)
library(plantnet)
library(curl)
# URL de l'API
API_URL <- "https://my-api.plantnet.org/v2/identify"

# Clé API
key <- "2b10fnqc31WV7nagO70MqEfe"  # Remplacez par votre clé API

# Projet (choisissez "weurope" ou "canada")
project <- "weurope"

# Langue (français)
lang <- "fr"

# Inclure ou non des images liées
includeRelatedImages <- FALSE

# URL de l'API complète
URL <- paste0(API_URL, "/", project, "?",
              "lang=", lang,
              "&include-related-images=", includeRelatedImages,
              "&api-key=", key)

# Dossier contenant les images
image_folder <- "/Volumes/Untitled1/ARAVOPHTO/couper"

# Récupérer la liste des fichiers d'image (par exemple, fichiers jpg et png)
image_files <- list.files(image_folder, pattern = "\\.jpg$|\\.png$", full.names = TRUE)

# Créer un dataframe pour stocker tous les résultats
all_results <- list()

# Parcourir chaque image du dossier
for (image_1 in image_files) {
  # Préparer les données pour la requête
  data <- list(
    "images" = httr::upload_file(image_1),
    "organs" = c("leaf")  # Vous pouvez ajouter plusieurs organes si nécessaire
  )
  Sys.sleep(0.2)  # Petite pause pour éviter d'envoyer trop de requêtes rapidement
  
  # Faire la requête POST à l'API
  response <- httr::POST(URL, body = data, encode = "multipart", config = httr::config(http_version = 1))
  
  # Vérifier le code de statut
  status <- response$status_code
  message("Status for image ", image_1, ": ", status)
  
  # Si la requête a réussi (code 200)
  if (status == 200) {
    # Analyser la réponse JSON
    result <- httr::content(response, as = 'parsed')
    
    # Extraire les résultats
    results <- result$results
    
    # Créer un dataframe propre à partir des résultats
    plant_data <- lapply(results, function(plant) {
      scientific_name <- ifelse(!is.null(plant$species$scientificName), plant$species$scientificName, NA)
      genus <- ifelse(!is.null(plant$species$genus$scientificName), plant$species$genus$scientificName, NA)
      family <- ifelse(!is.null(plant$species$family$scientificName), plant$species$family$scientificName, NA)
      common_names <- ifelse(!is.null(plant$species$commonNames), paste(plant$species$commonNames, collapse = ", "), NA)
      gbif_id <- ifelse(!is.null(plant$gbif$id), plant$gbif$id, NA)
      powo_id <- ifelse(!is.null(plant$powo$id), plant$powo$id, NA)
      iucn_category <- ifelse(!is.null(plant$iucn), plant$iucn$category, NA)
      
      # Créer un dataframe pour chaque plante avec les données vérifiées
      data.frame(
        image = image_1,  # Ajouter le nom du fichier image
        score = plant$score,
        scientific_name = scientific_name,
        genus = genus,
        family = family,
        common_names = common_names,
        gbif_id = gbif_id,
        powo_id = powo_id,
        iucn_category = iucn_category,
        stringsAsFactors = FALSE
      )
    })
    
    # Ajouter les résultats de cette image à la liste principale
    all_results <- c(all_results, plant_data)
    
  } else {
    # En cas d'erreur (par exemple, HTTP 400), ajouter NA pour cette image
    message("Failed to process image: ", image_1)
    
    # Ajouter une ligne avec NA pour cette image
    all_results <- c(all_results, list(data.frame(
      image = image_1,
      score = NA,
      scientific_name = NA,
      genus = NA,
      family = NA,
      common_names = NA,
      gbif_id = NA,
      powo_id = NA,
      iucn_category = NA,
      stringsAsFactors = FALSE
    )))
  }
}

# Convertir la liste en un seul dataframe
df <- bind_rows(all_results)
setwd("~/Desktop/M2/Lautaret")
write.csv(df, file = "plantnet.csv", row.names = FALSE)
# Afficher le dataframe final
species_to_keep <- c(
  "Sol nu", "Caillou", "Andro", "Anten", "Alch", "Alch",
  "Alopec", "Carex", 
  "Fest", "Geum ", "Gnap ", "Minuart", "Nardus", "Pilosella",
  "Plan", "Poa", "Potent", "Potent", "Salix ",
  "Senecio", "Sibbal", "Trif alpina","Jacobea"
)

# Filtrer le dataframe en gardant seulement les lignes où scientific_name contient une des espèces de la liste
df_filtered <- df %>%
  filter(sapply(scientific_name, function(x) {
    any(sapply(species_to_keep, function(species) grepl(species, x, ignore.case = TRUE)))
  }))

# Afficher le dataframe filtré
print(df_filtered)

#ANNOTER LES PHOTOS COUPER####
library(parallel)
image_dir <- "/Volumes/Untitled1/ARAVOPHTO/couper"
annotated_dir <- sprintf("%s/annotatedcouper", image_dir)

# Créer le répertoire annoté s'il n'existe pas
if (!dir.exists(annotated_dir)) {
  dir.create(annotated_dir, recursive = TRUE)
}

# Lire les chemins des images dans les répertoires source et annoté
image_paths <- list.files(path = image_dir, full.names = TRUE)
annotated_paths <- list.files(path = annotated_dir, full.names = TRUE)

# Exclure les images déjà traitées
annotated_basenames <- basename(annotated_paths)
to_process <- image_paths[!basename(image_paths) %in% annotated_basenames]

# Vérification pour éviter de continuer sans images à traiter
if (length(to_process) == 0) {
  stop("Toutes les images ont déjà été traitées.")
}

# Définir le nombre de cœurs à utiliser
num_cores <- detectCores() - 1  # Garder un cœur libre

# Fonction pour convertir et annoter les images
convert_and_annotate <- function(image_path) {
  # Vérifier si le fichier est au format .heic
  if (tolower(tools::file_ext(image_path)) == "heic") {
    # Lire et convertir l'image HEIC en JPG
    img <- magick::image_read(image_path)
    img <- magick::image_convert(img, format = "jpg")
    # Remplacer l'extension pour le fichier de sortie
    output_path <- sprintf("%s/%s.jpg", annotated_dir, tools::file_path_sans_ext(basename(image_path)))
  } else {
    # Lire l'image si ce n'est pas un fichier HEIC
    img <- magick::image_read(image_path)
    output_path <- sprintf("%s/%s", annotated_dir, basename(image_path))
  }
  
  # Ajouter une annotation à l'image
  site_name <- substr(basename(image_path), 1, nchar(basename(image_path)) - 4)  # Nom du site sans l'extension
  annotated_img <- magick::image_annotate(img, text = site_name, size = 100,boxcolor = "white", color = "red", location = "+30+30")
  
  # Sauvegarder l'image annotée
  magick::image_write(annotated_img, path = output_path)
  return(invisible())
}

# Charger et annoter les images en parallèle
images_annotated <- mclapply(to_process, convert_and_annotate, mc.cores = num_cores)

cat("Traitement terminé pour toutes les images restantes.\n")




#Kriegage GDD ####

library(gstat)
library(sp)
library(sf)
library(ggplot2)
library(viridis)

# Charger les données
data <- read.csv("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/Sol_roche/output data/fichier avec coord/XRF_grille_placette.csv")
data <- data[-43, ]
data <- subset(data, Type.de.sol == "Surface")
head(data$Type.de.sol)

# Nettoyer les colonnes pour supprimer les virgules et convertir en numérique
data$hum.pond <- as.numeric(gsub(",", ".", data$MnO))
data <- data[data$MnO < quantile(data$MnO, 0.95), ]
data$x_l93 <- as.numeric(data$x_l93)
data$y_l93 <- as.numeric(data$y_l93)

# Créer un objet spatial à partir des données
coordinates(data) <- ~x_l93+y_l93
proj4string(data) <- CRS("+proj=utm +zone=31 +datum=WGS84")

# Visualiser les points échantillonnés
plot(data, main = "Sampled Data Points")
# Créer un modèle de variogramme pour l'humidité
variogram_model <- variogram(hum.pond ~ 1, data)
plot(variogram_model, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
fit_variogram <- fit.variogram(variogram_model, model = vgm(psill = 0.0005, model = "Sph", range = 550, nugget = 0.0009))
plot(variogram_model, fit_variogram, main = "Fitted Variogram for Oxyde")
# Définir la zone d'interpolation (grille)
grd <- expand.grid(x = seq(min(data$x_l93), max(data$x_l93), by = 1),
                   y = seq(min(data$y_l93), max(data$y_l93), by = 1))
coordinates(grd) <- ~x+y
gridded(grd) <- TRUE
proj4string(grd) <- proj4string(data)
# Réaliser le krigeage pour interpoler l'humidité du sol
krige_result <- krige(hum.pond ~ 1, data, grd, model = fit_variogram)

# Convertir le résultat en data frame pour ggplot
krige_df <- as.data.frame(krige_result)
colnames(krige_df)[3] <- "prediction"
ggplot(krige_df, aes(x = x, y = y, fill = prediction)) +
  geom_tile() +
  scale_fill_viridis(option = "C") +
  labs(title = "Kriged Interpolation of MnO", fill = "MnO") +
  theme_minimal()

library(raster)

krige_raster <- raster(krige_result)
kmnO<-krige_result
writeRaster(krige_raster, filename = "kmnO.tif", format = "GTiff", overwrite = TRUE)

data <- data[data$MnO < quantile(data$MnO, 0.95), ]

dd<-krige_result@data["var1.pred"]
min(dd)
max(dd)
#transformation data 
library(magick)
library(parallel)

# Charger les bibliothèques nécessaires
library(magick)
library(parallel)

# Chemin vers les images et où sauvegarder les images annotées
image_dir <- "/Volumes/Untitled1/ARAVOPHTO"
annotated_dir <- sprintf("%s/annotated", image_dir)

# Créer le répertoire annoté s'il n'existe pas
if (!dir.exists(annotated_dir)) {
  dir.create(annotated_dir, recursive = TRUE)
}

# Lire les chemins des images
image_paths <- list.files(path = image_dir, full.names = TRUE)

# Vérification pour éviter de continuer sans images
if (length(image_paths) == 0) {
  stop("Aucune image trouvée dans le répertoire spécifié")
}

# Définir le nombre de cœurs à utiliser
num_cores <- detectCores() - 1  # Garder un cœur libre

# Fonction pour convertir et annoter les images
convert_and_annotate <- function(image_path) {
  # Vérifier si le fichier est au format .heic
  if (tolower(tools::file_ext(image_path)) == "heic") {
    # Lire et convertir l'image HEIC en JPG
    img <- magick::image_read(image_path)
    img <- magick::image_convert(img, format = "jpg")
    # Remplacer l'extension pour le fichier de sortie
    output_path <- sprintf("%s/%s.jpg", annotated_dir, tools::file_path_sans_ext(basename(image_path)))
  } else {
    # Lire l'image si ce n'est pas un fichier HEIC
    img <- magick::image_read(image_path)
    output_path <- sprintf("%s/%s", annotated_dir, basename(image_path))
  }
  
  # Ajouter une annotation à l'image
  site_name <- substr(basename(image_path), 1, nchar(basename(image_path)) - 4)  # Nom du site sans l'extension
  annotated_img <- magick::image_annotate(img, text = site_name, size = 100, color = "red", location = "+30+30")
  
  # Sauvegarder l'image annotée
  magick::image_write(annotated_img, path = output_path)
  return(invisible())
}

# Charger et annoter les images en parallèle
images_annotated <- mclapply(image_paths, convert_and_annotate, mc.cores = num_cores)

cat("Traitement terminé pour toutes les images.\n")




image_dir <- "/Volumes/Untitled1/ARAVOPHTO"
annotated_dir <- sprintf("%s/annotated", image_dir)

# Créer le répertoire annoté s'il n'existe pas
if (!dir.exists(annotated_dir)) {
  dir.create(annotated_dir, recursive = TRUE)
}

# Lire les chemins des images dans les répertoires source et annoté
image_paths <- list.files(path = image_dir, full.names = TRUE)
annotated_paths <- list.files(path = annotated_dir, full.names = TRUE)

# Exclure les images déjà traitées
annotated_basenames <- basename(annotated_paths)
to_process <- image_paths[!basename(image_paths) %in% annotated_basenames]

# Vérification pour éviter de continuer sans images à traiter
if (length(to_process) == 0) {
  stop("Toutes les images ont déjà été traitées.")
}

# Définir le nombre de cœurs à utiliser
num_cores <- detectCores() - 1  # Garder un cœur libre

# Fonction pour convertir et annoter les images
convert_and_annotate <- function(image_path) {
  # Vérifier si le fichier est au format .heic
  if (tolower(tools::file_ext(image_path)) == "heic") {
    # Lire et convertir l'image HEIC en JPG
    img <- magick::image_read(image_path)
    img <- magick::image_convert(img, format = "jpg")
    # Remplacer l'extension pour le fichier de sortie
    output_path <- sprintf("%s/%s.jpg", annotated_dir, tools::file_path_sans_ext(basename(image_path)))
  } else {
    # Lire l'image si ce n'est pas un fichier HEIC
    img <- magick::image_read(image_path)
    output_path <- sprintf("%s/%s", annotated_dir, basename(image_path))
  }
  
  # Ajouter une annotation à l'image
  site_name <- substr(basename(image_path), 1, nchar(basename(image_path)) - 4)  # Nom du site sans l'extension
  annotated_img <- magick::image_annotate(img, text = site_name, size = 40, color = "red", location = "+30+30")
  
  # Sauvegarder l'image annotée
  magick::image_write(annotated_img, path = output_path)
  return(invisible())
}

# Charger et annoter les images en parallèle
images_annotated <- mclapply(to_process, convert_and_annotate, mc.cores = num_cores)

cat("Traitement terminé pour toutes les images restantes.\n")

#### Install using the devtools package

API_URL <- "https://my-api.plantnet.org/v2/identify"

key <- "2b10fnqc31WV7nagO70MqEfe" # Your API key here
project <- "weurope" # try "weurope" or "canada"

lang <- "fr"
includeRelatedImages <- FALSE # try TRUE

URL <- paste0(API_URL,
              "/", project, "?",
              "lang=", lang,
              "&include-related-images=", includeRelatedImages,
              "&api-key=", key)

image_1 <- "~/Google Drive/Mon Drive/PHOTO/L12C00_CC.jpg"


data <- list(
  "images" = httr::upload_file(image_1),
  "organs" = "flower"
  
)

response <- httr::POST(URL, body=data, encode="multipart")

status <- response$status_code
message(status)

result <- httr::content(response, as = 'parsed')

message(result)
library(httr)
library(jsonlite)
library(dplyr)
results <- result$results

# Créer un dataframe propre à partir des résultats
plant_data <- lapply(results, function(plant) {
  
  # Vérification si les données existent pour chaque plante
  scientific_name <- ifelse(!is.null(plant$species$scientificName), plant$species$scientificName, NA)
  genus <- ifelse(!is.null(plant$species$genus$scientificName), plant$species$genus$scientificName, NA)
  family <- ifelse(!is.null(plant$species$family$scientificName), plant$species$family$scientificName, NA)
  common_names <- ifelse(!is.null(plant$species$commonNames), paste(plant$species$commonNames, collapse = ", "), NA)
  gbif_id <- ifelse(!is.null(plant$gbif$id), plant$gbif$id, NA)
  powo_id <- ifelse(!is.null(plant$powo$id), plant$powo$id, NA)
  iucn_category <- ifelse(!is.null(plant$iucn), plant$iucn$category, NA)
  
  # Créer un dataframe pour chaque plante avec les données vérifiées
  data.frame(
    score = plant$score,
    scientific_name = scientific_name,
    genus = genus,
    family = family,
    common_names = common_names,
    gbif_id = gbif_id,
    powo_id = powo_id,
    iucn_category = iucn_category,
    stringsAsFactors = FALSE
  )
})

# Convertir la liste en un seul dataframe
df <- bind_rows(plant_data)


