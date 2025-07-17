#DATA####

library(readxl)
library(readr)
library(ade4)
library(FactoMineR)
library(factoextra)
library(vegan)
library(dplyr)
a2001 <- read_excel("~/Downloads/ARAVO_RESURVEYxlsx.xlsx", 
                                 sheet = "2001_GRID.SPEC", col_types = c("text", 
                                                                         "text", "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "numeric", "numeric", "numeric", 
                                                                         "skip"))

##Présence-Absence####
a2022 <- read_csv("~/Downloads/ARAVO_RESURVEYxlsx.xlsx - 2022_GRID.SPEC.csv")
##traits#####
data("aravo")
traits<-aravo$traits

#Etre sur que ca match bien####
# Noms des espèces dans la matrice `traits`
names_traits <- rownames(traits)

# Noms des espèces dans `a2022`
names_a2022 <- colnames(a2022)[3:ncol(a2022)]  # Supposant que les espèces commencent après les 2 premières colonnes
not_in_a2022 <- setdiff(names_traits, names_a2022)

# Espèces dans `a2022` mais pas dans `traits`
not_in_traits <- setdiff(names_a2022, names_traits)

# Afficher les résultats
cat("Espèces dans `traits` mais pas dans `a2022` :\n", not_in_a2022, "\n")
cat("Espèces dans `a2022` mais pas dans `traits` :\n", not_in_traits, "\n")
# Renommer les colonnes ou lignes si nécessaire
colnames(a2022) <- gsub(" ", ".", colnames(a2022))  # Remplace les espaces par des points
rownames(traits) <- gsub(" ", ".", rownames(traits))  # Fait de même dans `traits`
# Garder uniquement les espèces communes
# Corriger les noms dans `a2022` (remplacez "ancien_nom" par le bon nom correspondant)
colnames(a2022) <- gsub("Carex.curvula.subsp..rosea", "Care.rosa", colnames(a2022))
colnames(a2022) <- gsub("Potent.aurea", "Pote.aure", colnames(a2022))
colnames(a2022) <- gsub("Plan.alpi", "Plan.alpi", colnames(a2022))
colnames(a2022) <- gsub("Trif.alpina", "Trif.alpi", colnames(a2022))
colnames(a2022) <- gsub("Alopec.alpi", "Alop.alpi", colnames(a2022))
colnames(a2022) <- gsub("Care.myos", "Kobr.myos", colnames(a2022))
colnames(a2022) <- gsub("Andro.vital", "Andr.vita", colnames(a2022))
colnames(a2022) <- gsub("Senecio.incanus/Jacobea.incana", "Sene.inca", colnames(a2022))
colnames(a2022) <- gsub("Sali.serpy", "Sali.serp", colnames(a2022))
colnames(a2022) <- gsub("Sibbal.proc", "Sibb.proc", colnames(a2022))
colnames(a2022) <- gsub("Potent.crantzii", "Pote.cran", colnames(a2022))
colnames(a2022) <- gsub("Poa.alpina", "Poa.alpi", colnames(a2022))
colnames(a2022) <- gsub("Pilosella.officinarum", "Hier.pili", colnames(a2022))
colnames(a2022) <- gsub("Minuart.sed", "Minu.sedo", colnames(a2022))
colnames(a2022) <- gsub("Gnap.supi", "Omal.supi", colnames(a2022))

# Répétez pour tous les noms divergents identifiés
names_a2022 <- colnames(a2022)[3:ncol(a2022)]
common_species <- intersect(names_traits, names_a2022)

# Filtrer les données pour ne conserver que les espèces communes
traits_aligned <- traits[common_species, ]
a2022_aligned <- a2022[, c(1:2, which(colnames(a2022) %in% common_species))]
acp_result <- PCA(traits_aligned, scale.unit = TRUE, graph = FALSE)
# Graphique des individus
# Biplot combinant individus et variables
fviz_pca_biplot(acp_result, repel = TRUE) +
  ggtitle("ACP")
#MERGING GRID AND OBS####
ARAVOGRID2022 <- read_csv("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/spatialisation/data loc/grille/ARAVOGRID2022.csv")
a2022
ARAVOGRID2022$lig<-ARAVOGRID2022$ROW
ARAVOGRID2022$col<-ARAVOGRID2022$COLUMN


a2022_aligned$lig <- as.numeric(a2022_aligned$lig)
a2022_aligned$col <- as.numeric(a2022_aligned$col)
ARAVOGRID2022$lig <- as.numeric(ARAVOGRID2022$lig)
ARAVOGRID2022$col <- as.numeric(ARAVOGRID2022$col)

# Fusionner les deux datasets sur les colonnes lig et col
a2022_enriched <- merge(a2022_aligned, ARAVOGRID2022, by = c("lig", "col"), all.x = TRUE)

# Afficher un aperçu des données enrichies
print(head(a2022_enriched))
# Supprimer les lignes contenant un 0 dans toutes les colonnes à partir de la 4ème colonne



a2022_enriched$X[a2022_enriched$X == 0] <- NA
a2022_filtered <- a2022_enriched[!is.na(a2022_enriched$X), ]
prepRLQ <- a2022_enriched[!is.na(a2022_enriched$X), ]
a2022_filtered[is.na(a2022_filtered)] <- 0
# Vérifier que les noms des espèces dans traits_aligned correspondent aux colonnes de a2022_filtered
traits_aligned$Species <- rownames(traits_aligned)  # Ajouter une colonne "Species" dans traits_aligned
species_columns <- intersect(colnames(a2022_filtered)[3:ncol(a2022_filtered)], traits_aligned$Species)

# Ajouter les traits fonctionnels pondérés pour chaque ligne
for (trait in colnames(traits_aligned)[1:(ncol(traits_aligned) - 1)]) {  # Boucle sur chaque trait
  # Calculer la moyenne pondérée des traits pour chaque ligne
  a2022_filtered[[trait]] <- apply(a2022_filtered[, species_columns], 1, function(row) {
    species_present <- names(row)[row > 0]  # Identifier les espèces présentes (valeur > 0)
    if (length(species_present) == 1) {
      # Si une seule espèce est présente, retourner directement la valeur du trait
      traits_aligned[traits_aligned$Species == species_present, trait]
    } else if (length(species_present) > 1) {
      # Moyenne des traits pour les espèces présentes
      mean(traits_aligned[traits_aligned$Species %in% species_present, trait], na.rm = TRUE)
    } else {
      NA  # Si aucune espèce n'est présente, retourner NA
    }
  })
}

# Vérifier un aperçu des données finales
print(head(a2022_filtered))
a2022_filtered <- na.omit(a2022_filtered)

prepRLQ<-a2022_filtered
prepRLQ<-as.data.frame(prepRLQ)



#RLQ####
##MERG####
###HUMI####

str(khumidit)
prepRLQ
krige_raster <- raster(khumidit, layer = "var1.pred") # Utiliser la couche de prédiction
# Créer un objet SpatialPoints avec les coordonnées de prepRLQ
coords <- prepRLQ[, c("X", "Y")]

# Extraire les valeurs du raster aux points donnés par prepRLQ
prepRLQ$humi <- extract(krige_raster, coords)

###DRY####
str(kdry)
prepRLQ
krige_raster <- raster(kdry, layer = "var1.pred") # Utiliser la couche de prédiction
# Créer un objet SpatialPoints avec les coordonnées de prepRLQ
coords <- prepRLQ[, c("X", "Y")]

# Extraire les valeurs du raster aux points donnés par prepRLQ
prepRLQ$kdry <- extract(krige_raster, coords)

###MNO####
str(kmnO)
prepRLQ
krige_raster <- raster(kmnO, layer = "var1.pred") # Utiliser la couche de prédiction
# Créer un objet SpatialPoints avec les coordonnées de prepRLQ
coords <- prepRLQ[, c("X", "Y")]

# Extraire les valeurs du raster aux points donnés par prepRLQ
prepRLQ$kmnO <- extract(krige_raster, coords)

###MO####
str(kmo)
prepRLQ
krige_raster <- raster(kmo, layer = "var1.pred") # Utiliser la couche de prédiction
# Créer un objet SpatialPoints avec les coordonnées de prepRLQ
coords <- prepRLQ[, c("X", "Y")]

# Extraire les valeurs du raster aux points donnés par prepRLQ
prepRLQ$mo <- extract(krige_raster, coords)

###ph####
kph<-raster("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/spatialisation/krigeage/krig_pH.tif")
str(kph)
 # Utiliser la couche de prédiction
# Créer un objet SpatialPoints avec les coordonnées de prepRLQ
coords <- prepRLQ[, c("X", "Y")]

# Extraire les valeurs du raster aux points donnés par prepRLQ
prepRLQ$ph <- extract(kph, coords)

###SWI####
kswi<-raster("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/spatialisation/krigeage/script phillipe/SWI.tif")

prepRLQ$swi <- extract(kswi, coords)

###resit####
kresist<-raster("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/spatialisation/krigeage/script phillipe/resistivity_cmd_aravo.tif")

prepRLQ$resist <- extract(kresist, coords)

###convex####
kconvex<-raster("~/Google Drive/Mon Drive/aravo - stage clementine Bellet/spatialisation/krigeage/script phillipe/ARA_CONVEX_1M.tif")
prepRLQ$convex <- extract(kconvex, coords)

###GDD####
kgdd<-raster("~/Downloads/krige_GDD.tif")
prepRLQ$gdd <- extract(kgdd, coords)


#Doing####

str(prepRLQ)
prepRLQ$Kobr.myos <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Kobr.myos))
prepRLQ <- prepRLQ[, !names(prepRLQ) %in% "kmnO"]
prepRLQ<-na.omit(prepRLQ)


# Identifier les colonnes correspondant à chaque catégorie
# Espèces : Par exemple, colonnes 3 à 23
str(prepRLQ)
prepRLQ$Kobr.myos <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Kobr.myos))
prepRLQ$Care.semp <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Care.semp))
prepRLQ$Fest.viol <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Fest.viol))
prepRLQ$Hier.pili <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Hier.pili))
prepRLQ$Pote.aure <- as.numeric(gsub("[^0-9.-]", "", prepRLQ$Pote.aure))

prepRLQ <- prepRLQ[, !names(prepRLQ) %in% "kmnO"]
prepRLQ<-na.omit(prepRLQ)

species_columns <- prepRLQ[, 3:23]

# Traits fonctionnels : Colonnes des traits (par exemple, 34 à 42)
trait_columns <- prepRLQ[, c("Height", "Spread", "Angle", "Area", "Thick", "SLA", "N_mass", "Seed")]

# Variables environnementales : Les colonnes restantes correspondant à l'environnement
env_columns <- prepRLQ[, c("humi", "kdry", "mo","ph","swi","resist","convex","gdd")]
colnames(env_columns) <- c("Humidité", "Dry", "Matière organique", "pH", "SWI", "Resistivité", "Convexité", "Growing_Degree_Days")
# Mettre chaque groupe dans une liste
lists_prepRLQ <- list(
  Environmental = env_columns,
  Species = species_columns,
  Traits = traits_aligned[, -ncol(traits_aligned)]
)

coa.abu <- dudi.coa(lists_prepRLQ$Species, scannf = FALSE, nf=2)
pca.trait <- dudi.pca(lists_prepRLQ$Traits, scannf = FALSE, 
                      row.w = coa.abu$cw)
pca.env <- dudi.pca(lists_prepRLQ$Environmental, scannf = FALSE, 
                    row.w = coa.abu$lw)
rlqF <- rlq(pca.env, coa.abu, pca.trait, 
            scannf = FALSE)
str(rlqF)
plot(rlqF)
summary(rlqF)
str(rlqF)

nrepet <- 49999# the higher, the longer the run lasts
four.comb.aravo <- fourthcorner(lists_prepRLQ$Environmental, lists_prepRLQ$Species,
                                lists_prepRLQ$Traits, modeltype = 5, p.adjust.method.G = "none",
                                p.adjust.method.D = "none", nrepet = nrepet)

summary(four.comb.aravo)
par(mfrow=c(1,1))
fr<-summary(four.comb.aravo)
str(fr)
colnames(fr)
fr_significatif <- subset(fr, Pvalue <= 0.05)
r_significatif <- fr[as.numeric(fr[, "Pvalue"]) <= 0.05, ]
plot(four.comb.aravo, alpha = 0.05, stat = "D2")
testrlq.aravo <- randtest(rlqF, modeltype = 5, nrepet = nrepet)
plot(testrlq.aravo)

plot(four.comb.aravo, x.rlq = rlqF, alpha = 0.05,
     stat = "D2", type = "biplot")

##PLOTING####
par(mfrow = c(1, 3)) 
d <- data.frame(
  RS1 = c(0.4332638, -0.1778495, 0.1405079, -0.4434502, 0.2525230, 0.3459701, -0.4891563, -0.3761972),
  RS2 = c(0.5814719, -0.3491318, 0.1226737, 0.2747810, -0.1733238, -0.3103077, 0.4165839, -0.3867426),
  row.names = c("Humidité", "Dry", "Matière organique", "pH", "SWI", "Resistivité", "Convexité", "GDD")
)
s.arrow(d) 
s.arrow(rlqF$c1) 
s.label(rlqF$lQ, boxes = FALSE)
###test####
# Création des données

# Vérification des données
print(data)


#BUILDING .TIF FOR FONCTIONNAL TRAITS IN ARAVO ####
library(gstat)
library(sp)
##SLA####
prepRLQ
prepRLQ<-as.data.frame(prepRLQ)
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$SLA <- as.numeric(prepRLQ$SLA)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(SLA ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(SLA ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(SLA ~ 1, prepRLQ, grid, model = model_variogram)
Rsla <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "SLA") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de SLA",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")


plot(prepRLQ, main = "Sampled Data Points")

##





##Height####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Height")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Height <- as.numeric(prepRLQ$Height)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Height ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Height ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Height ~ 1, prepRLQ, grid, model = model_variogram)
Rheight<- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Height") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Height",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##








##Seed####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Seed")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Seed <- as.numeric(prepRLQ$Seed)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Seed ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph", 10, 10))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Seed ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Seed ~ 1, prepRLQ, grid, model = model_variogram)
Rseed <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Seed") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Seed",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##








##Spread####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Spread")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Spread <- as.numeric(prepRLQ$Spread)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Spread ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Gau", 10, 100))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Spread ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Spread ~ 1, prepRLQ, grid, model = model_variogram)
Rspread <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Spread") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Spread",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##








##N_mass####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "N_mass")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$N_mass <- as.numeric(prepRLQ$N_mass)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(N_mass ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(N_mass ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(N_mass ~ 1, prepRLQ, grid, model = model_variogram)

RN_mass <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "N_mass") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de N_mass",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##









##Thick####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Thick")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Thick <- as.numeric(prepRLQ$Thick)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Thick ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Thick ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Thick ~ 1, prepRLQ, grid, model = model_variogram)
Rthick <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Thick") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Thick",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##











##Area####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Area")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Area <- as.numeric(prepRLQ$Area)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Area ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Area ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)

# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Area ~ 1, prepRLQ, grid, model = model_variogram)
Rarea <- raster(kriging_result, layer = "var1.pred") 
# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Area") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Area",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")




##









##Angle####
str(prepRLQ)
prepRLQ<-as.data.frame(prepRLQ)
print(head(prepRLQ[, c("X", "Y", "Angle")]))  # Vérifiez X, Y, SLA
prepRLQ$X <- as.numeric(prepRLQ$X)
prepRLQ$Y <- as.numeric(prepRLQ$Y)
prepRLQ$Angle <- as.numeric(prepRLQ$Angle)

# Créer un objet spatial à partir de `a2022_filtered`
coordinates(prepRLQ) <- ~X + Y

# Étape 2 : Calcul du variogramme
variogram_SLA <- variogram(Angle ~ 1, data = prepRLQ)
print(variogram_SLA)

# Ajuster un modèle de variogramme
model_variogram <- fit.variogram(variogram_SLA, model = vgm("Sph"))
print(model_variogram)
# Créer un modèle de variogramme pour l'humidité
variogram_SLA <- variogram(Angle ~ 1, prepRLQ)
plot(variogram_SLA, main = "Empirical Variogram")

# Ajuster un modèle théorique au variogramme empirique (modèle sphérique)
plot(variogram_SLA, model_variogram, main = "Fitted Variogram for Oxyde")
# Étape 3 : Définir une grille pour l'interpolation
x_range <- seq(min(prepRLQ$X), max(prepRLQ$X), length.out = 1000)
y_range <- seq(min(prepRLQ$Y), max(prepRLQ$Y), length.out = 1000)
grid <- expand.grid(X = x_range, Y = y_range)
coordinates(grid) <- ~X + Y
gridded(grid) <- TRUE

# Étape 4 : Effectuer le krigeage
kriging_result <- krige(Angle ~ 1, prepRLQ, grid, model = model_variogram)

# Étape 5 : Visualiser les résultats
# Convertir en data frame pour visualisation
kriging_result_df <- as.data.frame(kriging_result)

# Visualiser avec ggplot2
library(ggplot2)
ggplot(data = kriging_result_df, aes(x = X, y = Y, fill = var1.pred)) +
  geom_tile() +
  scale_fill_gradient(low = "blue", high = "red", name = "Angle") +
  theme_minimal() +
  labs(title = "Kriging des valeurs de Angle",
       x = "X (Coordonnées)",
       y = "Y (Coordonnées)")
raster_result <- raster(kriging_result, layer = "var1.pred") 
Rangle <- raster(kriging_result, layer = "var1.pred") 
# 
writeRaster(raster_result, filename = "angle.tif", format = "GTiff", overwrite = TRUE)

##


#test####
# Convertir les rasters en une data.frame (exclut les valeurs NA automatiquement)
df <- as.data.frame(c(Rangle, RN_mass,Rseed,Rsla,Rheight,Rthick,Rarea,Rspread), na.rm = TRUE)
raster_stack <- stack(Rangle, RN_mass,Rseed,Rsla,Rheight,Rthick,Rarea,Rspread)
# Renommer les colonnes pour clarté
# Convertir les rasters en data.frame
df <- as.data.frame(raster_stack, na.rm = TRUE)

# Renommer les colonnes pour plus de clarté
colnames(df) <- c("Rangle", "RN_mass","Rseed","Rsla","Rheight","Rthick","Rarea","Rspread")

# Visualiser un échantillon des données
head(df)

# Visualiser un échantillon des données
head(df)
summary(lm(Rsla~RN_mass,data=df))
plot(Rsla~RN_mass,data=df)
summary(lm(Rsla~RN_mass,data=df))
# Charger ggplot2
library(ggplot2)

# Créer le graphique avec ggplot2
ggplot(df, aes(x = RN_mass, y = Rangle)) +
  geom_point(color = "blue", alpha = 0.6) +  # Points pour les données
  geom_smooth(method = "lm", formula = y ~ x, color = "red", se = TRUE) +  # Ligne de régression
  labs(title = "Relation linéaire entre Rangle et RN_mass",
       x = "RN_mass",
       y = "Rangle") +
  theme_minimal()
#fin test####
 
#VEGEDIST SP####
prepRLQ<-as.data.frame(prepRLQ)

poisson
species_data <- prepRLQ[, 3:23]

distance_matrix <- vegdist(species_data, method = "jaccard")
cluster <- hclust(distance_matrix, method = "ward.D2") #
  fviz_nbclust(species_data, kmeans, method = "wss") + 
  labs(title = "Méthode du coude pour déterminer le nombre de clusters") +
  theme_minimal() 
#K=8

plot(cluster, main = "Dendrogramme des communautés", xlab = "Sites")

rect.hclust(cluster, k = 8, border = "red") 
  nmds <- metaMDS(species_data, distance = "jaccard", k = 8, trymax = 100)
plot(nmds, type = "t", main = "NMDS des communautés")
groups <- cutree(cluster, k = 8)
ordiplot(nmds, type = "n", main = "NMDS des communautés par groupes")
points(nmds, display = "sites", col = groups, pch = 19)  
legend("topright", legend = unique(groups), col = unique(groups), pch = 15, title = "Groupes")

groups <- cutree(cluster, k = 4)

# Ajouter les groupes en tant que nouvelle colonne au dataframe prepRLQ
prepRLQ$community <- as.factor(groups)
summary(prepRLQ$community)
#Oui ont distingue 5 communauté piscicole pendant le suivit 

community_matrix <- matrix %>%
  pivot_wider(names_from = Espece, values_from = Abondance, values_fill = 0)
species_matrix <- community_matrix[, -c(1:6)]  

distance_matrix <- vegdist(species_matrix, method = "bray")
pcoa_res <- pcoa(distance_matrix)

pcoa_df <- as.data.frame(pcoa_res$vectors)
pcoa_df$M <- community_matrix$M 
ggplot(pcoa_df, aes(x = Axis.1, y = Axis.2, color = M)) +
  geom_point(size = 4) +
  labs(title = "Analyse PCoA des communautés piscicoles",
       x = "Axe 1", y = "Axe 2") +
  theme_minimal()

# PERMANOVA pour tester les différences entre groupes
adonis_res <- adonis2(distance_matrix ~ M, data = community_matrix)
print(adonis_res)

# et les communauté piscicole semble varier en fonction des mois (p-value<0.05 , PERMANOVA)

aravo$traits

#VEGEDIST Trait####
prepRLQ

poisson
species_data <- prepRLQ[, 32:39]

distance_matrix <- vegdist(species_data, method = "bray")
cluster <- hclust(distance_matrix, method = "ward.D2") #
fviz_nbclust(species_data, kmeans, method = "wss") + 
  labs(title = "Méthode du coude pour déterminer le nombre de clusters") +
  theme_minimal() 
#K=8

plot(cluster, main = "Dendrogramme des communautés", xlab = "Sites")

rect.hclust(cluster, k = 4, border = "red") 
nmds <- metaMDS(species_data, distance = "bray", k = 4, trymax = 100)
plot(nmds, type = "t", main = "NMDS des communautés")


nmds_sites <- as.data.frame(scores(nmds, display = "sites")) # Coordonnées des sites
nmds_sites$site <- rownames(nmds_sites) # Ajouter une colonne pour les noms des sites

nmds_species <- as.data.frame(scores(nmds, display = "species")) # Coordonnées des espèces
nmds_species$species <- rownames(nmds_species) # Ajouter une colonne pour les noms des espèces
# Graphique NMDS
ggplot() +
  # Points pour les sites
  geom_point(data = nmds_sites, aes(x = NMDS1, y = NMDS2, color = "black"), size = 3) +
  # Labels des sites
  geom_text(data = nmds_sites, aes(x = NMDS1, y = NMDS2, label = site), vjust = -1, size = 3) +
  # Vecteurs pour les espèces
  geom_segment(data = nmds_species, aes(x = 0, y = 0, xend = NMDS1, yend = NMDS2),
               arrow = arrow(length = unit(0.2, "cm")), color = "red") +
  # Labels des espèces
  geom_text(data = nmds_species, aes(x = NMDS1, y = NMDS2, label = species), color = "red", size = 3) +
  # Esthétique
  theme_minimal() +
  labs(title = "NMDS Ordination", x = "NMDS1", y = "NMDS2") +
  theme(plot.title = element_text(hjust = 0.5))

groups <- cutree(cluster, k = 4)
ordiplot(nmds, type = "n", main = "NMDS des communautés par groupes")
points(nmds, display = "sites", col = groups, pch = 19)  
legend("topright", legend = unique(groups), col = unique(groups), pch = 15, title = "Groupes")

groups <- cutree(cluster, k = 4)



# NMDS avec les communautés par groupes et les espèces
ordiplot(nmds, type = "n", main = "NMDS des communautés par groupes")

# Ajouter les points pour les sites (communautés)
points(nmds, display = "sites", col = groups, pch = 19)

# Ajouter les flèches ou noms pour les espèces
text(nmds, display = "species", col = "blue", cex = 0.8) # Ajustez `col` et `cex` si nécessaire

# Ajouter la légende
legend("topright", legend = unique(groups), col = unique(groups), pch = 15, title = "Groupes")

# Ajouter les groupes en tant que nouvelle colonne au dataframe prepRLQ
prepRLQ$community <- as.factor(groups)
summary(prepRLQ$community)
##plot####


# Liste des variables à analyser
variables <- c("ph", "gdd", "convex", "swi", "resist", "mo", "kdry", "humi")

# Boucle pour réaliser l'ANOVA et les visualisations
results <- list()

for (var in variables) {
  # Créer une formule dynamique
  formula <- as.formula(paste(var, "~ community"))
  
  # Ajuster le modèle linéaire et l'ANOVA
  lm_model <- lm(formula, data = prepRLQ)
  anova_model <- aov(lm_model)
  
  # Stocker les résultats dans une liste
  results[[var]] <- list(
    summary = summary(anova_model),
    tukey = TukeyHSD(anova_model)
  )
  
  # Afficher les résultats ANOVA et test de Tukey
  cat("\n### Résultats pour :", var, "###\n")
  print(summary(lm_model))
  print(summary(anova_model))
  print(TukeyHSD(anova_model))
  
  # Visualisation des boxplots
  # Imprime correctement les graphiques dans une boucle
}




# Charger les librairies nécessaires
library(ggplot2)
library(patchwork)
# Liste des noms lisibles
variable_labels <- c(
  ph = "pH",
  gdd = "Growing Degree Days",
  convex = "Convexité",
  swi = "SWI",
  resist = "Resistivité",
  mo = "Matière organique",
  kdry = "Dry",
  humi = "Humidité"
)
# Liste des variables et couleurs sobres
variables <- c("ph", "gdd", "convex", "swi", "resist", "mo", "kdry", "humi")
colors <- c("#4E79A7", "#F28E2B", "#E15759", "#76B7B2") # Palette sobre

# Fonction pour créer les graphiques avec ANOVA
create_plot <- function(variable, data, colors, labels) {
  # Formule et ANOVA
  formula <- as.formula(paste(variable, "~ community"))
  anova_model <- aov(formula, data = data)
  pval <- summary(anova_model)[[1]]["Pr(>F)"][1] # Extraire la p-valeur
  
  # Récupérer le label lisible pour la variable
  variable_label <- labels[[variable]]
  
  # Créer le graphique ggplot
  plot <- ggplot(data, aes(x = community, y = !!sym(variable), fill = community)) +
    geom_boxplot(outlier.shape = NA) + # Ajout de points
    scale_fill_manual(values = colors, name = "Communautés") + # Palette de couleurs sobres
    theme_minimal() +
    labs(
      title = variable_label, # Utiliser le label lisible
      subtitle = paste0("ANOVA p-value < ", format(round(pval[1,], 4), scientific = TRUE)),
      x = "Communautés",
      y = variable_label
    ) +
    theme(
      plot.title = element_text(hjust = 0.5, face = "bold"),
      plot.subtitle = element_text(hjust = 0.5, size = 10, face = "italic")
    )
  
  return(plot)
}

# Créer les graphiques pour toutes les variables avec leurs noms lisibles
plots <- lapply(variables, create_plot, data = prepRLQ, colors = colors, labels = variable_labels)

# Combiner les graphiques dans une seule figure
combined_plot <- wrap_plots(plots, ncol = 4)

# Afficher le graphique combiné
print(combined_plot)




#ACP TRAIT VS ENV####

prepRLQ


traits <- prepRLQ[, c("Height", "Spread", "Angle", "Area", "Thick", "SLA", "N_mass", "Seed")]
env_variables <- prepRLQ[, c("humi", "kdry", "mo", "ph", "swi", "resist", "convex", "gdd")]

# Vérifier les données
summary(traits)
summary(env_variables)

# Standardiser les données pour éviter que les échelles affectent les résultats
traits_std <- scale(traits)
env_variables_std <- scale(env_variables)

# Réalisation de l'ACP
res_pca <- PCA(traits_std, scale.unit = TRUE, graph = FALSE)  # Avec FactoMineR
res_rlq <- dudi.pca(traits_std, scannf = FALSE, nf = 2)  # Avec ade4 pour PCA basique

# Visualisation des résultats
# Graphique des individus
fviz_pca_ind(res_pca, col.ind = "cos2", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Graphique des variables (traits fonctionnels)
fviz_pca_var(res_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Graphique biplot
fviz_pca_biplot(res_pca, repel = TRUE, col.var = "#2E9FDF", col.ind = "#696969")

# Graphique pour variables environnementales (optionnel, analyse séparée)
res_env_pca <- PCA(env_variables_std, scale.unit = TRUE, graph = FALSE)
fviz_pca_var(res_env_pca, col.var = "contrib", gradient.cols = c("#00AFBB", "#E7B800", "#FC4E07"),
             repel = TRUE)

# Résultats ACP résumé
print(res_pca)

str(species_data)

#LM####
prepRLQ
##TRAIT####
summary(lm(SLA~N_mass,data=prepRLQ))
summary(lm(Seed~Height,data=prepRLQ))
summary(lm(Area~Spread,data=prepRLQ))
summary(lm(Angle~Thick,data=prepRLQ))


##ENV####
summary(lm(humi~gdd,data=prepRLQ))
summary(lm(convex~swi,data=prepRLQ))
summary(lm(convex~ph,data=prepRLQ))
summary(lm(swi~ph,data=prepRLQ))
##2####
summary(lm(Angle~ph*convex,data=prepRLQ))

