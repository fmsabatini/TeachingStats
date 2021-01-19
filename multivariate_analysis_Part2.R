## ---- eval=T-------------------------------------------------------------------------------------------
# download packages, if necessary
install.packages(c("vegan",  "FD"))


## ---- warning=F, message=F-----------------------------------------------------------------------------
#load packages
library(vegan)
library(FD)


## ---- eval=T-------------------------------------------------------------------------------------------
# clean workspace
rm(list=ls())
#reimport data
data(dune)
data(dune.env)
data(varespec)
data(varechem)
data(BCI)
data(BCI.env)


## ---- eval=T-------------------------------------------------------------------------------------------
#hints
?kmeans
?tapply
?hclust
?rect.hclust
?cophenetic
?cutree
?cluster::agnes
?FD::gowdis


## ---- eval=T, echo=F-----------------------------------------------------------------------------------
# Non-hierarchical classification ####
dune.dist <- vegdist(dune, method="bray") # the distance matrix
dune.kmeans <- kmeans(dune, centers=5, iter.max=10000) # k-means method

## mean A1 per group
a1.mean <- tapply(dune.env$A1, dune.kmeans$cluster, "mean")
boxplot(dune.env$A1 ~ dune.kmeans$cluster, col=1:5, ylab="A1 thickness",
        xlab="Group (k-means)")

# Hierarchical classification of species data
dune.single <- hclust(dune.dist, method="single")
dune.complete <- hclust(dune.dist, method="complete")
dune.average <- hclust(dune.dist, method="average")

# Compare the classifications graphically
par(mfrow=c(3,1)) ## arrange three plots into one graph
plot(dune.single)
plot(dune.complete)
plot(dune.average)
dev.off() ## deactivates multi-plotting
# compute means of A1 when cutting the dendrogram to 5 groups
plot(dune.average)
rect.hclust(tree = dune.average, 5)
groups <- cutree(dune.average, 5)
tapply(dune.env$A1, groups, "mean")
boxplot(dune.env$A1 ~ groups, col=1:5, xlab="groups (UPGMA)", ylab="A1 thickness")

## evaluate classifications
coph.single <- cophenetic(dune.single)
coph.complete <- cophenetic(dune.complete)
coph.average <- cophenetic(dune.average)
plot(dune.dist, coph.single) 
# correlates original dissimilarities with the dissimilarity of the tree
cor(dune.dist,coph.single)
cor(dune.dist,coph.complete)
cor(dune.dist,coph.average)

## dissimilarity with mixed data - use Gower dissimilarity
dune.env.gowdis <- gowdis(dune.env)
dune.env.average <- hclust(dune.env.gowdis, method="average")
par(mfrow=c(1,2)) # compare clustering based on sp vs. env. data
plot(dune.average)
plot(dune.env.average)
dev.off()

coph.env.average <- cophenetic(dune.env.average)
cor(coph.env.average, coph.average)



## ---- eval=T-------------------------------------------------------------------------------------------
#Hints
?rda
?screeplot
str(output.rda) #replace 'output.rda' with the name of your rda object 
#You can extract the eigenvalues as
output.rda$CA$eig
?cumsum
?plot.cca
?biplot
?metaMDS
?envfit


## ---- eval=T, echo=F-----------------------------------------------------------------------------------
# PCA - Principal Components Analysis ####
## 1) The basic method (no scaling !!! Why is this wrong?)
vare.pca <- rda(varechem) # Note the function is called "rda" in vegan
vare.pca
# Total intertia is the variance, see decreasing eigenvalues of axes
# Question: how much variance does the first axis explain?
# Answer: 64057/85655 = 74.8 %
screeplot(vare.pca)
cumsum(vare.pca$CA$eig) / sum(vare.pca$CA$eig)  *100
# The first 4 axes explain already 90% of total variation
plot(vare.pca, scaling = "species")
#the term 'species' is a confusing legacy. It should be intended as 'variables'
plot(vare.pca, scaling = "sites") 

## 2) PCA with scaling
vare.pca.st <- rda(varechem, scale = TRUE) # Now the variance of species is 
#                     standarized and inertia is correlation instead of variance
vare.pca.st    #note that numbers differ, and the inertia 
               #is maximized to the numnber of species
# "3" scales both species and sites; "2" by sites, "1" by species(=variables)
plot(vare.pca.st, scaling = 2)
biplot (vare.pca.st, display = 'species')
plot(vare.pca.st, scaling = 1)
biplot (vare.pca.st, display = 'sites')

plot(vare.pca.st, choices=c(1,3), scaling = 3) 
biplot(vare.pca.st, choices=c(1,3), scaling = 3, display=c("species", "sites"))


## 3) Non-metric Multidimensional Scaling
varespec.dist <- vegdist(varespec, method = "bray")
nms2 <- metaMDS(varespec.dist, k=2, try = 20)
nms3 <- metaMDS(varespec.dist, k=3, try = 20)
stressplot(nms3)


## 4) passively project environmental variables with envfit
envfit2 <- envfit(nms3, varechem)
plot(nms2, display="sites")
plot(envfit2)


## 5) Plot results of NMDS and color sites based on km classification
plot(nms2, display="sites", type="n")
points(nms2$points, col=varespec.km$cluster)


## [Advanced] Group colors in a PCA
varespec.km <- kmeans(varechem, centers=5) #Clustering on species!
ordiplot(vare.pca.st, display="sites", type="n")
points(vare.pca.st, display="sites", col=varespec.km$cluster)
points(vare.pca.st, display="species", pch="+")
text(vare.pca.st, display="species", labels = colnames(vare.pca.st$Ybar))




## ---- eval=T-------------------------------------------------------------------------------------------
#Hints
?varpart
?rda
?cca
?plot.varpart


## ---- eval=T, echo=F-----------------------------------------------------------------------------------
# RDA - Redundance Analysis ####
varespec.hel <- decostand(varespec, method="hellinger")
# Basic RDA with continuous variables using the "varespec" dataset
vare.rda <- rda(varespec.hel ~ pH + P + N, varechem)
vare.rda # Note the info for 3 constrained and 20 unconstrained axes
plot(vare.rda) # compare the interpretation with the nmds plot with envfit

# RDA with factors with the "dune" dataset
dune.hel <- decostand(dune, method="hellinger")
dune.rda <- rda(dune.hel ~., data = dune.env)
dune.rda # "Management" has four levels, CCA computes three contrasts (dummy variables)
plot(dune.rda, scaling=1, type="points")


# Model building: include/select what variables to include
# Model with all variables (not recommended)
mod1 <- rda(varespec.hel~., varechem) # "." uses all available variables
mod1 # This model uses all variables (Overfitted!!!)

# Model selection
mod0 <- rda(varespec.hel ~1, varechem) # creates a basis model equal to the
                                       # unconstrained version
mod <- step(mod0, scope=formula(mod1), test="perm") # Uses AIC for including variables
mod # shows the best model according to the AIC (these are supposed to be
    # the best variables)


# CCA - Canonical Correspondence Analysis ####
# Basically, the same as explained for RDA (model selection, etc.)
# Here just see the reduced model with four variables, as before
vare.cca <- cca(varespec.hel ~ Al + K + N, varechem)
vare.cca # You will see the results differ a little bit (but this depends on the data)
summary(vare.cca) # reports all results of the model

# Variance partitioning ####

# This function is useful to separate (partition) the effect of variables
# It is implemented only for RDA as a linear method (nor for CCA)

varp <- varpart (varespec.hel, ~ Al, ~ K, ~ N, data = varechem)
varp # Check for the marginal effects of each variable
plot (varp, digits = 2) # Creates a Venn's diagram. Digits define the number of decimals


