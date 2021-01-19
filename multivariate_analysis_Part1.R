## ---- eval=T-------------------------------------------------------------------------------------------
# download packages, if necessary
install.packages(c("vegan", "MASS", "cluster", "tree", 
                   "BiodiversityR", "gclus", "ecodist", "FD", "psych", "pheatmap"))


## ---- warning=F, message=F-----------------------------------------------------------------------------
#load packages
library(vegan)
library(MASS)
library(cluster)
library(tree)
library(gclus)
#library(ecodist)
library(FD)
library(psych) #for pairs-plot
library(pheatmap) #for heatmaps


## ---- eval=T-------------------------------------------------------------------------------------------
# Take a look at vegan's vignettes
browseVignettes("vegan")


## ---- eval=T-------------------------------------------------------------------------------------------
data(varespec) 
data(varechem)
?varespec




## ---- eval=T-------------------------------------------------------------------------------------------
data(dune) 
data(dune.env)
?dune


## ---- eval=T-------------------------------------------------------------------------------------------
# tips
?str
?class
?summary
?head
?tail
?dim
?nrows
?ncols
?rownames
?colnames
?range
?apply  #e.g. apply(varespec, MARGIN=1, "max")


## ---- eval=T-------------------------------------------------------------------------------------------
#Hint # Get help on the decostand() function
?decostand


## ---- eval=T, echo=T-----------------------------------------------------------------------------------
# Transformation and standardization of the species data 
## Simple transformations

# Partial view of the raw data (abundance codes)
varespec[1:5, 2:4]

## 1) Transform abundances to presence-absence (1-0)
varespec.pa <- decostand(varespec, method = "pa")
varespec.pa[1:5, 2:4]

## 2) Standardization by columns (species)
# Scale abundances by dividing them by the maximum value of each 
# species
# Note: MARGIN = 2 (column, default value) for argument "max"
varespec.scal <- decostand(varespec, "max")
varespec.scal[1:5, 2:4]
# Display the maximum in each transformed column
apply(varespec.scal, 2, max)

## 3) Standardization by rows (sites)
# Scale abundances by dividing them by the site totals
# (profiles of relative abundance by site)
varespec.rel <- decostand(varespec, "total") # default MARGIN = 1
varespec.rel[1:5, 2:4]
# Display the sum of row vectors to determine if the scaling worked
# properly
rowSums(varespec.rel) # equivalent to: apply(varespec.rel, 1, sum)

## 4) Standardization to zero mean and unit s.d.
varechem.st <- decostand(varechem, "standardize")
# verify it worked
colMeans(varechem.st)
apply(varechem.st, MARGIN=2, sd)


## ---- eval=T-------------------------------------------------------------------------------------------
#Hint # Get help on the following functions in vegan
?diversity
?specaccum
?rarefy
?rarecurve


## ---- echo=T, eval=T-----------------------------------------------------------------------------------
# Species richness
specnumber(BCI) # returns the species richness per plot
summary(specnumber(BCI)) # report descriptive statistics
hist(specnumber(BCI)) # have a look to the frequencies

# Diversity indices
diversity(BCI, index = "shannon") # Computes shannon diversity index per plot
hist(diversity(BCI, index = "shannon")) # see the histogram
diversity(BCI, index = "simpson") # Computes shannon diversity index per plot
hist(diversity(BCI, index = "simpson")) # see the histogram

# rarefaction (individual-based rarefaction)
rarefy(BCI, 20) # gives you the species per 20 individuals
rarecurve(BCI) # sample size reflects individuals

# rarefaction (sample-based)
spa <- specaccum(BCI) 
plot(spa) # Plot the rarefaction curve
plot(spa, ci.type="poly", col="blue", lwd=2, ci.lty=0, ci.col="lightblue") # just nicer




## ---- eval=T-------------------------------------------------------------------------------------------
#Hints # 
?vegdist
?as.matrix
?vegan::mantel #!! the FD package masks the mantel function from vegan!


## ---- echo=T, eval=T-----------------------------------------------------------------------------------
## 1) 
distance.spe.q1 <- vegdist(varespec)
distance.spe.q2 <- vegdist(varespec, method ="euclidean") #does not account for double 0s!
distance.spe.q3 <- vegdist(varespec, binary=FALSE) 
# binary = FALSE looks at the abundance; 
distance.spe.q4 <- vegdist(varespec, binary=TRUE) 
# TRUE looks at presence-absence (Sorenson's index)
# equivalent to distancce.spe.q1

## 2)
dim(distance.spe.q1)
length(distance.spe.q1) #How do you obtain this number?
distances <- as.matrix(distance.spe.q1)
dim(distances)

## 3.1) Q mode analysis - dissimilarity between sites based on species data
distance.spe.q1 <- vegdist(varespec) # "bray" is the default, not necessary to type

## 3.2) Q mode analysis - dissimilarity between sites based on environmental variables
distance.env.q1 <- vegdist(varechem, "euclidean") 
# even better would be to standardize the variables
varechem.st <- decostand(varechem, method = "standardize", MARGIN=2) #by column!
distance.env.q1 <- vegdist(varechem.st, "euclidean") 

## 3.3) R mode analysis - dissimilarity between species, 
## based on their presence absence in a site
## Need to Transpose - in R variables are conventionally in columns!
varespec.t <- t(varespec)
#transpose matrix of species abundancces
distance.spe.r1 <- dist(decostand(varespec, "chi.square")) 
# dist = euclidean distance in {base}

## 3.4) R mode analysis -  between environmental predictors
distance.env.r1 <- vegdist(t(varechem.st), method ="euclidean") 
# To compare env. variables it does not really make sense to calculate dissimilarities
# (even if technically possible). Tather one should calculate correlations
psych::pairs.panels(varechem[,1:5],
             method = "pearson", # correlation method
             hist.col = "#00AFBB",
             density = TRUE,  # show density plots
             ellipses = TRUE # show correlation ellipses
             )

## 4) Compare different dissimilarity matrices graphically
# dissimilarities based on species (Q-mode), with different metrics
pheatmap::pheatmap(distance.spe.q1, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(distance.spe.q2, cluster_rows = F, cluster_cols = F)
pheatmap::pheatmap(distance.spe.q4, cluster_rows = F, cluster_cols = F)

# The Mantel test ####
# Mantel test calculates correlations between dissimilarities
# The higher is the result, the more similar the groups compared
vegan::mantel(distance.spe.q1, distance.spe.q2)

## 5) find the most similar plots 
#transform to matrix, and replace diagonal with ones 
#[to avoid considering distances of plots to themselves]
totest <- as.matrix(distance.spe.q1) + diag(nrow=nrow(varespec)) 
minn <- which(totest==min(totest), arr.ind=T)
totest[minn]


## ---- eval=T, echo=T-----------------------------------------------------------------------------------
# my solution
myjac <- function(x,y){ 
  x.pa <- rbind(x,y)>0
  return(1 - sum(colSums(x.pa)==2) / (sum(colSums(x.pa)==2) + sum(colSums(x.pa)==1)))
  }


