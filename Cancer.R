library(readxl)
CancerData <- read_excel("~/Desktop/Clustering Cancer/CancerData.xlsx", 
                         sheet = "Final")
View(CancerData)

Cancer <- CancerData
Cancer.num <- Cancer[2:13]

names(Cancer) <- c("Country", "Incidence", "Mortality", "SR_breast", "SR_lung", "SR_colon",
                       "Alcohol", "Tobacco", "HPV",
                       "Sc_colorectal", "Sc_mammography", "Equipment", "Drug")
names(Cancer.num) <- c("Incidence", "Mortality", "SR_breast", "SR_lung", "SR_colon",
                       "Alcohol", "Tobacco", "HPV",
                       "Sc_colorectal", "Sc_mammography", "Equipment", "Drug")
rownames(Cancer.num) <- Cancer$Country

#standardized and weighted data 
Cancer.scaled <- scale(Cancer.num)
Cancer.scaled <- as.data.frame(Cancer.scaled)

weights <- c(1, 1, 0.3333, 0.3333, 0.3333, 1, 1, 1, 0.5, 0.5, 1, 1)

Cancer.weighted <- mapply(FUN = function(arg1, arg2) ({arg1 *arg2}), 
                          arg1 = Cancer.scaled, 
                          arg2 = weights, 
                          SIMPLIFY = FALSE)
Cancer.weighted <- as.data.frame(Cancer.weighted)
rownames(Cancer.weighted) <- Cancer$Country

#standardized and weighted
Cancer.scaled <- scale(Cancer.num)
Cancer.scaled <- as.data.frame(Cancer.scaled)

weights <- c(1, 1, 0.3333, 0.3333, 0.3333, 1, 1, 1, 0.5, 0.5, 1, 1)

Cancer.weighted <- mapply(FUN = function(arg1, arg2) ({arg1 *arg2}), 
                          arg1 = Cancer.scaled, 
                          arg2 = weights, 
                          SIMPLIFY = FALSE)
Cancer.weighted <- as.data.frame(Cancer.weighted)
rownames(Cancer.weighted) <- Cancer$Country

# Hierarchical clustering 
res.dist <- daisy(Cancer.num, metric = "gower", 
                  weights = c(1, 1, 0.3333, 0.3333, 0.3333, 1, 1, 1, 0.5, 0.5, 1, 1))
res.hc <- hclust(res.dist, method = "ward.D2") 
fviz_dend(res.hc, cex = 0.5)
res.coph <- cophenetic(res.hc)
cor(res.dist, res.coph)

grp <- cutree(res.hc, k = 3)
head(grp, n = 3)
rownames(Cancer.num)[grp == 1]
rownames(Cancer.num)[grp == 2]
rownames(Cancer.num)[grp == 3]


#define linkage methods
m <- c("average", "single", "complete", "ward")
names(m) <- c("average", "single", "complete", "ward")

#function to compute agglomerative coefficient
ac <- function(x) {
  agnes(Cancer.weighted, method = x)$ac
}

#calculate agglomerative coefficient for each clustering linkage method
sapply(m, ac)
#ward's minimum variance method produces the highest coefficient 
#use that as the method

#perform hierarchical clustering using Ward's minimum variance
clust <- agnes(Cancer.weighted, method = "ward")

#produce dendrogram
pltree(clust, cex = 0.6, hang = -1, main = "Dendrogram") 

#compute distance matrix
res.dist <- daisy(Cancer.num, metric = "gower", 
                  weights = c(1, 1, 0.3333, 0.3333, 0.3333, 1, 1, 1, 0.5, 0.5, 1, 1))

#perform hierarchical clustering using Ward's method
final_clust <- hclust(res.dist, method = "ward.D2")
#cut the dendrogram into 3 clusters
groups <- cutree(final_clust, k=3)

#find number of observations in each cluster
table(groups)
fviz_dend(final_clust, cex = 0.5)
res.coph <- cophenetic(final_clust)
cor(res.dist, res.coph)

#append cluster labels to original data
final_data <- cbind(Cancer.num, cluster = groups)
#display first six rows of final data
head(final_data)

#name of countries in each cluster
grp <- cutree(final_clust, k = 3)
head(grp, n = 3)
rownames(Cancer.num)[grp == 1]
rownames(Cancer.num)[grp == 2]
rownames(Cancer.num)[grp == 3]

#find mean values for each cluster
aggregate(final_data, by=list(cluster=final_data$cluster), mean, na.rm = TRUE)

#dendrogram
fviz_dend(final_clust, k = 3, # Cut in four groups
          cex = 0.5, # label size
          k_colors = c("#2E9FDF", "#E7B800", "#FC4E07"),
          color_labels_by_k = TRUE, # color labels by groups
          rect = TRUE # Add rectangle around groups
)

# boxplot by indicators 
ggplot(final_data, aes(x = cluster, y = Incidence, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Mortality, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = SR_breast, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = SR_lung, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = SR_colon, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Alcohol, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Tobacco, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = HPV, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Sc_colorectal, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Sc_mammography, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Equipment, fill = cluster, group = cluster)) + 
  geom_boxplot()
ggplot(final_data, aes(x = cluster, y = Drug, fill = cluster, group = cluster)) + 
  geom_boxplot()

# heatmap 
numeric  <- data.matrix(Cancer.weighted)
heatmap(numeric, scale = "none")

# For assessing clustering tendency
# statistical method - hopkins statistics
hopkins <- get_clust_tendency(Cancer.weighted, n = nrow(Cancer.weighted)-1, graph = FALSE)
hopkins$hopkins_stat
# visual method
# Random data generated from the data set
random_df <- apply(Cancer.weighted, 2,
                   function(x){runif(length(x), min(x), (max(x)))})
random_df <- as.data.frame(random_df)
random_df <- scale(random_df)
rownames(Cancer.weighted) <- Cancer$Country
fviz_dist(dist(Cancer.weighted), show_labels = FALSE)+
  labs(title = "Cancer data")
fviz_dist(dist(random_df), show_labels = FALSE)+
  labs(title = "Random data")

# Compute clValid
clmethods <- c("hierarchical","pam")
country <- clValid(Cancer.weighted, nClust = 3:5,
                   clMethods = clmethods, validation = "internal")
# Summary
summary(country)

# Optimal number of clusters
# Elbow method
# function to compute total within-cluster sum of squares
fviz_nbclust(Cancer.weighted, pam, method = "wss", k.max = 5) + 
  theme_minimal() + ggtitle("K-medoids Elbow Method")

fviz_nbclust(Cancer.weighted, hcut, method = "wss", k.max = 5) + 
  theme_minimal() + ggtitle("Hierarchical Elbow Method")

# Silhouette
fviz_nbclust(Cancer.weighted, pam, method = "silhouette", k.max = 5) + 
  theme_minimal() + ggtitle("K-medoids Silhouette Plot")

fviz_nbclust(Cancer.weighted, hcut, method = "silhouette", k.max = 5) + 
  theme_minimal() + ggtitle("Hierarchical Silhouette Plot")

# Visualization of dendrograms (hierarchical)
hc <- hclust(res.dist, method = "ward.D2")
table(grp)
plot(hc, cex = 0.6)
rect.hclust(hc, k = 3, border = 2:5) #k=3
rect.hclust(hc, k = 4, border = 2:5) #k=4
rect.hclust(hc, k = 5, border = 2:5) #k=5

# Distribution of each indicator
# ase("indicator name")) 
ggplot(Cancer, aes(Drug)) +
  geom_histogram(color = "#000000", fill = "#0099F8") +
  ggtitle("Variable distribution") +
  theme_classic() +
  theme(plot.title = element_text(size = 18))

# visualizing missing values (missing value map)
vis_miss(Cancer)

