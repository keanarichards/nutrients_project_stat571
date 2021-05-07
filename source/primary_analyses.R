## ----include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
options(scipen = 0, digits = 3, tibble.print_max = 50) 

knitr::opts_chunk$set(echo = FALSE, message=FALSE, warning = F)

if(!require("pacman")) install.packages("pacman")
pacman::p_load(tidyverse, here, glmnet, car, data.table, summarytools, corrplot, GGally, varhandle, gtsummary, pROC, stargazer, sjPlot, report,tm,SnowballC,wordcloud,RColorBrewer, imputeTS,  clValid, cluster, factoextra, fpc) 


plot_aes = theme_minimal() +
  theme(legend.position = "top",
        text = element_text(size = 15, family = "Futura Medium"),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.ticks.y = element_blank())

# load data ---------------------------------------------------------------

food_name <- read.csv(here("data", "FOOD_NAME.csv"))
food_group <- read.csv(here("data", "FOOD_GROUP.csv"))
nutrient_amount <- read.csv(here("data", "NUTRIENT_AMOUNT.csv"))
nutrient_name <- read.csv(here("data", "NUTRIENT_NAME.csv"))



## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## creating wide version of data for PCA

## merging datasets on identifiers

## start by identifying nutrients - merge nutrient amount & nutrient name 

nutrient_amount_merged <- left_join(nutrient_amount %>% select(NutrientID, NutrientValue, FoodID), nutrient_name %>% select(NutrientID, NutrientName, NutrientUnit), by = "NutrientID") 

## then merging add food group labels to food name dataset

food_name_merged <- left_join(food_name %>% select(FoodID, FoodGroupID, FoodDescription), food_group %>% select(FoodGroupID, FoodGroupName), by = "FoodGroupID") 

## creating wide version of nutrient_amount dataset to be merged with food_name

nutrient_amount_wide <- nutrient_amount_merged %>% select(-c(NutrientID, NutrientUnit)) %>%  pivot_wider(names_from = NutrientName, values_from = NutrientValue)

## merging wide nutrient_amount with food_name_merged

wide_data <- left_join(food_name_merged, nutrient_amount_wide, by = "FoodID") 



## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## creating long version of data for logistic regression

long_data <- wide_data %>% pivot_longer(cols = PROTEIN:`OXALIC ACID`, names_to = "NutrientName", values_to = "NutrientValue")



## ----dealing with NAs, include=FALSE-------------------------------------------------------------------------------------------------------------------------
######## data with mean imputation
data_mean_imp_temp<-wide_data %>% select(-c(FoodID,FoodGroupID,FoodDescription,FoodGroupName))
data_mean_imp<-imputeTS::na_mean(data_mean_imp_temp,option="mean")
data_mean_imp_standardized<-data.frame(scale(data_mean_imp,center=T, scale=T))

data_mean_imp_standardized<- data_mean_imp_standardized[, -c(146, 152)]



# head(data_mean_imp)
# NA inspection per column: sapply(wide_data, function(x) sum(is.na(x)))

# sum(is.na(data_mean_imp))


## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
# optimal number of clusters : according to silhouette method, 8 clusters is optimal  
set.seed(0)
silh_ID1 = factoextra::fviz_nbclust(data_mean_imp_standardized, kmeans, method = "silhouette")

##### run k-means clustering
data_mean_imp.kmeans <- data_mean_imp_standardized %>% kmeans(centers = 8)  # centers: number of cluster
print(c('Size of each cluster is', data_mean_imp.kmeans$size))



## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
centers=as.data.frame(data_mean_imp.kmeans$centers)

##### randomly choose two variables and plot clustering results (hard coded)


clusterplot1 = data.frame(Carb = data_mean_imp_standardized$CARBOHYDRATE..TOTAL..BY.DIFFERENCE., 
          Kcal =data_mean_imp_standardized$ENERGY..KILOCALORIES.,
          group = as.factor(data_mean_imp.kmeans$cluster)) %>% 
  ggplot()+
  geom_point(aes(x=Carb,y=Kcal,col=group))+theme_bw()+labs(color="Cluster") +ggtitle("Clustering over randomly chosen two variables")+xlab("Carbohydrate")+ylab("Energy(Kcal)")+
	geom_point(data=centers,aes(x=`CARBOHYDRATE..TOTAL..BY.DIFFERENCE.`,y=`ENERGY..KILOCALORIES.`,cex=20),show.legend=F)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[1]+.1,y=centers$ENERGY..KILOCALORIES.[1],label = paste0("Center", 1), hjust = "left", fontface = "bold", size = 4)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[2]+.1,y=centers$ENERGY..KILOCALORIES.[2],label = paste0("Center", 2), hjust = "left", fontface = "bold", size = 4)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[3:4]+.1,y=centers$ENERGY..KILOCALORIES.[3:4],label = paste0("Center", 3:4), hjust = "left", fontface = "bold", size = 4)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[6:7]+.1,y=centers$ENERGY..KILOCALORIES.[6:7],label = paste0("Center", 6:7), hjust = "left", fontface = "bold", size = 4)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[5]+.1,y=centers$ENERGY..KILOCALORIES.[5],label = paste0("Center", 5), hjust = "left", fontface = "bold", size = 4)+ 
	ggplot2::annotate(geom = "text",x=centers$CARBOHYDRATE..TOTAL..BY.DIFFERENCE.[8]+.1,y=centers$ENERGY..KILOCALORIES.[8],label = paste0("Center", 8), hjust = "left", fontface = "bold", size = 4)

# + ggrepel::geom_text_repel(aes(label = team))



## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
###### data table with food name + clustering result
result_regular_kmeans<-wide_data %>% select (FoodDescription) %>%
        mutate(group = data_mean_imp.kmeans$cluster) %>% arrange(group)

###### Cluster 1 
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==1]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with","dehydrated")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc1.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()
###### Cluster 2 
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==2]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc2.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()


###### Cluster 3 
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==3]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc3.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()


###### Cluster 4 
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==4]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc4.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()



###### Cluster 5
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==5]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc5.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()


###### Cluster 6
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==6]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)

dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc6.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()



###### Cluster 7
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==7]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc7.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()



###### Cluster 8
docs<-Corpus(VectorSource(result_regular_kmeans$FoodDescription[result_regular_kmeans$group==8]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc8.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()

# drained, dehydrated, boiled, raw, frozen, prepared, ready, unprepared, made



## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
###### average calory count per cluster
avg_calories=data_mean_imp_standardized %>% mutate(group = data_mean_imp.kmeans$cluster) %>% group_by(group) %>% summarize(mean.calories=mean(ENERGY..KILOCALORIES.),se.calories=sd(ENERGY..KILOCALORIES.) /sqrt(length(ENERGY..KILOCALORIES.)),std.calories=sd(ENERGY..KILOCALORIES.))

kcal_cluster1 = ggplot(avg_calories) +
    geom_bar(aes(x=group, y=mean.calories), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean.calories-se.calories, ymax=mean.calories+se.calories), width=0.4, colour="orange", alpha=0.9, size=1.3)+xlab("Cluster")+ylab("Average Calories (Kcal)") + scale_x_continuous(breaks = seq(1,8, by = 1))+ggtitle("Average Kcal of Foods per Cluster")+ 
  theme(plot.title = element_text(hjust = 0.5))


## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
# PCA of nutrients 
data_mean_imp_temp =data_mean_imp_standardized %>% mutate(group = data_mean_imp.kmeans$cluster) %>% select(-c("FATTY ACIDS, MONOUNSATURATED, 12:1, LAUROLEIC","NA"))

data_mean_imp.c1<-data_mean_imp_temp[data_mean_imp_temp$group==1,]
data_mean_imp.c2<-data_mean_imp_temp[data_mean_imp_temp$group==2,]
data_mean_imp.c3<-data_mean_imp_temp[data_mean_imp_temp$group==3,]
data_mean_imp.c4<-data_mean_imp_temp[data_mean_imp_temp$group==4,]
data_mean_imp.c5<-data_mean_imp_temp[data_mean_imp_temp$group==5,]
data_mean_imp.c6<-data_mean_imp_temp[data_mean_imp_temp$group==6,]
data_mean_imp.c7<-data_mean_imp_temp[data_mean_imp_temp$group==7,]
data_mean_imp.c8<-data_mean_imp_temp[data_mean_imp_temp$group==8,]

# Cluster 1 
sdzero=data.frame(apply(data_mean_imp.c1,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c1=data_mean_imp.c1%>%select(-c(names(data_mean_imp.c1)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c1,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")

# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients1 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 1") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 2
sdzero=data.frame(apply(data_mean_imp.c2,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c2=data_mean_imp.c2%>%select(-c(names(data_mean_imp.c2)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c2,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients2 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 2") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))



# Cluster 3
sdzero=data.frame(apply(data_mean_imp.c3,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c3=data_mean_imp.c3%>%select(-c(names(data_mean_imp.c3)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c3,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients3 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 3") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))



# Cluster 4
sdzero=data.frame(apply(data_mean_imp.c4,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c4=data_mean_imp.c4%>%select(-c(names(data_mean_imp.c4)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c4,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients4 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 4") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 5
sdzero=data.frame(apply(data_mean_imp.c5,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c5=data_mean_imp.c5%>%select(-c(names(data_mean_imp.c5)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c5,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients5 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 5") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 6
sdzero=data.frame(apply(data_mean_imp.c6,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c6=data_mean_imp.c6%>%select(-c(names(data_mean_imp.c6)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c6,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients6 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 6") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 7
sdzero=data.frame(apply(data_mean_imp.c7,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c7=data_mean_imp.c7%>%select(-c(names(data_mean_imp.c7)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c7,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients7 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 7") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 8
sdzero=data.frame(apply(data_mean_imp.c8,2,sd))
whichone=which(sdzero == 0)
data_mean_imp.c8=data_mean_imp.c8%>%select(-c(names(data_mean_imp.c8)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c8,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =  abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients8 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 8") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = element_text(hjust = 0.5))


## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
## RUN PCA of nutrient data

pca_unscaled <- prcomp(data_mean_imp,center=TRUE,scale=FALSE)

# names(pca)
pca_unscaled.loading<-pca_unscaled$rotation
# knitr::kable(pca.loading)
# summary(pca)$importance 

# How many PCs to use? --> 6~7 PCs seem to explain most variance according to elbow rule
plot(summary(pca_unscaled)$importance[2, ],  # PVE
     ylab="PVE",
     xlab="Number of PC's",
     pch = 16, 
     main="Scree Plot of PVE for Nutrients_UNSCALEDdata")
pve_unscaled <- summary(pca_unscaled)$importance[2, 1:10]

plot(summary(pca_unscaled)$importance[3, ], pch=16,
     ylab="Cumulative PVE",
     xlab="Number of PC's",
     main="Scree Plot of Cumulative PVE for Nutrients_UNSCALED")

# plot(pve_unscaled, type="b", pch = 19, frame = FALSE)
# plot(summary(pca)$importance[3, ], pch=16,
 #    ylab="Cumulative PVE",
   #  xlab="Number of PC's",
  #   main="Scree Plot of Cumulative PVE for AFQT")


## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
## RUN PCA of nutrient data
## In order to scale the data, needed to remove variables with 0 standard deviaion (to solve this error: Error in prcomp.default(data_mean_imp, center = TRUE, scale. = TRUE) : cannot rescale a constant/zero column to unit variance)
##  "FATTY ACIDS, MONOUNSATURATED, 12:1, LAUROLEIC" and "NA" were removed for this reason

#te=data.frame(apply(data_mean_imp,2,sd))
#which(te == 0)
#names(data_mean_imp)[146]
#names(data_mean_imp)[152]

data_mean_imp_scale=data_mean_imp %>% select(-c("FATTY ACIDS, MONOUNSATURATED, 12:1, LAUROLEIC","NA"))

# data_mean_imp_scale<-scale(as.matrix(data_mean_imp_scale),center=T, scale=T)

pca_scaled <- prcomp(data_mean_imp_scale,center=TRUE,scale=TRUE)
# pca_scaled <- prcomp(data_mean_imp_scale,center=TRUE,scale.=TRUE)

# names(pca)
pca_scaled.loading<-pca_scaled $rotation
# knitr::kable(pca.loading)
# summary(pca)$importance 

# How many PCs to use? -
png(here("figs", "scree_p1.png"), width=12, height=8, units="in", res=300)
plot(summary(pca_scaled )$importance[2, ],  # PVE
     ylab="PVE",
     xlab="Number of PC's",
     pch = 16, 
     main="Scree Plot of PVE for Nutrients_SCALED_data")
dev.off()

pve_scaled  <- summary(pca_scaled)$importance[2, 1:10]
plot(pca_scaled)
plot(summary(pca_scaled)$importance[3, ], pch=16,
     ylab="Cumulative PVE",
     xlab="Number of PC's",
     main="Scree Plot of Cumulative PVE for Nutrients_SCALED")
# sanity check : round(cov(pca_scaled$x), 4) ; round(t(pca_scaled$rotation) %*% pca_scaled$rotation)


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# PC1 and PC2
pca <- prcomp(data_mean_imp_scale,center=TRUE,scale.=FALSE)

# kmeans 
data.mean.imp.spectrum.kmeans<- kmeans(x = pca$x[,1:10], 8)

p_unscaled <- data.table(x = pca$x[,1], 
                y = pca$x[,2],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster)) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC1 scores") + ylab("PC2 scores")+ggtitle("UNSCALED DATA") + theme(plot.title =element_text(hjust = 0.5))

# PC1 and PC2
pca <- prcomp(data_mean_imp_scale,center=TRUE,scale.=TRUE)

# kmeans 
data.mean.imp.spectrum.kmeans<- kmeans(x = pca$x[,1:7], 8)

p_scaled <- data.table(x = pca$x[,1], 
                y = pca$x[,2],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster)) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC1 scores") + ylab("PC2 scores")+ggtitle("SCALED DATA") + theme(plot.title =element_text(hjust = 0.5))

p_unscaled
p_scaled



## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------

stdzero=data.frame(apply(data_mean_imp,2,sd))
whichone=which(stdzero== 0)
data_mean_imp_scale=data_mean_imp %>% select(-c(names(data_mean_imp)[whichone]))
pca <- prcomp(data_mean_imp_scale,center=TRUE,scale=TRUE)


## plot top 20 loadings
top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(pca$rotation[,1]),
                  gene = rownames(pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(pca$rotation[,2]),
                  gene = rownames(pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top loadings") +
  xlab("Nutrient") + ylab("Absolute values of PC loadings")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
# optimal number of clusters : according to silhouette method, it's 8.
set.seed(0)
silh_ID3 = factoextra::fviz_nbclust(data_mean_imp_scale, kmeans, method = "silhouette")

# kmeans 
data.mean.imp.spectrum.kmeans<- kmeans(x = pca$x[,1:10], 8)



## ---- include=F----------------------------------------------------------------------------------------------------------------------------------------------
# Get row numbers for each cluster 
temp1<-data_mean_imp_scale %>% mutate(cluster=data.mean.imp.spectrum.kmeans$cluster)
rownum1=as.numeric(rownames(temp1[temp1$cluster==1,]))
rownum2=as.numeric(rownames(temp1[temp1$cluster==2,]))
rownum3=as.numeric(rownames(temp1[temp1$cluster==3,]))
rownum4=as.numeric(rownames(temp1[temp1$cluster==4,]))
rownum5=as.numeric(rownames(temp1[temp1$cluster==5,]))
rownum6=as.numeric(rownames(temp1[temp1$cluster==6,]))
rownum7=as.numeric(rownames(temp1[temp1$cluster==7,]))
rownum8=as.numeric(rownames(temp1[temp1$cluster==8,]))



## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
# PC1 and PC2, label= one representative food with highest PC score
max_pcscore_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,1]))],
    rownum2[which.max(abs(pca$x[rownum2,1]))],
    rownum3[which.max(abs(pca$x[rownum3,1]))],
    rownum4[which.max(abs(pca$x[rownum4,1]))],
    rownum5[which.max(abs(pca$x[rownum5,1]))],
    rownum6[which.max(abs(pca$x[rownum6,1]))],
    rownum7[which.max(abs(pca$x[rownum7,1]))],
    rownum8[which.max(abs(pca$x[rownum8,1]))])

max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcscore_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_scale)[1])
for (i in 1:length(max_pcscore_rownum)){max_pc_food[max_pcscore_rownum[i]]=as.character(max_pc_score_df[i])}


temp2<-data_mean_imp_scale%>% mutate(max_pc_food=max_pc_food)

clusterplot3 <- data.table(x = pca$x[,1], 
                y = pca$x[,2],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8))+theme_bw() + labs(color = "Food Clusters") +xlab("PC1 scores") + ylab("PC2 scores") + ggrepel::geom_text_repel(aes(x = pca$x[,1],y = pca$x[,2], label=name, color = col)) +ggtitle("8 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest absolute PC1 score") + theme(plot.title =element_text(hjust = 0.5))

#### PC3 and PC4, label= one representative food with highest absolute PC3 score
max_pcscore_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,3]))],
    rownum2[which.max(abs(pca$x[rownum2,3]))],
    rownum3[which.max(abs(pca$x[rownum3,3]))],
    rownum4[which.max(abs(pca$x[rownum4,3]))],
    rownum5[which.max(abs(pca$x[rownum5,3]))],
    rownum6[which.max(abs(pca$x[rownum6,3]))],
    rownum7[which.max(abs(pca$x[rownum7,3]))],
    rownum8[which.max(abs(pca$x[rownum8,3]))])

max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcscore_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_scale)[1])
for (i in 1:length(max_pcscore_rownum)){max_pc_food[max_pcscore_rownum[i]]=as.character(max_pc_score_df[i])}


temp2<-data_mean_imp_scale%>% mutate(max_pc_food=max_pc_food)

p1 <- data.table(x = pca$x[,3], 
                y = pca$x[,4],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8))+theme_bw() + labs(color = "Food Clusters") +xlab("PC3 scores") + ylab("PC4 scores") + ggrepel::geom_text_repel(aes(x = pca$x[,3],y = pca$x[,4], label=name, color = col)) +ggtitle("8 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest absolute PC3 score") + theme(plot.title =element_text(hjust = 0.5))

#### PC1 and PC5, label= one representative food with highest absolute PC1 score
max_pcscore_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,1]))],
    rownum2[which.max(abs(pca$x[rownum2,1]))],
    rownum3[which.max(abs(pca$x[rownum3,1]))],
    rownum4[which.max(abs(pca$x[rownum4,1]))],
    rownum5[which.max(abs(pca$x[rownum5,1]))],
    rownum6[which.max(abs(pca$x[rownum6,1]))],
    rownum7[which.max(abs(pca$x[rownum7,1]))],
    rownum8[which.max(abs(pca$x[rownum8,1]))])

max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcscore_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_scale)[1])
for (i in 1:length(max_pcscore_rownum)){max_pc_food[max_pcscore_rownum[i]]=as.character(max_pc_score_df[i])}

temp2<-data_mean_imp_scale%>% mutate(max_pc_food=max_pc_food)

p2 <- data.table(x = pca$x[,1], 
                y = pca$x[,5],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC1 scores") + ylab("PC5 scores") + ggrepel::geom_text_repel(aes(x = pca$x[,1],y = pca$x[,5], label=name, color = col, )) +ggtitle("8 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest PC1 score") + theme(plot.title =element_text(hjust = 0.5))

p1
p2



## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
###### data table with food name + clustering result
result_mean_imp_spectrum_kmeans<-wide_data %>% select (FoodDescription) %>%
        mutate(group = data.mean.imp.spectrum.kmeans$cluster) %>% arrange(group)

###### Cluster 1 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==1]))

# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))

###### Cluster 2 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==2]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))


###### Cluster 3 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==3]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))


###### Cluster 4 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==4]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Set1"))



###### Cluster 5
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==5]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))


###### Cluster 6
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==6]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Set2"))


###### Cluster 7
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==7]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))



###### Cluster 8
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==8]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeWords,c("and","with")) %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))




## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
###### average calory count per cluster
avg_calories=data_mean_imp_scale %>% mutate(group=data.mean.imp.spectrum.kmeans$cluster) %>% group_by(group) %>% summarize(mean.calories=mean(`ENERGY (KILOCALORIES)`),se.calories=sd(`ENERGY (KILOCALORIES)`) /sqrt(length(`ENERGY (KILOCALORIES)`)),std.calories=sd(`ENERGY (KILOCALORIES)`))

ggplot(avg_calories) +
    geom_bar(aes(x=group, y=mean.calories), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean.calories-se.calories, ymax=mean.calories+se.calories), width=0.4, colour="orange", alpha=0.9, size=1.3)+xlab("Cluster")+ylab("Average Calories (Kcal)") + scale_x_continuous(breaks = seq(1,8, by = 1))+ggtitle("Average Kcal of Foods per Cluster")+ 
  theme(plot.title = element_text(hjust = 0.5))




## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
# PCA of nutrients 
data_mean_imp_temp1 =data_mean_imp_scale %>% mutate(group = data.mean.imp.spectrum.kmeans$cluster)
data_mean_imp.c1<-data_mean_imp_temp1[data_mean_imp_temp1$group==1,]%>% select(c(-group))
data_mean_imp.c2<-data_mean_imp_temp1[data_mean_imp_temp1$group==2,]%>% select(c(-group))
data_mean_imp.c3<-data_mean_imp_temp1[data_mean_imp_temp1$group==3,]%>% select(c(-group))
data_mean_imp.c4<-data_mean_imp_temp1[data_mean_imp_temp1$group==4,]%>% select(c(-group))
data_mean_imp.c5<-data_mean_imp_temp1[data_mean_imp_temp1$group==5,]%>% select(c(-group))
data_mean_imp.c6<-data_mean_imp_temp1[data_mean_imp_temp1$group==6,]%>% select(c(-group))
data_mean_imp.c7<-data_mean_imp_temp1[data_mean_imp_temp1$group==7,]%>% select(c(-group))
data_mean_imp.c8<-data_mean_imp_temp1[data_mean_imp_temp1$group==8,]%>% select(c(-group))

# Cluster 1 

sdzero=data.frame(apply(data_mean_imp.c1,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c1=data_mean_imp.c1%>%select(-c(names(data_mean_imp.c1)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c1,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading =abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 1") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 2
sdzero=data.frame(apply(data_mean_imp.c2,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c2=data_mean_imp.c2%>%select(-c(names(data_mean_imp.c2)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c2,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 2") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))



# Cluster 3
sdzero=data.frame(apply(data_mean_imp.c3,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c3=data_mean_imp.c3%>%select(-c(names(data_mean_imp.c3)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c3,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 3") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))



# Cluster 4

sdzero=data.frame(apply(data_mean_imp.c4,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c4=data_mean_imp.c4%>%select(-c(names(data_mean_imp.c4)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c4,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading =abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 4") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 5
sdzero=data.frame(apply(data_mean_imp.c5,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c5=data_mean_imp.c5%>%select(-c(names(data_mean_imp.c5)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c5,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 5") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 6
sdzero=data.frame(apply(data_mean_imp.c6,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c6=data_mean_imp.c6%>%select(-c(names(data_mean_imp.c6)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c6,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 6") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 7
sdzero=data.frame(apply(data_mean_imp.c7,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c7=data_mean_imp.c7%>%select(-c(names(data_mean_imp.c7)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c7,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 7") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))


# Cluster 8
sdzero=data.frame(apply(data_mean_imp.c8,2,sd)) # removing 0 variance variables to run PCA
whichone=which(sdzero == 0)
data_mean_imp.c8=data_mean_imp.c8%>%select(-c(names(data_mean_imp.c8)[whichone]))
data_mean_imp.pca <- prcomp(data_mean_imp.c8,center=T,scale=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,1]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(data_mean_imp.pca$rotation[,2]),
                  gene = rownames(data_mean_imp.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 8") +
  xlab("Nutrient") + ylab("Absolute value of PC loading")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))






## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## mean silhouette for reg k means clusterse
ss <- cluster::silhouette(data_mean_imp.kmeans$cluster, dist(data_mean_imp_standardized))
# Statistics for k-means clustering
km_stats <- cluster.stats(dist(data_mean_imp_standardized),  data_mean_imp.kmeans$cluster)

mean_ss <- mean(ss[, 3])
plot_silhouette <- fviz_silhouette(ss)

## mean silhouette for spectrum clustering
ss1 <- cluster::silhouette(data.mean.imp.spectrum.kmeans$cluster, dist(pca$x[,1:10]))
km_stats1 <- cluster.stats(dist(data_mean_imp_standardized), data.mean.imp.spectrum.kmeans$cluster)

mean_ss1 <- mean(ss1[, 3])

plot_silhouette1 <- fviz_silhouette(ss1)

## ----include=FALSE-------------------------------------------------------------------------------------------------------------------------------------------
# create nutrient data for Mark

scenario=data.frame(Nutrient=names(data_mean_imp_standardized),Average_Amount=colMeans(data_mean_imp_standardized), STD=apply(data_mean_imp_standardized,2,sd)) 

scenario['Scenario1']=scenario$Average_Amount

scenario[which(scenario$Nutrient=="ENERGY (KILOCALORIES)"),]$Scenario1 <-scenario[which(scenario$Nutrient=="ENERGY (KILOCALORIES)"),]$Average_Amount*1.32


scenario[which(scenario$Nutrient=="CARBOHYDRATE, TOTAL (BY DIFFERENCE)"),]$Scenario1<-scenario[which(scenario$Nutrient=="CARBOHYDRATE, TOTAL (BY DIFFERENCE)"),]$Average_Amount*1.2

scenario[which(scenario$Nutrient=="PROTEIN"),]$Scenario1<-scenario[which(scenario$Nutrient=="PROTEIN"),]$Average_Amount*4.15

scenario[which(scenario$Nutrient=="FATTY ACIDS, SATURATED, TOTAL"),]$Scenario1<-scenario[which(scenario$Nutrient=="FATTY ACIDS, SATURATED, TOTAL"),]$Average_Amount*270

scenario[which(scenario$Nutrient=="FAT (TOTAL LIPIDS)"),]$Scenario1<-scenario[which(scenario$Nutrient=="FAT (TOTAL LIPIDS)"),]$Average_Amount*10.8



## ----include=F-----------------------------------------------------------------------------------------------------------------------------------------------
tempp=data.frame(t(scenario$Scenario1))
oldn=names(tempp)
newn=names(data_mean_imp_standardized)
s1_temp=tempp %>% rename_at(vars(oldn), ~ newn)

# data including scenario 1 
s1<-rbind(data_mean_imp_standardized,s1_temp)

set.seed(0)
# spectrum kmeans 
pca <- prcomp(s1,center=TRUE,scale.=TRUE)
data.mean.imp.spectrum.kmeans<- kmeans(x = pca$x[,1:10], 8)

data.mean.imp.spectrum.kmeans$cluster[5691] 

# regular kmeans
data.mean.imp.kmeans<- kmeans(x=s1, 8)
data.mean.imp.kmeans$cluster[5691] 

# Plot 
centers=as.data.frame(data.mean.imp.spectrum.kmeans$centers)

p <- data.table(x = pca$x[,1], 
                y = pca$x[,2],
                col = as.factor(data.mean.imp.spectrum.kmeans$cluster)) %>%
  ggplot() + geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(8))+theme_bw() + labs(color = "Food Clusters") +xlab("PC1 scores") + ylab("PC2 scores") +ggtitle("8 Food Clusters Based on 10 PCs") + theme(plot.title =element_text(hjust =0.5))+geom_point(data=centers,aes(x="PC1",y="PC2",cex=20),show.legend=F)+  ggplot2::annotate(geom = "text",x=centers$PC1,y=centers$PC2,label = paste0("Center", 8), hjust = "left", fontface = "bold", size = 4) 



## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# Get row numbers for each cluster 
temp1<-s1%>% mutate(cluster=data.mean.imp.kmeans$cluster)
clst_results<-data.mean.imp.kmeans$cluster[5691] 
clst_center<-data.frame(data.mean.imp.kmeans$centers[data.mean.imp.kmeans$cluster[5691],])

distance_from_centers=data.frame()

for (i in 1:dim(s1)[1]){distance_from_centers[i,1]=dist(rbind(clst_center[,1],data.frame(t(s1[i,]))[,1]),method="euclidean")}

foodname=data.frame(wide_data$FoodDescription)
foodname[5691,1]<-"NA"

distance_from_centers %>% mutate("Food"=foodname[,1]) %>% arrange(V1)




## ------------------------------------------------------------------------------------------------------------------------------------------------------------
# select crucial nutrients from mean imputed data
data_mean_imp_crucial=data_mean_imp %>% select(c(starts_with("VIT"),MAGNESIUM,CALCIUM,PHOSPHORUS,POTASSIUM,`CHOLINE, TOTAL`,IRON,SELENIUM,ZINC,MANGANESE,COPPER,MOISTURE,PROTEIN,`CARBOHYDRATE, TOTAL (BY DIFFERENCE)`,`FAT (TOTAL LIPIDS)`,`FATTY ACIDS, TRANS, TOTAL`,`FATTY ACIDS, SATURATED, TOTAL`,`FATTY ACIDS, MONOUNSATURATED, TOTAL`,`FATTY ACIDS, POLYUNSATURATED, TOTAL`,`FATTY ACIDS, TOTAL TRANS-MONOENOIC`,`FATTY ACIDS, TOTAL TRANS-POLYENOIC`,`ENERGY (KILOCALORIES)`))

oldnames=names(data_mean_imp_crucial%>%select(-c(contains("VIT"),contains("PROTEIN"),contains("CARBO"),contains("FAT"),contains("ENERGY"))))

newnames=c("Minerals_MAGNESIUM","Minerals_CALCIUM","Minerals_PHOSPHORUS","Minerals_POTASSIUM", "Minerals_CHOLINE, TOTAL" ,"Minerals_IRON","Minerals_SELENIUM","Minerals_ZINC","Minerals_MANGANESE","COPPER","Water_MOISTURE")

data_mean_imp_crucial=data_mean_imp_crucial %>% rename_at(vars(oldnames), ~ newnames)


## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
## RUN PCA of nutrient data
pca <- prcomp(data_mean_imp_crucial,center=TRUE,scale=TRUE)

pve <- summary(pca)$importance[2, 1:15]
plot(pve, type="b", pch = 19, frame = FALSE)
png(here("figs", "scree_p2.png"), width=12, height=8, units="in", res=300)
scree_p2 = plot(summary(pca)$importance[3, ], pch=16,
     ylab="Cumulative PVE",
     xlab="Number of PC's",
     main="Scree Plot of Cumulative PVE for Nutrients_SCALED")
dev.off()

# will use 10 PCs


## ------------------------------------------------------------------------------------------------------------------------------------------------------------
## plot top 20 loadings
top_k <- 20

## get pc1 and pc2
pc1 <- data.frame(loading = abs(pca$rotation[,1]),
                  gene = rownames(pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = abs(pca$rotation[,2]),
                  gene = rownames(pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top loadings") +
  xlab("Nutrient") + ylab("Absolute values of PC loadings")+
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1),plot.title = 
element_text(hjust = 0.5))




## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
# optimal number of clusters : according to silhouette method, it's 2
set.seed(0)
silh_ID2 = factoextra::fviz_nbclust(data_mean_imp_crucial, kmeans, method = "silhouette")

# kmeans 
data.crucial.mean.imp.spectrum.kmeans<- kmeans(x = pca$x[, 1:10], 2)


## ---- include=F----------------------------------------------------------------------------------------------------------------------------------------------
# Get row numbers for each cluster 
temp1<-data_mean_imp_crucial %>% mutate(cluster=data.crucial.mean.imp.spectrum.kmeans$cluster)
rownum1=as.numeric(rownames(temp1[temp1$cluster==1,]))
rownum2=as.numeric(rownames(temp1[temp1$cluster==2,]))


## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
# PC1 and PC2, label= one representative food with highest PC score
max_pcloading_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,1]))],
    rownum2[which.max(abs(pca$x[rownum2,1]))])


max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcloading_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_crucial)[1])
for (i in 1:length(max_pcloading_rownum)){max_pc_food[max_pcloading_rownum[i]]=as.character(max_pc_score_df[i])}

temp2<-data_mean_imp_crucial%>% mutate(max_pc_food=max_pc_food)

clusterplot2 <- data.table(x = pca$x[,1], 
                y = pca$x[,2],
                col = as.factor(data.crucial.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(2)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC1") + ylab("PC2") + ggrepel::geom_text_repel(aes(x = pca$x[,1],y = pca$x[,2], label=name, color = col, )) +ggtitle("2 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest PC1 score")+ theme(plot.title =element_text(hjust = 0.5))

#### PC3 and PC4, label= one representative food with highest PC3 score
max_pcloading_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,3]))],
    rownum2[which.max(abs(pca$x[rownum2,3]))])

max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcloading_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_crucial)[1])
for (i in 1:length(max_pcloading_rownum)){max_pc_food[max_pcloading_rownum[i]]=as.character(max_pc_score_df[i])}

temp2<-data_mean_imp_crucial%>% mutate(max_pc_food=max_pc_food)

p1 <- data.table(x = pca$x[,3], 
                y = pca$x[,4],
                col = as.factor(data.crucial.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(2)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC3") + ylab("PC4") + ggrepel::geom_text_repel(aes(x = pca$x[,3],y = pca$x[,4], label=name, color = col, )) +ggtitle("2 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest PC3 score")+ theme(plot.title =element_text(hjust = 0.5))

#### PC1 and PC6, label= one representative food with highest PC1 score
max_pcloading_rownum=
  c(rownum1[which.max(abs(pca$x[rownum1,1]))],
    rownum2[which.max(abs(pca$x[rownum2,1]))])

max_pc_score_df<-as.factor(wide_data$FoodDescription[max_pcloading_rownum])
max_pc_food <- rep(NA, dim(data_mean_imp_crucial)[1])
for (i in 1:length(max_pcloading_rownum)){max_pc_food[max_pcloading_rownum[i]]=as.character(max_pc_score_df[i])}

temp2<-data_mean_imp_crucial%>% mutate(max_pc_food=max_pc_food)

p2 <- data.table(x = pca$x[,1], 
                y = pca$x[,6],
                col = as.factor(data.crucial.mean.imp.spectrum.kmeans$cluster),name=temp2$max_pc_food) %>%
  ggplot() + 
  geom_point(aes(x = x, y = y, col = col))+
  scale_color_manual(values = scales::hue_pal()(2)) +
  theme_bw() + labs(color = "Food Clusters") +xlab("PC1") + ylab("PC6") + ggrepel::geom_text_repel(aes(x = pca$x[,1],y = pca$x[,6], label=name, color = col, )) +ggtitle("2 Food Clusters Based on 10 PCs",subtitle = "label=food per cluster with highest PC1 score")+ theme(plot.title =element_text(hjust = 0.5))

p1
p2



## ----wordcloud_crucial,echo=F--------------------------------------------------------------------------------------------------------------------------------
###### data table with food name + clustering result
result_mean_imp_spectrum_kmeans<-wide_data %>% select (FoodDescription) %>%
        mutate(group = data.crucial.mean.imp.spectrum.kmeans$cluster) %>% arrange(group)

###### Cluster 1 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==1]))

# inspect(docs)
docs <- docs %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc9.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()

###### Cluster 2 
docs<-Corpus(VectorSource(result_mean_imp_spectrum_kmeans$FoodDescription[result_mean_imp_spectrum_kmeans$group==2]))
# inspect(docs)
docs <- docs %>%
  tm_map(removeNumbers) %>%
  tm_map(removePunctuation) %>%
  tm_map(stripWhitespace)
dtm <- TermDocumentMatrix(docs) 
matrix <- as.matrix(dtm) 
words <- sort(rowSums(matrix),decreasing=TRUE) 
df <- data.frame(word = names(words),freq=words)
set.seed(1234) # for reproducibility 
png(here("figs", "wc10.png"), width=12, height=8, units="in", res=300)
wordcloud(words = df$word, freq = df$freq, min.freq = 1, max.words=100, random.order=FALSE, rot.per=0.35, colors=brewer.pal(8, "Dark2"))
dev.off()




## ----echo=F--------------------------------------------------------------------------------------------------------------------------------------------------
###### average calory count per cluster
avg_calories=data_mean_imp_crucial %>% mutate(group=data.crucial.mean.imp.spectrum.kmeans$cluster) %>% group_by(group) %>% summarize(mean.calories=mean(`ENERGY (KILOCALORIES)`),se.calories=sd(`ENERGY (KILOCALORIES)`) /sqrt(length(`ENERGY (KILOCALORIES)`)),std.calories=sd(`ENERGY (KILOCALORIES)`))

kcal_cluster2 = ggplot(avg_calories) +
    geom_bar(aes(x=group, y=mean.calories), stat="identity", fill="skyblue", alpha=0.7) +
    geom_errorbar( aes(x=group, ymin=mean.calories-se.calories, ymax=mean.calories+se.calories), width=0.4, colour="orange", alpha=0.9, size=1.3)


## ---- echo=F-------------------------------------------------------------------------------------------------------------------------------------------------
# PCA of nutrients 
data_mean_imp_crucial_temp1 =data_mean_imp_crucial %>% mutate(group = data.crucial.mean.imp.spectrum.kmeans$cluster) 
data_mean_imp_crucial.c1<-data_mean_imp_crucial_temp1[data_mean_imp_crucial_temp1$group==1,]
data_mean_imp_crucial.c2<-data_mean_imp_crucial_temp1[data_mean_imp_crucial_temp1$group==2,]


# Cluster 1 
sdzero=data.frame(apply(data_mean_imp_crucial.c1,2,sd))
whichone=which(sdzero == 0)
data_mean_imp_crucial.c1=data_mean_imp_crucial.c1%>%select(-c(names(data_mean_imp_crucial.c1)[whichone]))

data_mean_imp_crucial.pca <- prcomp(data_mean_imp_crucial.c1,scale=T,center=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = data_mean_imp_crucial.pca$rotation[,1],
                  gene = rownames(data_mean_imp_crucial.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = data_mean_imp_crucial.pca$rotation[,2],
                  gene = rownames(data_mean_imp_crucial.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients9 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 1") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))

# Cluster 2
sdzero=data.frame(apply(data_mean_imp_crucial.c2,2,sd))
whichone=which(sdzero == 0)
data_mean_imp_crucial.c2=data_mean_imp_crucial.c2%>%select(-c(names(data_mean_imp_crucial.c2)[whichone]))
data_mean_imp_crucial.pca <- prcomp(data_mean_imp_crucial.c2,scale=T,center=T)

top_k <- 20
## get pc1 and pc2
pc1 <- data.frame(loading = data_mean_imp_crucial.pca$rotation[,1],
                  gene = rownames(data_mean_imp_crucial.pca$rotation),
                  pc = "PC1")
pc2 <- data.frame(loading = data_mean_imp_crucial.pca$rotation[,2],
                  gene = rownames(data_mean_imp_crucial.pca$rotation),
                  pc = "PC2")
# get top_k of pc1 and pc2
pc1_top <- pc1 %>% arrange(-loading) %>% slice(1:top_k)
pc2_top <- pc2 %>% arrange(-loading) %>% slice(1:top_k)
important_nutrients10 = rbind(pc1_top, pc2_top) %>%
  ggplot(aes(x = reorder(gene, -loading), y = loading)) +
  geom_point() +
  ggtitle("Top Nutrients of Cluster 2") +
  xlab("Nutrient") +
  facet_wrap(~pc, nrow = 1, scales = "free_x") +
  theme_bw() + 
  theme(axis.text.x = element_text(angle = -45, hjust = 0, vjust = 1))


