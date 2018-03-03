### TEST!
rm(list=ls())

library(dendextend)
source("R/hclust_gen_ward.R")

#Define data
my.data <- rbind(  c(0,0) , c(0,1), c(0.5,0.5) ,
                   c(20,20), c(20,21), c(20.5,20.5)
)

plot(my.data)

a <- hclust(dist(my.data),method="ward.D")

my.hclust <- generalWard(dist(my.data))

my.hclust %>% as.dendrogram() %>%
  set("leaves_pch", 19) %>%
  set("labels_cex", 0.5) %>%
  set("branches_k_color", k = 2) %>%
  plot()
