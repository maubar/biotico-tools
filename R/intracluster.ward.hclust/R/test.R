### TEST!
rm(list=ls())

library(dendextend)
source("R/hclust_gen_ward.R")

#Define data
my.data <- rbind(  c(0,0) , c(0,1), c(1,0) , c(1,1), c(0.5,0.5),
                   c(10,10), c(10,11), c(11,10), c(11,11)
)

rownames(my.data) <- c("a1","a2","a3","a4","a5", "b1","b2","b3","b4")


plot(my.data)

a <- hclust(dist(my.data),method="ward.D")

my.hclust <- generalWard(dist(my.data))

my.hclust %>% as.dendrogram() %>%
  set("leaves_pch", 19) %>%
  set("labels_cex", 0.5) %>%
  set("branches_k_color", k = 2) %>%
  plot()
