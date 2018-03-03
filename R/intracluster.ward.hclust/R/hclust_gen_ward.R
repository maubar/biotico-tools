library(purrr)

#Helper functions
merge_two_clusters <- function(c1,c2,dist.m,iteration){
  new_cluster <- list(id=iteration,
                      members=c(c1$members,c2$members),
                      cum_dist=c1$cum_dist+c2$cum_dist)
  # Update cumulative distances
  for(c1.m in c1$members){
    for(c2.m in c2$members){
      new_cluster$cum_dist <- new_cluster$cum_dist + dist.m[c1.m,c2.m]
    }
  }
  return(new_cluster)
}

intracluster_var <- function(my.cluster){
  return(my.cluster$cum_dist / length(my.cluster$members))
}

intraclust_var_change <- function(c1,c2,c1_c2_merge){
  return(intracluster_var(c1_c2_merge) - intracluster_var(c1) - intracluster_var(c2))
}

generalWard <- function(dist.obj){
  # Check if it's a distance object, transform into distance matrix
  dist.m <- dist.obj
  if(class(dist.obj)=="dist"){
    dist.m <- as.matrix(dist.obj)
  }

  # Objects to cluster
  n.objects <- nrow(dist.m)

  #Define hclust object to return
  my.cluster <- list(
    merge=matrix(0,n.objects-1,2),
    height=rep(0.0,n.objects-1),
    order=rep(0.0,n.objects-1),
    labels=NULL,
    call=NULL,
    dist.method=ifelse(class(dist.obj)=="dist",attr(dist.obj,"method"),"custom")
  )
  class(my.cluster) <- "hclust"

  #Prepare initial cluster list to run the algorithm
  cluster_list <- map(1:n.objects,function(obj)list(id=-obj,members=obj,cum_dist=0))

  # Run the iteration of the hiearchical clustering. For N objects, N-1 iterations
  for(iteration in 1:(n.objects-1)){
    best_new_cluster <- NULL
    best_cluster_merge_indexes <- NULL
    best_intracluster_var_delta <- NULL

    #Brute force all possible mixes
    for(c1_idx in 1:(length(cluster_list)-1)){
      for(c2_idx in (c1_idx+1):length(cluster_list)){
        candidate_merge <- merge_two_clusters(cluster_list[[c1_idx]],cluster_list[[c2_idx]],dist.m,iteration)
        candidate_var_delta <- intraclust_var_change(cluster_list[[c1_idx]],cluster_list[[c2_idx]],candidate_merge)

        if(is.null(best_new_cluster) || (candidate_var_delta < best_intracluster_var_delta)){
          best_new_cluster <- candidate_merge
          best_cluster_merge_indexes <- c(c2_idx,c1_idx)
          best_intracluster_var_delta <- candidate_var_delta
        }
      }
    }

    # Update hclust object data structures
    my.cluster$merge[iteration,] <- c(cluster_list[[best_cluster_merge_indexes[[1]]]]$id,cluster_list[[best_cluster_merge_indexes[[2]]]]$id )
    my.cluster$height[iteration] <- intracluster_var(best_new_cluster)

    #Add new cluster to list, and delete merged clusters
    cluster_list[[length(cluster_list)+1]] <- best_new_cluster
    cluster_list <- cluster_list[-best_cluster_merge_indexes]
  }

  my.cluster$order <- best_new_cluster$members


  return(my.cluster)
}
