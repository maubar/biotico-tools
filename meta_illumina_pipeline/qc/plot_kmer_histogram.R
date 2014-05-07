k17 <- read.table("k17.hist",header=F,sep=' ',col.names=c("kmer_freq","count"))
k31 <- read.table("k31.hist",header=F,sep=' ',col.names=c("kmer_freq","count"))
plot(k17$kmer_freq, k17$count,xlim=c(0,100),ylim=c(0,1000000))
plot(k31$kmer_freq, k31$count,xlim=c(0,100),ylim=c(0,1000000))
