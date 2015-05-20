#!/usr/bin/Rscript --vanilla

require(tools)

for (file.name in list.files(include.dirs=FALSE)) {
  print(file.name)
  name <- file_path_sans_ext(file.name)
  print(name)
  dat <- read.table(file.name,header=T,sep="\t", na.strings="na")
  dat$mean <- rowMeans(dat[,3:11])
  dat$mean2 <- apply(dat[,3:11],1,mean)
  dat$median <- apply(dat[,3:11],1,median)
  write.table(dat, file=paste(name, ".mean", sep=""), append=FALSE, sep="\t", eol="\n", col.names=TRUE, row.names=FALSE, quote=FALSE, qmethod="double")
}
