head(rowData(se))$id[order(start(head(rowData(se))))]

current <- as.character(unique(seqnames(rowData(se))))[1]
si <- which(seqnames(rowData(se))==current)

subr <- rowData(se)[si]

orn <- subr$id[order(start(subr))]

assays(se)[[1]]

rfh <- assays(se)[[1]][orn,]

ng <- dim(rfh)[1]

gs <- split(1:ng,round(as.numeric(cut(1:ng,20))))

smaller <- c()
for(i in 1:length(gs)){
  new <- apply(rfh[gs[[i]],],2,mean)
  smaller <- rbind(smaller,new)
}
rownames(smaller) <- 1:length(gs)


start(subr)[order(start(subr))]

rnames <- c()
for(i in 1:length(gs)){
  rnames <- c(rnames,paste(
  start(subr)[gs[[i]][1]],
  " - ",
  start(subr)[tail(gs[[i]],1)],
  sep=""))
}