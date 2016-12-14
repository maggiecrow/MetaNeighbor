neighbor.voting.LeaveOneExpOut <- function(exp.labels,cell.labels,network,means=T){
  
  # genes.label : needs to be in 1s and 0s
  l <- dim(cell.labels)[2]
  c <- dim(cell.labels)[1]
  e <- unique(exp.labels)
  
  #print("Make genes label CV matrix")
  test.cell.labels = matrix(cell.labels, nrow=c, ncol=length(e)*l)
  
  exp.cols=rep(e,each=l)
  
  for (i in 1:length(e)){
    d<-which(exp.labels==i)
    a<-which(exp.cols==i)
    test.cell.labels[d,a]<-0
  }
  
  #print("Get sums - mat. mul.")
  #sumin    = ( t(network) %*% test.genes.labels)
  sumin    = ( (network) %*% test.cell.labels)
  
  #print("Get sums - calc sumall")
  sumall   = matrix(apply(network,2,sum), ncol = dim(sumin)[2], nrow=dim(sumin)[1])
  
  #print("Get sums - calc predicts")
  predicts = sumin/sumall
  
  #print("Hide training data")
  nans = which(test.cell.labels == 1, arr.ind=T)
  predicts[nans] <- NA
  
  #Hide other experiment data 
  for (i in 1:length(e)){
    d<-which(exp.labels!=i)
    a<-which(exp.cols==i)
    predicts[d,a]<-NA
  }
  
  #print("Rank test data")
  predicts = apply(abs(predicts), 2, rank,na.last="keep",ties.method="average")
  
  filter = matrix(cell.labels, nrow=c, ncol=length(e)*l)
  for (i in 1:length(e)){
    d<-which(exp.labels!=i)
    a<-which(exp.cols==i)
    filter[d,a]<-NA
  }
  negatives = which(filter == 0, arr.ind=T)
  positives = which(filter == 1, arr.ind=T)
  
  predicts[negatives] <- 0
  
  #print("Calculate ROC - np")
  np = colSums(filter,na.rm=T) # Postives
  
  #print("Calculate ROC - nn")
  nn = apply(filter,2,function(x) sum(x==0,na.rm=T))     # Negatives
  
  #print("Calculate ROC - p")
  p =  apply(predicts,2,sum,na.rm=T)
  
  #print("Calculate ROC - rocN")
  rocNV = (p/np - (np+1)/2)/nn
  rocNV = matrix(rocNV, ncol=length(e), nrow=l)
  colnames(rocNV)=e
  rownames(rocNV)=colnames(cell.labels)
  
  if(means==T){
    scores=list(rowMeans(rocNV,na.rm=T))
  }
  else {
    scores = list(rocNV)      
  }
}

run_MetaNeighbor <- function(data, experiment_labels, celltype_labels, genesets, file_ext) {

    ROCs<-vector("list",length=length(genesets))
    names(ROCs)=names(genesets)
    nv.mat<-matrix(0,ncol=dim(celltype_labels)[2],nrow=length(genesets))
    rownames(nv.mat)=names(genesets)
    colnames(nv.mat)=colnames(celltype_labels)
    
    for (l in 1:length(genesets)){
      print(l)
      geneset=genesets[[l]]
      m<-match(rownames(data),geneset)
      dat.sub=data[!is.na(m),]
      dat.sub=cor(dat.sub,method="s")
      dat.sub=as.matrix(dat.sub)
      rank.dat=dat.sub
      rank.dat[]=rank(dat.sub,ties.method="average",na.last="keep")
      rank.dat[is.na(rank.dat)]=0
      rank.dat=rank.dat/max(rank.dat)
      ROCs[[l]]=neighbor.voting.LeaveOneExpOut(experiment_labels,celltype_labels,rank.dat,means=F)
    }
    for(i in 1:length(ROCs)){
      nv.mat[i,]=rowMeans(ROCs[[i]][[1]],na.rm=T)
    }
    save(ROCs,file=paste(file_ext,"IDscore.list.Rdata",sep="."))
    save(nv.mat, file=paste(file_ext,"IDscore.matrix.Rdata",sep="."))
    return(rowMeans(nv.mat,na.rm=T))
  }