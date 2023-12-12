## first biogeographic (DEC) model (still in progress)

fitBiogeog<-function(tree,x,model="DEC",...){
  if(hasArg(min.q)) min.q<-list(...)$min.q
  else min.q<-1e-12
  if(hasArg(max.q)) max.q<-list(...)$max.q
  else max.q<-max(nodeHeights(tree))*100
  if(model=="DEC") model<-"ARD"
  if(hasArg(quiet)) quiet<-list(...)$quiet
  else quiet<-FALSE
  if(is.factor(x)) x<-setNames(as.character(x),names(x))
  if(is.matrix(x)) X<-strsplit(colnames(x),"+",fixed=TRUE)
  else X<-strsplit(x,"+",fixed=TRUE)
  ns<-sapply(X,length)
  ## get the regions
  regions<-sort(unique(unlist(X)))
  nn<-length(regions)
  ## fix the order of the input data
  if(is.matrix(x)){
    Levs<-sapply(X,function(x) paste(sort(x),collapse="+"))
    colnames(x)<-Levs
  } else x<-sapply(X,function(x) paste(sort(x),collapse="+"))
  ss<-vector()
  for(i in 1:length(regions))
    ss<-c(ss,apply(Combinations(length(regions),i,regions),
      1,paste,collapse="+"))
  if(!any(ss=="0")) ss<-c(0,ss)
  tmodel<-matrix(0,length(ss),length(ss),dimnames=list(ss,ss))
  poly<-strsplit(ss,"+",fixed=TRUE)
  index<-0
  tmodel[1,]<-rep(0,ncol(tmodel))
  E<-setNames(rep(0,nn),regions)
  D<-matrix(rep(0,nn*nn),nn,nn,dimnames=list(regions,regions))
  if(model=="ER"){
    E[]<-1
    D[]<-2
    diag(D)<-0
  } else if(model=="ARD"){
    E[]<-1:nn
    ind<-(nn+1)
    for(i in 1:nn) for(j in 1:nn){
      if(i!=j){
        D[i,j]<-ind
        ind<-ind+1
      }
    }
  }
  for(i in 2:nrow(tmodel)){
    for(j in 1:ncol(tmodel)){
      if(i==j) tmodel[i,j]<-0
      else {
        INT<-intersect(poly[[i]],poly[[j]])
        SDij<-setdiff(poly[[i]],poly[[j]])
        SDji<-setdiff(poly[[j]],poly[[i]])
        if(length(INT)==0){
          if(length(SDji)==1&&length(SDij)==1&&SDji=="0") tmodel[i,j]<-E[SDij]
        } else {
          if(length(SDji)==0&&length(SDij)==1) tmodel[i,j]<-E[SDij]
          else if(length(SDij)==0&&length(SDji)==1){ 
            print(paste(paste(poly[[i]],collapse="+"),"->",paste(poly[[j]],collapse="+")))
            tmp<-D[poly[[i]],SDji]
            tmodel[i,j]<-if(length(tmp)==1) tmp else paste(tmp,collapse="+")
          }
        }
      }
    }
  }
  if(!quiet){
    cat("\nThis is the design matrix of the fitted model.\nDoes it make sense?\n\n")
    print(tmodel)
    cat("\n")
    flush.console()
  }
  if(is.matrix(x)){
    X<-matrix(0,nrow(x),length(ss),dimnames=list(rownames(x),ss))
    X[rownames(x),colnames(x)]<-x
  } else X<-to.matrix(x,ss)
  ## initialize q
  q.init<-rexp(n=max(D))
  ## optimize the likelihood
  pw<-reorder(tree,"postorder")
  lik_func<-function(p) -biogeog_pruning(exp(p),tree=pw,x=X,model=tmodel)
  dec_fit<-nlminb(log(q.init),lik_func)
  print(dec_fit)
  return(dec_fit)
}

biogeog_pruning<-function(q,tree,x,model=NULL,...){
  if(hasArg(return)) return<-list(...)$return
  else return<-"likelihood"
  pw<-if(!is.null(attr(tree,"order"))&&
      attr(tree,"order")=="postorder") tree else 
        reorder(tree,"postorder")
  k<-ncol(x)
  if(hasArg(pi)) pi<-list(...)$pi
  else pi<-rep(1/k,k)
  Q<-matrix(0,k,k)
  colnames(Q)<-rownames(Q)<-colnames(model)
  for(i in 1:nrow(model)){
    for(j in 1:ncol(model)){
      if(model[i,j]!="0"){
        ind<-as.numeric(strsplit(model[i,j],"+",fixed=TRUE)[[1]])
        Q[i,j]<-sum(q[ind])
      }
    }
  }
  diag(Q)<--rowSums(Q)
  print(Q)
  L<-rbind(x[pw$tip.label,],
    matrix(0,tree$Nnode,k,
      dimnames=list(1:tree$Nnode+Ntip(tree))))
  nn<-unique(pw$edge[,1])
  pp<-vector(mode="numeric",length=length(nn))
  root<-min(nn)
  for(i in 1:length(nn)){
    ee<-which(pw$edge[,1]==nn[i])
    PP<-matrix(NA,length(ee),k)
    for(j in 1:length(ee)){
      P<-expm(Q*pw$edge.length[ee[j]])
      PP[j,]<-P%*%L[pw$edge[ee[j],2],]
    }
    L[nn[i],]<-apply(PP,2,prod)
    if(nn[i]==root){
      if(pi[1]=="fitzjohn") pi<-L[nn[i],]/sum(L[nn[i],])
      L[nn[i],]<-pi*L[nn[i],]
    }
    pp[i]<-sum(L[nn[i],])
    L[nn[i],]<-L[nn[i],]/pp[i]
  }
  prob<-sum(log(pp))
  print(prob)
  if(return=="likelihood") 
    if(is.na(prob)||is.nan(prob)) 
      return(-Inf) else return(prob)
  else if(return=="conditional") L
  else if(return=="pi") pi
}
