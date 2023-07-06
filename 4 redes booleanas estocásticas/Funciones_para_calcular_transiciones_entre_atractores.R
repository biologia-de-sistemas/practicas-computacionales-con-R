# instalar paquetes requeridos
list_of_packages <- c("BoolNet","expm","igraph","markovchain","igraph",
                      "RColorBrewer")

new_packages <- list_of_packages[!(list_of_packages %in%  installed.packages()[,"Package"])]
if (length(new_packages)) install.packages(new_packages)

#Libraries
library(BoolNet)
library(expm)
library(igraph)
library(markovchain)
library(igraph)
library(RColorBrewer)
#############################################################################################
# FUNCTIONS
#############################################################################################
AttrTable <- function(Net) {
  # Graphs and returns attractor table
  #Extracts the attractors in a table form from the attractor object created by BoolNet function "getAttractors"
  attrs <- getAttractors(Net)
  Attrs <- print(plotAttractors(attrs))
  Attrs <- Attrs[[1]]
  return(Attrs)
}
####################################################
 Get.Attractors.Landscape <- function(Network) {
   attrs <- getAttractors(Network)
   TransTable <- getTransitionTable(attrs)

   Inits <- TransTable[[1]]
   for(i in 2:length(Network$genes)) Inits <- cbind(Inits, TransTable[[i]])

   AttractsClusterVector <- TransTable$attractorAssignment

   return(list(Inits, AttractsClusterVector))
 }
####################################################
Estimate.Basins.Distribution <- function(Network) {
  AttrsLandscape <- Get.Attractors.Landscape(Network)
  Basins <- AttrsLandscape[[2]]
  TablaAttr <- AttrTable(Network)
  Dist.StateSpace <- numeric(ncol(TablaAttr))
  names(Dist.StateSpace) <- names(table(Basins))
  Dist.StateSpace[match(names(table(Basins)), names(Dist.StateSpace))] <- as.numeric(table(Basins)/length(Basins))
  return(Dist.StateSpace)
}
####################################################
 StateSpace.OneStep.Mapping <- function(Network) {
   MappingStructure <- list()
   MappingStructure[[1]] <- Get.Attractors.Landscape(Network)[[1]]
   MappingStructure[[2]] <- t(apply(MappingStructure[[1]], 1,  function(i) stateTransition(Network, i)))
   colnames(MappingStructure[[1]]) <- colnames(MappingStructure[[2]])
   return(MappingStructure)
 }
####################################################
Collapse <- function(MappingStructure) {
  MappingStructure.Collapsed <- list()
  MappingStructure.Collapsed[[1]] <- apply(MappingStructure[[1]], 1,  function(i) paste(i, collapse=""))
  MappingStructure.Collapsed[[2]] <- apply(MappingStructure[[2]], 1,  function(i) paste(i, collapse=""))
  return(MappingStructure.Collapsed)
}
####################################################
Implicit.InterAttractor.Simulation <- function(Network, P.error, Nreps) {
  # Simulates binomial mutations in all space
  #Get State Space and x(t+1) = f(x) for each state
  AttrsLandscape <- Get.Attractors.Landscape(Network)
  Ngenes <- ncol(AttrsLandscape[[1]])
  Nattractors <- length(names(table(AttrsLandscape[[2]])))
  StateSpace <- AttrsLandscape[[1]]
  Next.T.StateSpace <- StateSpace*0
  Next.T.StateSpace <- t(sapply(1:nrow(Next.T.StateSpace), function(i) Next.T.StateSpace[i,] <- stateTransition(Network, StateSpace[i,])))
  colnames(StateSpace) <- colnames(Next.T.StateSpace)
  Char.State.Space <- apply(StateSpace, 1, function(i) paste(i, collapse=""))

  #Create Attractors Transition Probability Matrix
  T.Prob.Mat <- matrix(0, Nattractors, Nattractors)
  rownames(T.Prob.Mat) <- 1:Nattractors
  colnames(T.Prob.Mat) <- 1:Nattractors
  AttrsInd <- as.numeric(colnames(T.Prob.Mat))

  #Create Muation Indicator Vector
  MutMatrix <- rbinom(Nreps*Ngenes*nrow(StateSpace), 1, P.error)

  #Create concatenated vector of X(t+1)
  NextSs <- as.numeric(apply(Next.T.StateSpace,1, function(i) rep(i, Nreps)))

  #Simulate "errors" in X(t+1)
  Mutind <- which(MutMatrix==1)
  ZeroInd <- Mutind[which(NextSs[Mutind]==1)]
  UnoInd <- Mutind[which(NextSs[Mutind]==0)]
  NextSs[ZeroInd] <- 0
  NextSs[UnoInd] <- 1

  #Split and count states in X(t+1). Match them with basins.
  NextSs <- apply(matrix(NextSs, Nreps*nrow(StateSpace), Ngenes, byrow=TRUE), 1, function(i) paste(i, collapse=""))
  #NextSs.L <- split(NextSs, rep(1:nrow(StateSpace), each=Nreps))
  NextSs.LL <- lapply(split(NextSs, rep(1:nrow(StateSpace), each=Nreps)), table)
  Basins.L <-  lapply(NextSs.LL, function(i) AttrsLandscape[[2]][match(names(i), Char.State.Space)])

  #Update Attractors Transition Probability Matrix
  for(j in 1:nrow(StateSpace)) T.Prob.Mat[AttrsLandscape[[2]][j], ] <- T.Prob.Mat[AttrsLandscape[[2]][j], ] + sapply(1:length(AttrsInd), function(i) sum(NextSs.LL[[j]][which(AttrsInd[i]==Basins.L[[j]])]))

  #Normalize Attractors Transition Probability Matrix
  #SUMS <- apply(T.Prob.Mat, 1, sum)
  #for(i in 1:nrow(T.Prob.Mat)) T.Prob.Mat[i,] <- T.Prob.Mat[i,]/SUMS[i]
  T.Prob.Mat <- t(sapply(1:nrow(T.Prob.Mat), function(i) T.Prob.Mat[i,]/sum(T.Prob.Mat[i,])))
  return(T.Prob.Mat)
}

####################################################
Stochastic.InterAttractor.Simulation <- function(Network, error, Nsteps) {

  Attractors.Landscape <- Get.Attractors.Landscape(Network)
  AttractorsLandscape.Splited <- split(as.data.frame(Attractors.Landscape[[1]]), Attractors.Landscape[[2]])
  Mapping.Structure <- StateSpace.OneStep.Mapping(Network)
  Attractors.Mapping <- Get.Attractors.Mapping(Mapping.Structure, Attractors.Landscape)

  Simulation <- numeric(Nsteps)

  i <- sample(1:nrow(Mapping.Structure[[1]]),1)
  Simulation[1] <- Attractors.Mapping[i,3]
  Si <- Mapping.Structure[[2]][i,]

  for(Sims in 2:length(Simulation)) {
    Si.i <- Si
    Mut <- sample(c(0,1), length(Si), prob=c(1-error,error), rep=TRUE)

    if(sum(Mut==1)>0) {
      Cuales <- which(Mut==1)
      for(j in 1:length(Cuales)) {
	if(Si[Cuales[j]]==0) Si.i[Cuales[j]] <- 1
	if(Si[Cuales[j]]==1) Si.i[Cuales[j]] <- 0
      }
    }

    Simulation[Sims] <- Attractors.Mapping[[3]][which(paste(Si.i, collapse="")==Attractors.Mapping[[1]])]

    i <- sample(1:nrow(AttractorsLandscape.Splited[[Simulation[Sims]]]),1)
    Si <- AttractorsLandscape.Splited[[Simulation[Sims]]][i,]
    print(Sims)
    print("################################################")
  }

  FittedMLE <- markovchainFit(data = as.character(Simulation), method = "mle", name = "IAT MLE")
  FittedMLE.TPM <- as.matrix(FittedMLE$estimate[1:length(unique(Simulation))])

  return(list(Simulation, FittedMLE.TPM))
}

####################################################
BooleanNet.Stochastic.Matrix.Model <- function(Network, Prob) {
  #Network <- loadNetwork(Net)
  MappingStructure <- StateSpace.OneStep.Mapping(Network)
  MappingStructure.Collapsed <- Collapse(MappingStructure)
  MappingStructure.Map <- cbind(MappingStructure.Collapsed[[1]], MappingStructure.Collapsed[[2]])
  T.Graph <- graph.edgelist(MappingStructure.Map)
  T.mat <- as.matrix(get.adjacency(T.Graph))

  P <- T.mat*0

  for(S in 1:nrow(T.mat)) {
    Ti <- as.numeric(unlist(strsplit(rownames(T.mat)[S],"")))

    Tis <- t(sapply(colnames(T.mat), function(i) as.numeric(unlist(strsplit(i,"")))))
    Hammings <- apply(Tis, 1, function(i) sum(!Ti==i))

    P[S,] <- dbinom(x=Hammings, size=length(Ti), prob=Prob)
  }

  diag(P) <- 0

  for(i in 1:nrow(P)) P[i,] <- P[i,]/sum(P[i,])

  T.tild <- (1-Prob)^length(Ti) * T.mat + P
  for(i in 1:nrow(T.tild)) T.tild[i,] <- T.tild[i,]/sum(T.tild[i,])

  return(list(Matrix.T = T.mat, Matrix.P = P, Matrix.Ttilde = T.tild))
}

####################################################
BooleanNet.Stochastic.Matrix.Model.N <- function(Net, Prob) {

  MappingStructure <- StateSpace.OneStep.Mapping(Net)
  MappingStructure.Collapsed <- Collapse(MappingStructure)
  MappingStructure.Map <- cbind(MappingStructure.Collapsed[[1]], MappingStructure.Collapsed[[2]])
  T.Graph <- graph.edgelist(MappingStructure.Map)
  T.mat <- as.matrix(get.adjacency(T.Graph))

  P <- T.mat*0

  #   for(S in 1:nrow(T.mat)) {
  #     Ti <- as.numeric(unlist(strsplit(rownames(T.mat)[S],"")))
  #
  #     Tis <- t(sapply(colnames(T.mat), function(i) as.numeric(unlist(strsplit(i,"")))))
  #     Hammings <- apply(Tis, 1, function(i) sum(!Ti==i))
  #
  #     P[S,] <- dbinom(x=Hammings, size=length(Ti), prob=Prob)
  #   }
  #
  StatesRows <- t(sapply(1:nrow(T.mat), function(i) as.numeric(unlist(strsplit(rownames(T.mat)[i],"")))))
  StatesCols <- t(sapply(colnames(T.mat), function(i) as.numeric(unlist(strsplit(i,"")))))
  Hammings <- t(sapply(1:nrow(StatesRows), function(j) apply(StatesCols, 1, function(i) sum(!StatesRows[j]==i))))
  P <- t(sapply(1:nrow(Hammings), function(i) dbinom(x=Hammings[i,], size=ncol(StatesRows), prob=Prob)))
  rownames(P) <- rownames(T.mat)
  colnames(P) <- colnames(T.mat)

  diag(P) <- 0

  #for(i in 1:nrow(P)) P[i,] <- P[i,]/sum(P[i,])
  P <- t(sapply(1:nrow(P), function(i) P[i,]/sum(P[i,])))

  T.tild <- (1-Prob)^ncol(StatesRows) * T.mat + P

  #for(i in 1:nrow(T.tild)) T.tild[i,] <- T.tild[i,]/sum(T.tild[i,])
  T.tild <- t(sapply(1:nrow(T.tild), function(i) T.tild[i,]/sum(T.tild[i,])))
  rownames(P) <- rownames(T.mat)
  colnames(P) <- colnames(T.mat)
  rownames(T.tild) <- rownames(T.mat)
  colnames(T.tild) <- colnames(T.mat)

  return(list(Matrix.T = T.mat, Matrix.P = P, Matrix.Ttilde = T.tild))
}

####################################################
Calculate.MFPT.Matrix <- function(P, AttrsNames=1:ncol(P)) {

  for(i in 1:nrow(P)){
    x<-1:nrow(P)
    z<-matrix(0,nrow(P)-1,1)
    unos<-matrix(1,nrow(P)-1,1)
    id<-diag(ncol(P)-1)
    x1<-cbind(z,id)
    x2<-cbind(unos,P[x[x!=i],x[x!=i]])
    G<-rbind(x1,x2)   ##matrix G start process

    for(j in 1:(nrow(P)-1)){
      T6<-G[1:nrow(G)-1,1:ncol(G)-1]
      U6<-as.matrix(G[1:nrow(G)-1,ncol(G)])
      R6<-t(as.matrix(G[nrow(G),1:ncol(G)-1]))
      Q6<-G[nrow(G),ncol(G)]
      G5<-T6+U6%*%(1-Q6)^(-1)%*%R6
      G<-G5
    }

    H<-vector()
    H[i]<-0

    for(s in 1:nrow(P) ){
      if((s==1)&(s!=i)) H[s]<-G[s]
      if((s<i)&(s>1)) H[s]<-G[s]
      if(s>i) H[s]<-G[s-1]
    }

    H<-as.matrix(H)

    if(i>1) M<-cbind(M,H)
    else M<-H
  }
  rownames(M) <- AttrsNames
  colnames(M) <- AttrsNames
  return(M)
}
####################################################
Calc.Stationary.Distribution <- function(TPM, AttrsNames) {
  e=eigen(t(TPM))$vectors[,1]
  Stationary <- e/sum(e)
  names(Stationary) <- AttrsNames
  return(Stationary)
}
####################################################
MFPT.Transition.Rates <- function(MFPTsMAT) {
  Mat <- (1/MFPTsMAT) - 1/t(MFPTsMAT)
  diag(Mat) <- 0
  return(Mat)
}
####################################################
Plot.Attractor.Global.Ordering <- function(TransitionRates) {
  Graph <- graph.adjacency(TransitionRates, mode="directed", add.rownames=TRUE, weighted=TRUE)
  igraph.options(vertex.size=25)
  E(Graph)[E(Graph)$weight>0]$color <- "red"
  plot(Graph)
}
####################################################
Probability.Distribution.Temporal.Evolution <- function(P, Time=100) {
  Times <- 1:Time
  Pinit <- P[1,]*0
  Pinit <- rep(1/length(Pinit), length(Pinit))

  ProbEvolution <- matrix(0, length(Times)+1, ncol(P))
  ProbEvolution[1,] <- Pinit

  for(i in Times) ProbEvolution[i+1,] <- Pinit %*% (P%^%i)
  names(ProbEvolution) <- colnames(P)
  return(ProbEvolution)
}
####################################################
Sub.Attractor.Transition.Matrix <- function(CompleteTPM, net) {
  AL <- Get.Attractors.Landscape(net)
  AL[[1]] <- apply(AL[[1]], 1,  function(i) paste(i, collapse=""))
  AL[[1]][AL[[2]]==1]

  Atts.ind <- as.numeric(names(table(AL[[2]])))

  Sub.Atts.TPM <- matrix(0, length(Atts.ind),length(Atts.ind))
  row.names(Sub.Atts.TPM) <- Atts.ind
  colnames(Sub.Atts.TPM) <- Atts.ind

  for(i in 1:length(Atts.ind)) {
    for(j in 1:length(Atts.ind)) {
      Subi <- CompleteTPM[which(rownames(CompleteTPM)%in%AL[[1]][AL[[2]]==i]),which(rownames(CompleteTPM)%in%AL[[1]][AL[[2]]==j])]
      Sub.Atts.TPM[i,j] <- sum(apply(Subi, 1, sum))
    }
  }

  SUMS <- apply(Sub.Atts.TPM, 1, sum)
  for(i in 1:nrow(Sub.Atts.TPM)) Sub.Atts.TPM[i,] <- Sub.Atts.TPM[i,]/SUMS[i]
  return(Sub.Atts.TPM)
}

####################################################
Plot.Probability.Time.Order <- function(TPM, Initial, AttrsNames, timeF) {
  Times <- 1:timeF
  P <- as.matrix(TPM)
  Pinit <- P[1,]*0
  Pinit[which(AttrsNames==Initial)] <- 1

  #barplot(-log(Pinit %*% (ProbMat%^%10000)))
  ProbEvolution <- matrix(0, length(Times)+1, ncol(P))
  ProbEvolution[1,] <- Pinit

  for(i in Times) ProbEvolution[i+1,] <- ProbEvolution[i,] %*% (P)

  ORD <- 1:length(AttrsNames)
  names(ORD) <- AttrsNames[order(apply(ProbEvolution, 2, which.max))]
  print(apply(ProbEvolution, 2, function(i) which.max(round(i,3))))
  print(ORD)
  barplot(ORD, col=match(names(ORD), AttrsNames), main=paste(Initial,"->", "Temporal sequence of cell-fate attainment"), axes=FALSE)
}
####################################################
PlotTemporalProbabilityPDF <- function(TPM, Initial, AttrsNames, timeF, File) {
  #pdf(File)
  #par(mfrow=c(2,1))
  Plot.Probability.Evolution(TPM, Initial, AttrsNames, timeF)
  #Plot.Probability.Time.Order(TPM, Initial, AttrsNames, timeF)
  #par(mfrow=c(1,1))
  #dev.off()
}
####################################################
Plot.Probability.Evolution <- function(TPM, Initial, AttrsNames, timeF) {
  Times <- 1:timeF
  P <- as.matrix(TPM)
  Pinit <- P[1,]*0
  Pinit[which(AttrsNames==Initial)] <- 1
  
  #barplot(-log(Pinit %*% (ProbMat%^%10000)))
  ProbEvolution <- matrix(0, length(Times)+1, ncol(P))
  ProbEvolution[1,] <- Pinit
  
  for(i in Times) ProbEvolution[i+1,] <- ProbEvolution[i,] %*% (P)
  #######################
  #NormProbEvolution <- ProbEvolution*0
  #print(ProbEvolution)
  #for(i in 2:nrow(ProbEvolution)) NormProbEvolution[i,] <- ProbEvolution[i,]/max(ProbEvolution[i,])
  #print(ProbEvolution)
  #ProbEvolution <- NormProbEvolution
  #print(ProbEvolution)
  #print(ProbEvolution)
  #print(apply(ProbEvolution, 2, function(i) round(i, 3)))
  #ProbEvolution <- apply(ProbEvolution, 2, function(i) round(i, 3))
  #NormProbEvolution <- ProbEvolution*0
  #print(ProbEvolution)
  #for(i in 2:nrow(ProbEvolution)) NormProbEvolution[i,] <- ProbEvolution[i,]/max(ProbEvolution[i,])
  #print(ProbEvolution)
  #ProbEvolution <- NormProbEvolution
  #######################
  
  matplot(ProbEvolution, type="l", lwd=2.5, col=1:length(AttrsNames), main=paste(Initial,"->", "Probability Temporal Evolution"), xlab="t", ylab="P(A)")
  #matplot(NormProbEvolution, type="l", lwd=2)
  #AttrsNames[apply(ProbEvolution, 1, which.max)]
  print(AttrsNames[order(apply(ProbEvolution, 2, which.max))])
  abline(v=apply(ProbEvolution, 2, function(i) which(max(round(i, 3))==round(i, 3))[1]), col=1:length(AttrsNames), lwd=2.5)
  #legend(40, 0.8, AttrsNames, col=1:length(AttrsNames), lwd=2)
  colnames(ProbEvolution) <- AttrsNames
  return(ProbEvolution)
}
#############################################################################################
#############################################################################################