# Título: Redes booleanas estocásticas

# Nombre: Biología de sistemas

# Fecha: febrero 2021
#############################################################################

# objetivo : replicar las figuras de Méndez-López, L.F. et al., 2017. 


################## Preludio ##########################
rm(list=ls(all=TRUE)) #Limpiamos el Workspace (espacio de trabajo); ls() enlista 
#las variables, y rm() las borra "remove""

# fijar directorio de trabajo
#setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Bio Mates PCBIOL 2021 1/Prácticas computacionales/4 Redes booleanas estocásticas")

# cargar BoolNet
library(BoolNet) 

# cargar la red
net <- loadNetwork("Mendez_lopez2017_BMCSystemsBiology.txt") 

################## Figura 1 c #####################
# obtener atractores
attr_wt <- getAttractors(net)

# graficarlos
plotAttractors(attr_wt)

# obtener las cuencas de atracción, para cada uno de los 3 atractores

cuenca_epi <- attr_wt$attractors[[1]]$basinSize/2^9
cuenca_sen <- attr_wt$attractors[[2]]$basinSize/2^9
cuenca_mes <- attr_wt$attractors[[3]]$basinSize/2^9

# desplegar el porcentaje del tamaño de las cuencas
cuenca_epi
cuenca_sen
cuenca_mes

#verificar que sumen 1
sum(cuenca_epi, cuenca_sen, cuenca_mes)

## evaluar los efectos de NFkB
NFkB = fixGenes(net, "NFkB", 1)
NFkB_attr = getAttractors(NFkB)

getAttractorSequence(NFkB_attr,1)
getAttractorSequence(attr_wt,3)
NFkB_attr$attractors[[1]]$basinSize/2^8
attr_wt$attractors[[3]]$basinSize/2^9

getAttractorSequence(NFkB_attr,2)
getAttractorSequence(attr_wt,1)
NFkB_attr$attractors[[2]]$basinSize/2^8
attr_wt$attractors[[1]]$basinSize/2^9

getAttractorSequence(NFkB_attr,3)
getAttractorSequence(attr_wt,2)
NFkB_attr$attractors[[3]]$basinSize/2^8
attr_wt$attractors[[2]]$basinSize/2^9

# no cambian los atractores, solamente los tamaños de sus cuencas

######### Figura 2 (un mutante de su elección) #####

# .. re-declarar red, con mutantes
net_mut <-fixGenes(net,"ESE2", 1)
attr_mut <- getAttractors(net_mut)
plotAttractors(attr_mut)

# sacar cuencas de atracción
cuenca_mes_mut <- attr_mut$attractors[[1]]$basinSize/2^8
cuenca_mes_mut

# hacer lo propio con mutante 2
net_mut2 <- fixGenes(net,"Snai2", 0)
attr_mut2 <- getAttractors(net_mut2)
plotAttractors(attr_mut2)

# cuencas de atracción
cuenca_epi_sen <- attr_mut2$attractors[[1]]$basinSize/2^8
cuenca_epi <- attr_mut2$attractors[[2]]$basinSize/2^8

cuenca_epi_sen
cuenca_epi

###### Figura 3 a, izquierda, #####################

# a) Calcula la matriz de transición Pi, usando el código de José Dávila

# corre las funciones
source('Funciones_para_calcular_transiciones_entre_atractores.R')
# De:. Davila-velderrain J, Juarez-Ramiro L, Martinez-Garcia JC, Alvarez-Buylla ER. Methods for Characterizing the Epigenetic Attractors Landscape Associated with Boolean Gene Regulatory Networks. 2015;1-23. 
# https://github.com/JoseDDesoj/Epigenetic-Attractors-Landscape-R.

# Figure
# calcula la TMP
TPM <- Implicit.InterAttractor.Simulation(net, P.error = 0.05, Nreps = 100) # tarda un rato

TPM

## b) (para practicar): simular 1 cadena de markov, 100 pasos
Nsteps <- 100 # number of steps
pi0 = c(1,0,0) # Initial probability distributions: x(0)=E
v = vector("numeric", Nsteps) # creates an emplty vecor of size n (preallocate)
r = length(pi0) # size of the sample of the initial distribution 
v[1] = 1 #condición inicial

# Una sola realización de la cadena: iniciamos en un atractor, y luego....
for (i in 2:Nsteps) {
        
  v[i]=sample(r, 1, prob=TPM[v[i-1],]) # sample the new value: select the row in the probability matrix that gves the vector of probabiltiies according to current state
}

#graficar
matplot(v, type="l", lwd=2.5, col=3,  xlab="t", ylab="Attractor")

## ahora, vamos a iterar varias veces, 
## para sacar el promedio de las veces que 
# cada atractor fue visitado en cada paso de tiempo (), 
# es decir, la probabilidad (vista como frecuencia) 
# del atractor i al tiempo N, para i =1,2,3

# ojo, no vamos a guardar toda la cadena de markov; vamos a ir generando 
# vectores que, iterativamente, vayan actualizando el número de
# veces que cada cuenca de atracción fue visitada, en cada tiempo.

iterations = 100
# initialize number of times we have 1,2,3
numE = vector("numeric", Nsteps)# preallocate: create emplty vector of size N, number of time steps, for the probabilty of epithel
numS = vector("numeric", Nsteps)# preallocate: create emplty vector of size N, number of time steps, for the probabilty of epithel
numM = vector("numeric", Nsteps)# preallocate: create emplty vector of size N, number of time steps, for the probabilty of epithel

for (jj in 1:iterations) {
        
  for (i in 2:Nsteps) {
          
    v[i] = sample(r, 1, prob=TPM[v[i-1],]) # sample the new value: select the row in the probability matrix that gves the vector of probabiltiies according to current state
    
    if (v[i]==1) {
            
        numE[i]=numE[i] + 1
        
    } else if (v[i]==2) {
            
      numS[i]=numS[i]+1
      
    } else 
            
      numM[i]=numM[i]+1
  }
}

#para sacar frecuencias dividimos número de veces que atractor fue visitado 
# entre número de iteraciones
numE = numE/iterations
numS = numS/iterations
numM = numM/iterations


numE[1] = 1 #corregimos para la condición inicial (que calculamos fuera del
# loop, contador empezó en 2, no en 1)


## comparamos esto con resultados analíticos

Times <- 1:Nsteps # construye un vector con los pasos de tiempo

Pinit <- c(1,0,0) # distribución de probabilidad inicial: E=1

ProbEvolution <- matrix(0, length(Times)+1, ncol(TPM)) # creates a matrix where the evolution of the  vector will be stored
ProbEvolution[1,] <- Pinit # the first element is the initial probability

# itara el cálculo de la distribución de probabilidad de los atratores:

# x(t)=x(t-1)*TPM para todo t de 2 a numSteps

for(i in Times) ProbEvolution[i+1,] <- ProbEvolution[i,] %*% (TPM) 

# grafica y compara resultados analíticos con numéricos
par(new = FALSE)
matplot(numE, type="p", lwd=3, col=5, ylab=NULL, xlab=NULL)
par(new=TRUE)
matplot(ProbEvolution[,1], type="l", lwd=3, col=6, ylab="P(E)", xlab="time steps")

par(new=FALSE)
matplot(numS, type="p", lwd=3, col=5, ylab=NULL, xlab=NULL)
par(new=TRUE)
matplot(ProbEvolution[,2], type="l", lwd=3, col=6, ylab="P(S)", xlab="time steps")

par(new=FALSE)
matplot(numM, type="p", lwd=3, col=5, ylab=NULL, xlab=NULL)
par(new=TRUE)
matplot(ProbEvolution[,3], type="l", lwd=3, col=6, ylab="P(M)", xlab="time steps")


### Calcular el tiempo promedio mínimo 
# ("Mean first passage time" )
# para pasar de un atractor al otro; todas las 3 por 3 combinaciones
#entre los tres atractores

MFPT=matrix(0, nrow = 3, ncol = 3) #crea matrix en la que vaciaremos este tiempo
numTimeSteps=5000 # sufficiently large to allow the passage from one attractor to th othr

# for all the different initial conditions, stratring from attractor 1, 2 or 3
for (initial_condition in 1:3){
  
  v[1]=initial_condition
  #iterate many times, to take the average
  for (jj in 1:iterations){
    # simulate a markov chain
    for (i in 2:numTimeSteps){
      v[i]=sample(r, 1, prob=TPM[v[i-1],]) # sample the new value: select the row in the probability matrix that gves the vector of probabiltiies according to current state
    }
    # in this simulated markov chain, find the first time that the attractor j is attained
    for (attrnum in 1:3){
      timetoattractor=min(which(v==attrnum))
      if (is.finite(timetoattractor)){
        # we will sum this time to then take the average
        MFPT[initial_condition, attrnum]=MFPT[initial_condition, attrnum]+timetoattractor
        
      }
    }
  }
}

MFPT = (MFPT-iterations)/iterations
MFPT

library(fields)

x = 1:3
y = 1:3

image.plot(y, x,MFPT, ylab = "Atractor i ", xlab = "Atractor j",legend.lab="MFPT")

