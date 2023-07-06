# Título: Optimización paramétrica

# Nombre: Biología de sistemas

# Fecha: febrero 2021
######################################################################################3#
## parameter estimation in R
# http://tbb.bio.uu.nl/rdb/grindR/

# borrar el contenido en memoria
rm(list=ls())

# establecer directorio de trabajo
# setwd("C:/Users/Elisa/Dropbox/Docencia adscripción biomédicas/Posgrado sostenibilidad 2020/Modulo EDOs Elisa Domínguez Hüttinger/Prácticas_Computacionales_ODEs_en_R")

# (1) cargar Grind.R
source("Grind.R")

# (2) establecer modelo
model <- function(time, state, parms) {
        with(as.list(c(state,parms)),{
                
                dx <- xpre*kprodEprod - x*kdegEdeg - x*y*kdimEdim + (ytot-y)*kdisEdis
                dy <- - x*y*kdimEdim + (ytot-y)*kdisEdis
                dxy <- - dy
                return(list(c(dx, dy, dxy)))
        })
}

# (3) Declarar primer estimación para los parámetros
p <- c(xpre = 10, kprodEprod = 0.5, kdegEdeg = 0.5, kdimEdim = 10, kdisEdis = 1, ytot = 6)

# Establecer condiciones iniciales (las cuales, efectivamente, son parámetros)
s <- c(x = 0, y = 5, xy = 1)

# (4) Introducir datos experimentales
t_exp <- c(0, 1, 3, 6)
xy_exp <- c(0.5, 2.5, 3, 1.8)/0.5 
# Guardar en un dataframe (necesario para usar la función fit()). 
# Los nombres de las variables en el dataframe deben ser igual a sus nombres en el modelo
data <- data.frame(time = t_exp, xy = xy_exp)

# (5) Correr la optimización
w <- c("kdimEdim", "ytot") # w contiene los nombres de los parámetros libres
f <- fit(legend = FALSE, free = w, tstep = 0.01) 


# Revisar error estándar, estadístico t, valor p, etc.
summary(f)

# Es posible revisar los resultados del ejemplo con f$par
f$par
f$ssr

# (6) some options
# assume you have more datasets
#data_1 <- data; # the original data
# invent a new dataset (only for practice purposes of course)
#data_2 <- data; # the original data
#data_2$XY_t <- c(1.4, 4, 5) # empty to fill it in again
#f <- fit(list(data_1,data_2),  free=w, tstep=0.01, differ = "Ytot") # which gives the same result
# differ are the parameters that are allowed to be different

p_opt=p
p_opt[4]=f$par[1]
p_opt[6]=f$par[2]

p_opt
par(pty="s") # axis square
run(parms=p_opt, legend=FALSE, tstep=0.001, tmax=3, main=c(paste("cost=", toString(round(f$ssr, digits = 3) ), sep=" ")))
points(t_exp, xy_exp, col="darkgreen", lty=5)
