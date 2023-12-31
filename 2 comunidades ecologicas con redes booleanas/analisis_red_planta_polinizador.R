# T�tulo: Comunidades ecologicas con redes booleanas

# Nombre: Biolog�a de sistemas

# Fecha: Febreo 2021
#######################################################################################

# cargar librer�as
library(BoolNet)

# establecer directorio de trabajo
setwd("C:/Users/Elisa/Dropbox/Docencia adscripci�n biom�dicas/EXAMENES Y RESPUESTAS A PR�CTICAS PCBIOL 2021 1/ensamblaje de comunidades")

# cargar red
net <- loadNetwork("Red_planta_polinizador.txt")
net

# mostrar la red
plotNetworkWiring(net)

# obtener atractores
attr <- getAttractors(net)
attr

# mostrar atractores
plotAttractors(attr)
plotStateGraph(attr) 

# empezamos en (0,1,1,1,1) y preguntamos a qu� atractor vamos a llegar
path <- getPathToAttractor(net, c(0,1,1,1,1))
path

plotSequence(sequence=path)
