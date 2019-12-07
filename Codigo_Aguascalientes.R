
##### ESTABLECIENDO DIRECTORIO DE TRABAJO #######
### Autores: Vicente Tapia Wenderoth y Jose Luis Tapia Ramirez
#### fecha 06 diciembre 2019
### Proyecciones de poblacion estado de Aguascalientes Mexico

setwd("E:/PROYECTO FINAL ADIII")

###### REMOVIENDO VARIABLES  EXISTENTES #########
rm(list=ls())                  

############ LIBERANDO MEMORIA ##################
gc()                            

###### CARGANDO PAQUETERIAS NECESARIAS ##########
require(forecast)
require(ggplot2)
require(mvtnorm)

###### CARGANDO INSUMOS NECESARIOS #############
load("Aguascalientes.RData")


#### CARGANDO TASAS CENTRALES DE MORTALIDAD ####
mx <- mx [,-c(1,2)]
colnames(mx) <- c(1970:2015)

########## DESAGREGANDO FECUNDIDAD #############

TGF <- 5*colSums(fx5[,-1])

Vx <- fx5[,-1]
FxF <- Vx 
## FxF es la prporcion de FxF
for(i in  2:47){
  Vx[,i-1] <- 5*cumsum(fx5[,i])
  FxF[,i-1] <- Vx[,i-1]/TGF[i-1]
}

x5 <- seq(17.5, 47.5, 5)
Yx.fec <- log(-log(FxF))
Yx.fec.lm <- list()
for(i in 1:46){
  Yx.fec.lm[[i]] <- lm(Yx.fec[-7,i] ~ x5[-7])
}
a <- vector(length = 46)
b <- vector(length = 46)

for(i in 1:46){
  a[i] <- Yx.fec.lm[[i]]$coefficients[1]
  b[i] <- Yx.fec.lm[[i]]$coefficients[2]
}

A <- exp(-exp(a))
B <- exp(b)

x <- c(15:50)
FxF.proy <- matrix(0,36,46)
for(i in 1:46){
  FxF.proy[,i] <- TGF[i]*A[i]^(B[i]^(x))
}
fx <- matrix(0,36,46)
fx[1,] <- FxF.proy[1,]
fx[2:36,] <- FxF.proy[2:36,] - FxF.proy[1:35,]
matplot(fx,type = "l")


colnames(fx)[1:46] <- c(1970:2015)

######## TASAS ESPECIFICAS DE FECUNDIDAD #########

fx <- fx[-36,]
fx <- fx[,11:46]                

######### TASAS ESPECÍFICAS DE MIGRACIÓN #########

ixt.F<-as.matrix(Inmx[Inmx$Sexo=="Mujeres",-c(1:2)])/  
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
ixt.M<-as.matrix(Inmx[Inmx$Sexo=="Hombres",-c(1:2)])/  
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])
ext.F<-as.matrix(Emigx[Emigx$Sexo=="Mujeres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Mujeres",-c(1:2)])
ext.M<-as.matrix(Emigx[Emigx$Sexo=="Hombres",-c(1:2)])/
  as.matrix(Nx[Nx$Sexo=="Hombres",-c(1:2)])

ixt.F<-ixt.F[1:(dim(ixt.F)[1]-30),41:dim(ixt.F)[2]]
ixt.M<-ixt.M[1:(dim(ixt.M)[1]-30),41:dim(ixt.M)[2]]
ext.F<-ext.F[1:(dim(ext.F)[1]-30),41:dim(ext.F)[2]]
ext.M<-ext.M[1:(dim(ext.M)[1]-30),41:dim(ext.M)[2]]


#### DEFINIENDO TIEMPOS, INICIOS Y FINALES ######

edades <- dim(mx) [1] 
edades.fec <- dim(fx)[1]
tiempo.mort <- dim(mx)[2]
añoini.mort <- 1970
añoini.fec <- 1980
tiempo.fec <- dim(fx)[2] 
añobase <- 2015
horizonte <- 35
añofin <- añobase+horizonte
tiempo.tot <- tiempo.mort+horizonte

edades.mig <- dim(ixt.F)[1]
tiempo.mig <- dim(ixt.F)[2]
añoini.mig <- 1995


###### INICIANDO ESTIMACION CON EL METODO #######
########## DE LEE - CARTER (1992) ###############

lc.svd <- function(m,edades,tiempo1,tiempo2,ln){
  if (ln == TRUE){
    lm <- log(m)
  } else{
    lm <- m
  }
  ax <- rowMeans(lm[,tiempo1:tiempo2])
  
  lm_a <- lm - ax
  
  d <- matrix(0, nr = min(edades,tiempo2),
              nc = min(edades,tiempo2))
  
  diag(d) <- svd(lm_a)$d
  
  kt <- (d%*%t(-svd(lm_a)$v))
  bx <- -svd(lm_a)$u
  
  lc.svd <- list(ax = ax, bx = bx, kt = kt, D=d)
  
}

tabmort <- function(m,edades,sex){
  
  mx <- m
  
  nax <- matrix(0.5,dim(mx)[1],dim(mx)[2])
  ## 1 MUJERES 2 HOMBRES
  if(sex==1){
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.01724){
        nax[1,i] <- 0.14903-2.05527*mx[1,i]
      }else if(mx[1,i]>=0.01724 & mx[1,i]<0.06891){
        nax[1,i] <- 0.04667+3.88089*mx[1,i]
      }else{nax[1,i] <- 0.31411}
    }
  }else{
    for(i in 1:dim(mx)[2]){
      if(mx[1,i]<0.023){
        nax[1,i] <- 0.14929-1.99545*mx[1,i]
      }else if(mx[1,i]>=0.023 & mx[1,i]<0.08307){
        nax[1,i] <- 0.02832+3.26021*mx[1,i]
      }else{nax[1,i] <- 0.29915}
    }
  }
  
  
  nax[edades,] <- 1/mx[edades,]
  
  qx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 1:(dim(mx)[1])){
    qx[i,]<-mx[i,]/(1+(1-nax[i,])*mx[i,])
  }
  
  px <- 1-qx
  
  lx<-matrix(1,dim(mx)[1],dim(mx)[2])
  
  for(i in 2:dim(mx)[1]){
    lx[i,] <- lx[i-1,]*px[i-1,]
  }
  
  dx <- matrix(0,dim(mx)[1],dim(mx)[2])
  dx[dim(mx)[1],] <- lx[dim(mx)[1],]
  for(i in 1:(dim(mx)[1]-1)){
    dx[i,]<-lx[i,]-lx[i+1,]
  }
  
  
  Lx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Lx[1,] <- dx[1,]/mx[1,]
  Lx[edades,] <- dx[edades,]/mx[edades,]
  for(i in 2:(edades-1)){
    Lx[i,]<-(lx[i,]+lx[i+1,])/2
  }
  
  Tx<-matrix(0,dim(mx)[1],dim(mx)[2])
  Tx[edades,]<-Lx[edades,]
  for(i in (edades-1):1){
    Tx[i,]<-Lx[i,]+Tx[i+1,]
  }
  
  ex <- Tx/lx
  
  Sx<-matrix(0,(dim(mx)[1]+1),dim(mx)[2])
  Sx[1,]<-Lx[1,]/lx[1,]
  Sx[(edades+1),] <- Tx[edades,]/Tx[(edades-1),]
  for(i in 2:edades){
    Sx[i,]<-Lx[i,]/Lx[i-1,]
  }
  
  tabmort <- list(Edad=c(0:(edades-1)),mx=mx, nax=nax, qx=qx, 
                  px=px, lx=lx, dx=dx, Lx=Lx, Tx=Tx, ex=ex, Sx=Sx)
}


lc.mort <- lc.svd(mx,edades,tiempo1 = 41, 
                  tiempo2=tiempo.mort, ln=TRUE)
lc.fec <- lc.svd(fx, edades = edades.fec,
                 tiempo1 = 35, 
                 tiempo2 = tiempo.fec,
                 ln = TRUE)
lc.inmF <- lc.svd(ixt.F, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)
lc.inmM <- lc.svd(ixt.M, edades = edades.mig,
                  tiempo1 = 1,
                  tiempo2 = tiempo.mig,
                  ln = TRUE)
lc.emigF <- lc.svd(ext.F, edades = edades.mig,
                   tiempo1 = 1,
                   tiempo2 = tiempo.mig,
                   ln = TRUE)
lc.emigM <- lc.svd(ext.M, edades = edades.mig,
                   tiempo1 = 1,
                   tiempo2 = tiempo.mig,
                   ln = TRUE)

dim(lc.mort$kt)


kt1.fit <- auto.arima(lc.mort$kt[1,], trace = TRUE, d=1)

ft1.fit <- auto.arima(lc.fec$kt[1,], trace = TRUE, d=1, allowdrift = T)

it1F.fit <- auto.arima(lc.inmF$kt[1,], trace = TRUE , allowdrift = F)

it1M.fit <- auto.arima(lc.inmM$kt[1,], trace = TRUE , allowdrift = F)

et1F.fit <- auto.arima(lc.emigF$kt[1,], trace = TRUE , allowdrift = F)

et1M.fit <- auto.arima(lc.emigM$kt[1,], trace = TRUE , allowdrift = F)

### ESTABLECIENDO HORIZONTE DE PROYECCIÓN #########

h<- 35


##### PROYECTANDO LOS COMPONENTES POR SEPARADO ####
kt.for <- forecast(kt1.fit, h = h, c(95))
ft.for <- forecast(ft1.fit, h = h, c(95))
itF.for <- forecast(it1F.fit, h = h, c(95))
itM.for <- forecast(it1M.fit, h = h, c(95))
etF.for <- forecast(et1F.fit, h = h, c(95))
etM.for <- forecast(et1M.fit, h = h, c(95))


####### HACIENDO EL PROCESO DETERMINISTA ###########

mx.for <- exp(lc.mort$ax + lc.mort$bx[,1]%*%t(kt.for$mean))  
SxF.for <- tabmort(mx.for[111:220,], edades = 110, sex = 1)$Sx
ExF.for <- tabmort(mx.for[111:220,], edades = 110, sex = 1)$ex
colnames(SxF.for) <- c(2016:2050)
SxM.for <- tabmort(mx.for[1:110,], edades = 110, sex = 2)$Sx 
colnames(SxM.for) <- c(2016:2050)

fx.for <- exp(lc.fec$ax + lc.fec$bx[,1]%*%t(ft.for$mean))

colnames(fx.for) <- c(2016:2050)
colSums(fx.for)


### CALCULANDO TASAS ESPECIFICAS PARA MIGRANTES #####

ixF.for <- rbind(exp(lc.inmF$ax + lc.inmF$bx[,1]%*%t(itF.for$mean)),
                 matrix(0,30,35))
ixM.for <- rbind(exp(lc.inmM$ax + lc.inmM$bx[,1]%*%t(itM.for$mean)),
                 matrix(0,30,35))
exF.for <- rbind(exp(lc.emigF$ax + lc.emigF$bx[,1]%*%t(etF.for$mean)),
                 matrix(0,30,35))
exM.for <- rbind(exp(lc.emigM$ax + lc.emigM$bx[,1]%*%t(etM.for$mean)),
                 matrix(0,30,35))


##### PROYECCIÓN POR EL METODO DE LAS COMPONENTES ####

PxF <- Px[Px$Sexo=="Mujeres", -c(1,2)] ### poblacion de mujeres conciliadas a mitad de año
PxM <- Px[Px$Sexo=="Hombres", -c(1,2)] ### poblacion de hombres conciliadas a mitad de año

NxF <- Nx[Nx$Sexo=="Mujeres", -c(1,2)] 
NxM <- Nx[Nx$Sexo=="Hombres", -c(1,2)] 

PxF.for <- matrix(0,110,36)
PxM.for <- matrix(0,110,36)

NxF.for <- matrix(0,110,36)
NxM.for <- matrix(0,110,36)

PxF.for[,1] <- PxF[,"2016"]
PxM.for[,1] <- PxM[,"2016"]

NxF.for[,1] <- NxF[,"2015"]
NxM.for[,1] <- NxM[,"2015"]

Bx <- matrix(0,35,35)
BF <- vector(length = 35) #### nacimientos totales
BM <- vector(length = 35) #### nacimientos totales masculinos

##### PROYECCIONES DE LAS MUJERES #################

for(i in 2:36){

  ### Ahora con las a mitad de año para edades intermedias de 1 a 108 años
  
  PxF.for[2:109,i] <- (PxF.for[1:108,i-1] +
                      0.5*NxF.for[1:108,i-1]*ixF.for[1:108,i-1]) * SxF.for[1:108,i-1]+
                      NxF.for[2:109,i-1]*0.5*ixF.for[2:109,i-1] -                               
                      NxF.for[1:108,i-1]*exF.for[1:108,i-1]
  
  
  #### Calculando el ultimo grupo de edad grupo abierto
  
  PxF.for[110,i] <- (PxF.for[109,i-1]+
                    0.5*NxF.for[109,i-1] * ixF.for[109,i-1])* SxF.for[109,i-1]-
                    NxF.for[109,i-1]*exF.for[109,i-1]+
                    (PxF.for[110,i-1]+
                    NxF.for[110,i-1]*0.5*ixF.for[110,i-1]) * SxF.for[110,i-1]+
                    NxF.for[110,i-1]*0.5*ixF.for[110,i-1]-
                    NxF.for[110,i-1]*exF.for[110,i-1]
  
  
  ### Calculando los nacimientos
  
  Bx[,i-1] <- fx.for[,i-1]*(PxF.for[16:50,i-1]+
                              0.5*NxF.for[16:50,i-1]*ixF.for[16:50,i-1]+
                              PxF.for[16:50,i])*0.5
  
  ### Calculando los nacimientos de mujeres
  
  BF[i-1] <- (1/2.05)*sum(Bx[,i-1])
  
  ### Calculando el primer grupo de edad
  
  PxF.for[1,i] <- BF[1]*SxF.for[1,i-1]+
    NxF.for[1,i-1]*0.5*ixF.for[1,i-1]-
    NxF.for[1,i-1]*exF.for[1,i-1]
  
  ## Poblacion a mitad de año
  NxF.for[,i] <- 0.5*(PxF.for[,i-1]+PxF.for[,i])
}

matplot(NxF.for,type="l")

##### PROYECCIONES DE LOS HOMBRES #################

for(i in 2:36){
  
  ### ahora con las a mitad de año edades intermedias 1 a 108 años
  PxM.for[2:109,i] <- (PxM.for[1:108,i-1] +
                      0.5*NxM.for[1:108,i-1]*ixM.for[1:108,i-1]) * SxM.for[1:108,i-1]+
                      NxM.for[2:109,i-1]*0.5*ixM.for[2:109,i-1] -                               
                      NxM.for[1:108,i-1]*exM.for[1:108,i-1]
  
  ####como calcular el ultimo gpo de edad grupo abierto
  
  PxM.for[110,i] <- (PxM.for[109,i-1]+
                    0.5*NxM.for[109,i-1] * ixM.for[109,i-1])* SxM.for[109,i-1]-
                    NxM.for[109,i-1]*exM.for[109,i-1]+
                    (PxM.for[110,i-1]+
       NxM.for[110,i-1]*0.5*ixM.for[110,i-1]) * SxM.for[110,i-1]+
    NxM.for[110,i-1]*0.5*ixM.for[110,i-1]-
    NxM.for[110,i-1]*exM.for[110,i-1]
 
  BM[i-1] <- (1.05/2.05)*sum(Bx[,i-1])
  ###PRIMER GRUPO DE EDAD
  PxM.for[1,i] <- BM[1]*SxM.for[1,i-1]+
    NxM.for[1,i-1]*0.5*ixM.for[1,i-1]-
    NxM.for[1,i-1]*exM.for[1,i-1]
  ## poblacion a mitad de año
  NxM.for[,i] <- 0.5*(PxM.for[,i-1]+PxM.for[,i])
}

PxF.for <- round(PxF.for,0)
PxM.for <- round(PxF.for,0)

matplot(NxM.for,type="l")


### Indicadores y totales

### Población de hombres ####

colSums(NxM.for)

### Poblacion de mujeres ####
colSums(NxF.for)

### Población total ######
colSums(NxM.for)+ colSums(NxF.for)

### Tasa global de fecundidad ###
colSums(fx.for)

# grafica TGF
v1 <- 1:35
v2 <- 2016:2050

matplot(colSums(fx.for), type = "l", xlab = "Año", ylab = "Tasa global de fecundidad", las=1,xaxt="n",
        main="Tasa global de fecundidad en Aguascalientes:  2016-2050")
axis(side = 1, at = v1, labels = v2,) 


# grafica Tasas especificas de fecundidad

v3 <- 1:35
v4 <- 15:49

matplot(fx.for, type = "l", xlab = "Edad", ylab = "Tasa especifica de fecundidad", las=1,xaxt="n",
        main="Tasas específicas de fecundidad por edades simples en Aguascalientes:  2015-2050") 
        axis(side = 1, at = v3, labels = v4,)    


### considerando migracion interna (contante)

## hombres
InmigrantesH <- NxM.for*IniH
EmigrantesH <- NxM.for*EmiH

pob_total_hombres <- NxM.for+InmigrantesH-EmigrantesH

colSums(pob_total_hombres)
colSums(NxM.for)

## mujeres
InmigrantesM <- NxF.for*IniM
EmigrantesM <- NxF.for*EmiM

pob_total_mujeres <- NxF.for+InmigrantesM-EmigrantesM

colSums(pob_total_mujeres)
colSums(NxF.for)

### total población

pob_tot <- colSums(pob_total_hombres)+colSums(pob_total_mujeres)
colSums(pob_total_hombres)+colSums(pob_total_mujeres)

### grafica de poblacion total ###

v5 <- 1:36
v6 <- 2015:2050

matplot((colSums(pob_total_hombres)+colSums(pob_total_mujeres)), type = "l", xlab = "Año", ylab = "", las=1,xaxt="n",
        main="Población total proyectada:  2016-2050") 
axis(side = 1, at = v5, labels = v6,)     

##### GRÁFICA EZPERANZA DE VIDA



e0xF.for <- tabmort(mx.for[111:220,], edades=110, sex=1)$ex
e0xM.for <- tabmort(mx.for[1:110,], edades=110, sex=2)$ex

colnames(e0xF.for) <- c(2016:2050)
colnames(e0xM.for) <- c(2016:2050)

e0xF.for[1,]
e0xM.for[1,]

matplot(e0xF.for[1,],type="l")

TGF.dat<-data.frame(year= c(1970:2050),
                    mean= c(TGF, colSums(fx.for)))

ggplot(TGF.dat, aes(x=year, y=mean))+  geom_line(aes(y=mean), col= "blue")


v7 <- 1:36
v8 <- 2015:2050

matplot(cbind(colSums(pob_total_hombres),colSums(pob_total_mujeres)),type="l",xaxt = 'n',
        main = "Proyección población total por sexo 2016-2050",xlab="años",ylab="población",col=c("black","blue"),lty=c(3,3),cex.main=1)
legend(1, 10200, legend=c("Hombres", "Mujeres"),
       col=c("black", "blue"), lty=1:2, cex=0.8)
axis(side=1, at = v7, labels = v8,)
     
qxF.for <- tabmort(mx.for[111:220,], edades=110, sex=1)$qx
qxM.for <- tabmort(mx.for[1:110,], edades=110, sex=2)$qx

###################################################################
##################################################################
##################### ESPERANZA DE VIDA ###########################

mx2<-data.matrix(mx)

e0xF <- tabmort(mx2[111:220,], edades=110, sex=1)$ex
e0xM <- tabmort(mx2[1:110,], edades=110, sex=2)$ex
e0xF.for <- tabmort(mx.for[111:220,], edades=110, sex=1)$ex
e0xM.for <- tabmort(mx.for[1:110,], edades=110, sex=2)$ex



colnames(e0xF.for) <- c(2016:2050)
colnames(e0xM.for) <- c(2016:2050)

e0xF.for[1,]
e0xM.for[1,]



E0F.dat<-data.frame(year= c(1970:2050),
                    mean= c(e0xF[1,], e0xF.for[1,]))
E0M.dat<-data.frame(year= c(1970:2050),
                    mean= c(e0xM[1,], e0xM.for[1,]))



ggplot(E0F.dat, aes(x=year, y=mean))+  
  geom_line(aes(y=mean), col= "blue")+
  labs (title = "Esperanza de vida al nacimiento mujeres, 1970-2050",  x = "Esperanza de vida al nacimiento", y = "año")


ggplot(E0M.dat, aes(x=year, y=mean))+  
  geom_line(aes(y=mean), col= "blue")+
  labs (title = "Esperanza de vida al nacimiento hombres, 1970-2050",  x = "Esperanza de vida al nacimiento", y = "año")





