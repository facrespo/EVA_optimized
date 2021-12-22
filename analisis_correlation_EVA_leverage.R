setwd("E:/investigacion/Ivan_araya/estructura_DEUDA");

library(xlsx);
library(lpSolve);
library(plot3D);
library(tseries);
library(dplyr);
library(tidyverse);
library(alabama);
library(pracma);
library(nloptr);

archivo <- read.xlsx("BD_EVA_16E.xlsx", 1, header=TRUE, colClasses=NA);

archivo$activo <- archivo$Deudat1 + archivo$Patrimoniot1;

nn <- dim(archivo);

nombre_empresas <- distinct(archivo, archivo$Empresa);

ne <- dim(nombre_empresas);

archivo_opt <- NULL;
for (i in 1:ne[1]){
   archivoi <- archivo[archivo$Empresa==nombre_empresas[i,1],];
   ndim <- dim(archivoi)
   r<-cor(archivoi$EVAtotal,archivoi$LeverageDP)
   test<-cor.test(archivoi$EVAtotal,archivoi$LeverageDP)
   l <- cbind(nombre_empresas[i,1],r)
   l[3] <- test$p.value;
   archivo_opt <- as.data.frame(archivo_opt);
   archivo_opt <- rbind(archivo_opt,l);
}
colnames(archivo_opt) <- c("Empresas","Correlation EVA(n) D(n)","p-value");

archivo_opt <- as.data.frame(archivo_opt);

write.xlsx(archivo_opt, file="EVA_Leverage_Correlation.xlsx", sheetName = "Hoja1", col.names = TRUE, row.names = FALSE, append = FALSE);
