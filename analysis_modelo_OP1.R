setwd("E:/investigacion/Ivan_araya/estructura_DEUDA");

library(openxlsx);
library(lpSolve);
library(plot3D);
library(tseries);
library(dplyr);
library(tidyverse);
library(scatterplot3d);

archivo <- read.xlsx("BD_EVA_16E.xlsx", sheet="BD_EVA_16E", colNames = TRUE);

archivo$assets <- archivo$Deudat1 + archivo$Patrimoniot1;

nn <- dim(archivo);

theta1 <- 3;
theta2 <- 0;

for (i in 1:nn[1]){
  coef <- c(archivo$Alfa[i],archivo$Beta[i]);
  A <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
  b <- c(archivo$assets[i], 0, 0);
  dir <- c("==","<=", "<=");
  solucion <- lp('max', coef, A, dir, b);
  archivo$EVA_OPT[i] <- solucion$objval;
  archivo$DOPT[i] <- solucion$solution[1];
  archivo$POPT[i] <- solucion$solution[2];
  archivo$LEVEST[i] <- solucion$solution[1]/solucion$solution[2];
  archivo$DvA[i] <- solucion$solution[1]/archivo$assets[i];
  archivo$PvA[i] <- solucion$solution[1]/archivo$assets[i];
  archivo$Lambda1[i] <-  archivo$EVA_OPT[i]/archivo$assets[i];
  archivo$Indicador[i] <-  abs(archivo$Alfa[i])/abs(archivo$Beta[i]);
}

write.xlsx(archivo, file="BD_EVA_Optimized_3_0.xlsx", sheetName = "Hoja1", col.names = TRUE, row.names = FALSE, append = FALSE);

archivo <- read.xlsx("BD_EVA_Optimized_3_0.xlsx", sheet="Hoja1", colNames = TRUE);

par(mar=c(1,1,1,1))
plot(archivo$LEVEST,archivo$Alfa,main="Apalancamiento estimado versus Alfa", xlab="Limite optimizacion");

plot(archivo$LEVEST,archivo$Lambda1,main="Apalancamiento estimado versus Lambda 1", xlab="Limite optimizacion");

plot(archivo$LEVEST,archivo$Beta,main="Apalancamiento estimado versus Beta", xlab="Limite optimizacion");

plot(archivo$LEVEST,archivo$EVA_OPT,main="leverage estimado versus EVA", xlab="Limite optimizacion");

plot(archivo$LeverageDP,archivo$EVAtotal,main="Leverage versus EVA", xlab="Limite optimizacion");

plot(archivo$Alfa,archivo$EVAtotal,main="Alfa versus EVA", xlab="alfa");

plot(archivo$Beta,archivo$EVAtotal,main="Beta versus EVA", xlab="beta");

plot(archivo$Alfa,archivo$LeverageDP,main="Alfa versus Leverage", xlab="alfa");

plot(archivo$Beta,archivo$LeverageDP,main="Beta versus Leverage", xlab="beta");

plot(archivo$LEVEST,archivo$EVAtotal,main="Leverage optimizado versus EVA total", xlab="Leverage optimizado");

plot(archivo$Alfa,archivo$EVA_OPT,main="alfa versus EVA optimo", xlab="alfa");

plot(archivo$Beta,archivo$EVA_OPT,main="beta versus EVA optimo", xlab="beta");


plot(archivo$EVAtotal,archivo$EVA_OPT,main="EVA real versus EVA optimized", xlab="EVA");

boxplot(EVAtotal ~ LEVEST,data=archivo,main="Apalancamiento estimado versus EVA real", xlab="Limite optimizacion");

boxplot(EVAtotal ~ LEVEST,data=archivo,main="Apalancamiento estimado versus EVA real", xlab="Limite optimizacion");



plot(archivo$assets,archivo$DOPT,main="Valor assets versus deuda Optima");

plot(archivo$assets,archivo$POPT,main="Valor assets versus Patrimonio Optimo");

plot(archivo$LEVEST,archivo$assets,main="Valor assets versus Patrimonio Optimo");

plot(archivo$EVAtotal,archivo$EVA_OPT,main="EVA real versus EVA Optimo");
lines(c(-4e+07,8e+07),c(-4e+07,8e+07))

archivo2<-archivo[archivo$LEVEST<=0.1,];

plot(archivo2$Alfa,archivo2$EVAtotal,main="Alfa versus EVA", xlab="alfa");

plot(archivo2$Beta,archivo2$EVAtotal,main="Beta versus EVA", xlab="beta");


plot(archivo$Alfa,archivo$DOPT,main="Alfa versus Monto Deuda Optima");

plot(archivo$DOPT,archivo$EVA_OPT,main="Deuda Optima versus EVA Optimo");

plot(archivo$Beta,archivo$POPT,main="Beta versus Monto Patrimonio Optimo");

plot(archivo$POPT,archivo$EVA_OPT,main="Patrimonio Optimo versus EVA Optimo");

plot(archivo$DvA,archivo$EVA_OPT,main="Deuda div assets Optimo versus EVA Optimo");

plot(archivo$PvA,archivo$EVA_OPT,main="Patrimonio div assets Optimo versus EVA Optimo");

scatter3D(x=archivo$DOPT, y=archivo$POPT, z=archivo$EVA_OPT , pch = ".", col = "red", bty = "f", cex = 2, colkey = FALSE);

surf3D(x=as.matrix(archivo$DOPT), y=as.matrix(archivo$POPT), z=as.matrix(archivo$EVA_OPT), colvar = as.matrix(archivo$EVA_OPT), colkey = FALSE)

scatter3D(x=archivo$DOPT, y=archivo$POPT, z=archivo$EVA_OPT/archivo$assets , pch = ".", col = "red", bty = "f", cex = 2, colkey = FALSE, ticktype = "detailed");

scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVA_OPT/archivo$assets , pch = ".", col = "red", bty = "f", cex = 2, colkey = FALSE, ticktype = "detailed")
par(mar=c(1,1,1,1))
scatter3D(x=archivo$DOPT, y=archivo$POPT, z=archivo$EVA_OPT, pch = 18, cex = 2, theta = 20, phi = 20, ticktype = "detailed",
xlab = "wt", ylab = "disp", zlab = "mpg",  surf = list(x = archivo$DOPT, y = archivo$POPT, z = archivo$EVA_OPT, facets = NA, fit = fitpoints), main = "mtcars")

scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVA_OPT, pch = 18, cex = 2, theta = 20, phi = 20, ticktype = "detailed",
          xlab = "alpha", ylab = "betha", zlab = " ",  surf = list(x = archivo$Alfa, y=archivo$Beta, z = matrix(archivo$EVA_OPT), facets = NA,  fit=archivo$EVA_OPT), main = " ")


scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVAtotal, pch = ".", col = "red", bty = "f", cex = 2, colkey = FALSE, ticktype = "detailed", xlab = "alpha", ylab = "betha")

scatter3D(xs, ys, zs, ticktype = "detailed", pch = 16, 
          sctt3D+    xlim = c(0, 2*pi), ylim = c(0, 0.3), zlim = c(-1.5, 1.5), 
          sctt3D+    CI = CI, theta = 20, phi = 30, cex = 2,
          sctt3D+    surf = list(x = M$x, y = M$y, z = z, border = "black", facets = NA)
          sctt3D+    )

scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVAtotal, pch = ".", col = "red", bty = "f", cex = 2, colkey = FALSE, theta=30, ticktype = "detailed", xlab = "alpha", ylab = "betha")
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVA_OPT, add = TRUE, colkey = FALSE, pch = ".", cex = 3, col = "black")

datos = cbind(archivo$Alfa, archivo$Beta, archivo$EVAtotal)
par(mar=c(1,1,1,1))
s3d <- scatterplot3d(datos, type = "h", color = "blue",
                     angle=55, pch = 16)


par(mar=c(1,1,1,1))
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$Deudat1, pch = ".", col = "gray", bty = "f", cex = 3, colkey = FALSE, theta=30, ticktype = "detailed", xlab = "alpha", ylab = "betha")
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$DOPT, add = TRUE, colkey = FALSE, pch = ".", cex = 3, col = "black")

par(mar=c(1,1,1,1))
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$Patrimoniot1, pch = ".", col = "gray", bty = "f", cex = 3, colkey = FALSE, theta=30, ticktype = "detailed", xlab = "alpha", ylab = "betha")
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$POPT, add = TRUE, colkey = FALSE, pch = ".", cex = 3, col = "black")



x <- archivo$Alfa
y <- archivo$Beta
z <- archivo$EVA_OPT
# Compute the linear regression (z = ax + by + d)
fit <- lm(z ~ x + y)
# predict values on regular xy grid
grid.lines = 26
x.pred <- seq(min(x), max(x), length.out = grid.lines)
y.pred <- seq(min(y), max(y), length.out = grid.lines)
xy <- expand.grid( x = x.pred, y = y.pred)
z.pred <- matrix(predict(fit, newdata = xy), 
                 nrow = grid.lines, ncol = grid.lines)
# fitted points for droplines to surface
fitpoints <- predict(fit)
# scatter plot with regression plane
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVAtotal, pch = 18, cex = 2, 
          theta = 20, phi = 20, ticktype = "detailed",
          xlab = "alpha", ylab = "betha", zlab = "EVA",  
          surf = list(x = x.pred, y = y.pred, z = z.pred,  
                      facets = NA, fit = fitpoints), main = "Comparation")
scatter3D(x=archivo$Alfa, y=archivo$Beta, z=archivo$EVA_OPT, add = TRUE, colkey = FALSE, pch = ".", cex = 3, col = "black")

