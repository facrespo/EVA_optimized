setwd("E:/investigacion/Ivan_araya/estructura_DEUDA");

library(openxlsx);
library(lpSolve);
library(plot3D);
library(tseries);
library(dplyr);
library(tidyverse);
library(alabama);
library(pracma);
library(nloptr);

archivo <- read.xlsx("BD_EVA_16E.xlsx", sheet="BD_EVA_16E", colNames = TRUE);

decision<-1;

archivo$assets <- archivo$Deudat1 + archivo$Patrimoniot1;

nombre_empresa="CAROZZI S.A.";

archivoi <- archivo[archivo$Empresa==nombre_empresa,];

nn <- dim(archivoi);

T <- 0.27;
Nint<-25;

theta1 <- 3.0;
theta2 <- 3.0;
theta <- 10000;
constante <- 1000000;
epsilon <- 0.05;
alpha <- 0.95; #nivel de confianza del intervalo
c1 <- 0.05;
cm1 <- 0.05;
xfx <- function(x,m,s) (x*dnorm(x, mean = m, s = s, log = FALSE));
fx <- function(x,m,s) (dnorm(x, mean = m, s = s, log = FALSE));

kd=seq(-1*Nint/2,1*Nint/2)/100;
ke=seq(-1*Nint/2,1*Nint/2)/100;

Rdm <- median(archivoi[1:(nn[1]-1)]$Rendimientot);
limitesa <- range(archivoi[1:(nn[1]-1)]$Alfa);
limitesb <- range(archivoi[1:(nn[1]-1)]$Beta);
limitesa <- limitesa-(Rdm*(1-T))
limitesb <- limitesb-(Rdm*(1-T))
limitesa <- limitesa/(2*T-1);
limitesb <- -1*limitesb;
rangea=limitesa[1]-limitesa[2];
rangeb=limitesb[1]-limitesb[2];

kd=as.data.frame(limitesa[2]+((0:(2*Nint+1))/(2*Nint))*rangea);
ke=as.data.frame(limitesb[2]+((0:(2*Nint+1))/(2*Nint))*rangeb);


Tasas <- zeros((2*Nint+2)^2,4);
l<-0;
for (i in 1:(2*Nint+2)){
   for (j in 1:(2*Nint+2)){
      Tasas[l+j,1]=kd[i,1];
      Tasas[l+j,2]=ke[j,1];
      Tasas[l+j,3]=Rdm*(1-T)+kd[i,1]*(2*T-1);
      Tasas[l+j,4]=Rdm*(1-T)-ke[j,1];
   }
   l<-(l+j);
}
Ntot<-l;

colnames(Tasas) <- c("kd","ke","Alpha","Beta");

ndim <- dim(archivoi)
DP <- cbind(archivoi$Deudat1,archivoi$Patrimoniot1);
AC <- archivoi$assets[ndim[1]-1];
Am <- mean(archivoi$assets[1:(ndim[1]-1)]);
As <- sd(archivoi$assets[1:(ndim[1]-1)]);
archivo_opt <- NULL;
for (i in 1:Ntot){
   R <- cbind(Tasas[i,3],Tasas[i,4]);
   if (decision==0){
      l <- cbind(AC,DP[ndim[1],1],DP[ndim[1],2],R[1],R[2],archivoi$EVAtotal[ndim[1]],archivoi$EVAtotal[ndim[1]]/AC,0,0,0,0,0,0,0,0,0,0,0,0,0,0);
      l <- as.data.frame(l);
      colnames(l) <- c("assets(n-1)","D(n)","P(n)","alpha(n-1)","beta(n-1)","EVA(n)","EVA(n)/assets(n-1)","EVA_OP1_3.0_3.0","EVA_OP1/A(n-1)_3.0_3.0","D_OP1_3.0_3.0","P_OP1_3.0_3.0","A_ORP3_3.0_3.0","D_ORP3_3.0_3.0","P_ORP3_3.0_3.0","EVA_ORP3_3.0_3.0","EVA/A_ORP3_3.0_3.0","A_ORP2_3.0_3.0","D_ORP2_3.0_3.0","P_ORP2_3.0_3.0","EVA_ORP2_3.0_3.0","EVA/A_ORP2_3.0_3.0");
   }
   else {
   A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
   b1 <- c(AC, 0, 0);
   dir1 <- c("==","<=", "<=");
   coef <- R;
   solucion <- lp('max', coef, A1, dir1, b1);
   l <- cbind(AC,DP[ndim[1],1],DP[ndim[1],2],R[1],R[2],archivoi$EVAtotal[ndim[1]],archivoi$EVAtotal[ndim[1]]/AC,solucion$objval,(solucion$objval)/AC);
   l <- as.data.frame(l);
   l$DOPT <- solucion$solution[1];
   l$POPT <- solucion$solution[2];
   A = Am - qnorm(alpha)*As;
   A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
   b1 <- c(A, 0, 0);
   dir1 <- c("==","<=", "<=");
   solucion1 <- lp('max', coef, A1, dir1, b1);
   l$AR <- A;
   l$DOPTR <- solucion1$solution[1];
   l$POPTR <- solucion1$solution[2];
   l$EVAR <- solucion1$objval;
   l$Lambda1R <- solucion1$objval/A;
   xfxx <- function(x) (xfx(x,m=Am,s=As));
   fxx <- function(x) (fx(x,m=Am,s=As));
   pedazo <- integral(fxx,10,Inf);
   fn <- function(x) (-1*((R[1]*x[1] + R[2]*x[2]) - (c1+cm1)*(integral(xfxx,x[1]+x[2],Inf)+(x[1]+x[2])*integral(fxx,x[1]+x[2],Inf)) - cm1*(x[1]+x[2]-Am)));
   #Equality constraints
   eval_g_eq <- function(x)
   {
      return ( (x[1] + x[2])*(1-integral(fxx,x[1]+x[2],Inf))-Am+integral(xfxx,x[1]+x[2],Inf))
   }
   # Inequality constraints
   eval_g_ineq <- function (x) {
      constr <- c(x[1] - theta1*x[2],
                  -x[1]+theta2*x[2])
      return (constr)
   }
   # Lower and upper bounds
   lb <- c(0, 0)
   ub <- c(AC, AC)
   # Initial values
   x0 <- c((AC*theta2)/(1+theta2), (AC)/(1+theta2));
   opts <- list( "algorithm"
                 = "NLOPT_LN_COBYLA", #NLOPT_GN_ISRES
                 "xtol_rel"
                 = 1.0e-5,
                 "maxeval"= 10000,
                 "tol_constraints_ineq" = rep( 1.0e-10, 2 ));
   res <- nloptr(
      x0          = x0,
      eval_f      = fn,
      lb          = lb,
      ub          = ub,
      eval_g_ineq = eval_g_ineq,
      opts        = opts )
   print(res)
   l$DOPTR1 <- res$solution[1];
   l$POPTR1 <- res$solution[2];
   l$APTR1 <- l$DOPTR1+l$POPTR1;
   l$EVATR1 <- R[1]*res$solution[1]+R[2]*res$solution[2];
   l$Lambda1TR1 <- l$EVATR1/l$APTR1;
   colnames(l) <- c("assets(n-1)","D(n)","P(n)","alpha(n-1)","beta(n-1)","EVA(n)","EVA(n)/assets(n-1)","EVA_OP1_3.0_3.0","EVA_OP1/A(n-1)_3.0_3.0","D_OP1_3.0_3.0","P_OP1_3.0_3.0","A_ORP3_3.0_3.0","D_ORP3_3.0_3.0","P_ORP3_3.0_3.0","EVA_ORP3_3.0_3.0","EVA/A_ORP3_3.0_3.0","D_ORP2_3.0_3.0","P_ORP2_3.0_3.0","A_ORP2_3.0_3.0","EVA_ORP2_3.0_3.0","EVA/A_ORP2_3.0_3.0");
   }
   archivo_opt <- as.data.frame(archivo_opt);
   archivo_opt <- rbind(archivo_opt,l);
}

#colnames(archivo_opt) <- c("Empresas","assets(n-1)","D(n)","P(n)","alpha(n-1)","beta(n-1)","EVA(n)","EVA(n)/assets(n-1)","EVA_OP1_3.0_3.0","EVA_OP1/A(n-1)_3.0_3.0","D_OP1_3.0_3.0","P_OP1_3.0_3.0","A_ORP3_3.0_3.0","D_ORP3_3.0_3.0","P_ORP3_3.0_3.0","EVA_ORP3_3.0_3.0","EVA/A_ORP3_3.0_3.0","A_ORP2_3.0_3.0","D_ORP2_3.0_3.0","P_ORP2_3.0_3.0","EVA_ORP2_3.0_3.0","EVA/A_ORP2_3.0_3.0");
archivo_opt <- as.data.frame(archivo_opt);
archivo_opt <- cbind(Tasas,archivo_opt);

theta1 <- 3.0;
theta2 <- 0;

archivo_opt2 <- NULL;
for (i in 1:Ntot){
   R <- cbind(Tasas[i,3],Tasas[i,4]);
   if (decision==0){
      l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0);
      l <- as.data.frame(l);
      colnames(l) <- c("EVA_OP1_3.0_0","EVA_OP1/A(n-1)_3.0_0","D_OP1_3.0_0","P_OP1_3.0_0","A_ORP3_3.0_0","D_ORP3_3.0_0","P_ORP3_3.0_0","EVA_ORP3_3.0_0","EVA/A_ORP3_3.0_0","D_ORP2_3.0_0","P_ORP2_3.0_0","A_ORP2_3.0_0","EVA_ORP2_3.0_0","EVA/A_ORP2_3.0_0");
      
   }
   else {
      A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
      b1 <- c(AC, 0, 0);
      dir1 <- c("==","<=", "<=");
      coef <- R;
      solucion <- lp('max', coef, A1, dir1, b1);
      l <- cbind(solucion$objval,(solucion$objval)/AC);
      l$DOPT <- solucion$solution[1];
      l$POPT <- solucion$solution[2];
      l <- as.data.frame(l);
      A = Am - qnorm(alpha)*As;
      A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
      b1 <- c(A, 0, 0);
      dir1 <- c("==","<=", "<=");
      solucion1 <- lp('max', coef, A1, dir1, b1);
      l$AR <- A;
      l$DOPTR <- solucion1$solution[1];
      l$POPTR <- solucion1$solution[2];
      l$EVAR <- solucion1$objval;
      l$Lambda1R <- solucion1$objval/A;
      xfxx <- function(x) (xfx(x,m=Am,s=As));
      fxx <- function(x) (fx(x,m=Am,s=As));
      pedazo <- integral(fxx,10,Inf);
      fn <- function(x) (-1*((R[1]*x[1] + R[2]*x[2]) - (c1+cm1)*(integral(xfxx,x[1]+x[2],Inf)+(x[1]+x[2])*integral(fxx,x[1]+x[2],Inf)) - cm1*(x[1]+x[2]-Am)));
      #Equality constraints
      eval_g_eq <- function(x)
      {
         return ( (x[1] + x[2])*(1-integral(fxx,x[1]+x[2],Inf))-Am+integral(xfxx,x[1]+x[2],Inf))
      }
      # Inequality constraints
      eval_g_ineq <- function (x) {
         constr <- c(x[1] - theta1*x[2],
                     -x[1]+theta2*x[2],
                     x[1]+x[2]-AC)
         return (constr)
      }
      # Lower and upper bounds
      lb <- c(0, 0);
      ub <- c(AC, AC);
      # Initial values
      x0 <- c((AC*theta2)/(1+theta2), (AC)/(1+theta2));
      opts <- list( "algorithm"
                    = "NLOPT_LN_COBYLA",
                    "xtol_rel"
                    = 1.0e-5,
                    "maxeval"= 5000,
                    "tol_constraints_ineq" = rep( 1.0e-10, 3));
      res <- nloptr(x0          = x0,
         eval_f      = fn,
         lb          = lb,
         ub          = ub,
         eval_g_ineq = eval_g_ineq,
         opts        = opts );
      print(res)
      l$DOPTR1 <- res$solution[1];
      l$POPTR1 <- res$solution[2];
      l$APTR1 <- res$solution[1]+res$solution[2];
      l$EVATR1 <- R[1]*res$solution[1]+R[2]*res$solution[2];
      l$Lambda1TR1 <- l$EVATR1/l$APTR1;
      colnames(l) <- c("EVA_OP1_3.0_0","EVA_OP1/A(n-1)_3.0_0","D_OP1_3.0_0","P_OP1_3.0_0","A_ORP3_3.0_0","D_ORP3_3.0_0","P_ORP3_3.0_0","EVA_ORP3_3.0_0","EVA/A_ORP3_3.0_0","D_ORP2_3.0_0","P_ORP2_3.0_0","A_ORP2_3.0_0","EVA_ORP2_3.0_0","EVA/A_ORP2_3.0_0");
   }
   archivo_opt2 <- rbind(archivo_opt2,l);
}

#colnames(archivo_opt2) <- c("EVA_OP1_3.0_0","EVA_OP1/A(n-1)_3.0_0","D_OP1_3.0_0","P_OP1_3.0_0","A_ORP3_3.0_0","D_ORP3_3.0_0","P_ORP3_3.0_0","EVA_ORP3_3.0_0","EVA/A_ORP3_3.0_0","A_ORP2_3.0_0","D_ORP2_3.0_0","P_ORP2_3.0_0","EVA_ORP2_3.0_0","EVA/A_ORP2_3.0_0");
archivo_opt2 <- as.data.frame(archivo_opt2);


archivo_opt <- cbind(archivo_opt, archivo_opt2);

theta1 <- 3.0;
theta2 <- 1.5;

archivo_opt3 <- NULL;
for (i in 1:Ntot){
   R <- cbind(Tasas[i,3],Tasas[i,4]);
   if (decision==0){
      l <- cbind(0,0,0,0,0,0,0,0,0,0,0,0,0,0);
      l <- as.data.frame(l);
      colnames(l) <- c("EVA_OP1_3.0_1.5","EVA_OP1/A(n-1)_3.0_1.5","D_OP1_3.0_1.5","P_OP1_3.0_1.5","A_ORP3_3.0_1.5","D_ORP3_3.0_1.5","P_ORP3_3.0_1.5","EVA_ORP3_3.0_1.5","EVA/A_ORP3_3.0_1.5","D_ORP2_3.0_1.5","P_ORP2_3.0_1.5","A_ORP2_3.0_1.5","EVA_ORP2_3.0_1.5","EVA/A_ORP2_3.0_1.5");
      
   }
   else {
      A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
      b1 <- c(AC, 0, 0);
      dir1 <- c("==","<=", "<=");
      coef <- R;
      solucion <- lp('max', coef, A1, dir1, b1);
      l <- cbind(solucion$objval,(solucion$objval)/AC);
      l$DOPT <- solucion$solution[1];
      l$POPT <- solucion$solution[2];
      l <- as.data.frame(l);
      A = Am - qnorm(alpha)*As;
      A1 <- matrix(c(1, 1, -1, 1, -1*theta1, theta2), ncol=2);
      b1 <- c(A, 0, 0);
      dir1 <- c("==","<=", "<=");
      solucion1 <- lp('max', coef, A1, dir1, b1);
      l$AR <- A;
      l$DOPTR <- solucion1$solution[1];
      l$POPTR <- solucion1$solution[2];
      l$EVAR <- solucion1$objval;
      l$Lambda1R <- solucion1$objval/A;
      xfxx <- function(x) (xfx(x,m=Am,s=As));
      fxx <- function(x) (fx(x,m=Am,s=As));
      pedazo <- integral(fxx,10,Inf);
      fn <- function(x)  (-1*((R[1]*x[1] + R[2]*x[2]) - (c1+cm1)*(integral(xfxx,x[1]+x[2],Inf)+(x[1]+x[2])*integral(fxx,x[1]+x[2],Inf)) - cm1*(x[1]+x[2]-Am)));
      #Equality constraints
      eval_g_eq <- function(x)
      {
         return ( (x[1] + x[2])*(1-integral(fxx,x[1]+x[2],Inf))-Am+integral(xfxx,x[1]+x[2],Inf))
      }
      # Inequality constraints
      eval_g_ineq <- function (x) {
         constr <- c(x[1] - theta1*x[2],
                     -x[1]+theta2*x[2])
         return (constr)
      }
      # Lower and upper bounds
      lb <- c(0, 0)
      ub <- c(AC, AC)
      # Initial values
      x0 <- c((AC*theta2)/(1+theta2), (AC)/(1+theta2));
      opts <- list( "algorithm"
                    = "NLOPT_LN_COBYLA",
                    "xtol_rel"
                    = 1.0e-5,
                    "maxeval"= 10000,
                    "tol_constraints_ineq" = rep( 1.0e-10, 2));
      res <- nloptr(
         x0          = x0,
         eval_f      = fn,
         lb          = lb,
         ub          = ub,
         eval_g_ineq = eval_g_ineq,
         opts        = opts )
      print(res)
      l$DOPTR1 <- res$solution[1];
      l$POPTR1 <- res$solution[2];
      l$APTR1 <- res$solution[1]+res$solution[2];
      l$EVATR1 <- R[1]*res$solution[1]+R[2]*res$solution[2];
      l$Lambda1TR1 <- l$EVATR1/l$APTR1;
      colnames(l) <- c("EVA_OP1_3.0_1.5","EVA_OP1/A(n-1)_3.0_1.5","D_OP1_3.0_1.5","P_OP1_3.0_1.5","A_ORP3_3.0_1.5","D_ORP3_3.0_1.5","P_ORP3_3.0_1.5","EVA_ORP3_3.0_1.5","EVA/A_ORP3_3.0_1.5","D_ORP2_3.0_1.5","P_ORP2_3.0_1.5","A_ORP2_3.0_1.5","EVA_ORP2_3.0_1.5","EVA/A_ORP2_3.0_1.5");
   }
   archivo_opt3 <- rbind(archivo_opt3,l);
}

#colnames(archivo_opt3) <- c("EVA_OP1_3.0_1.5","EVA_OP1/A(n-1)_3.0_1.5","D_OP1_3.0_1.5","P_OP1_3.0_1.5","A_ORP3_3.0_1.5","D_ORP3_3.0_1.5","P_ORP3_3.0_1.5","EVA_ORP3_3.0_1.5","EVA/A_ORP3_3.0_1.5","A_ORP2_3.0_1.5","D_ORP2_3.0_1.5","P_ORP2_3.0_1.5","EVA_ORP2_3.0_1.5","EVA/A_ORP2_3.0_1.5");
archivo_opt3 <- as.data.frame(archivo_opt3);

archivo_opt <- cbind(archivo_opt, archivo_opt3);


write.xlsx(archivo_opt, file="EVA_optimized_assets_random_Carozzi.xlsx", sheetName = "Hoja1", col.names = TRUE, row.names = FALSE, append = FALSE, overwrite = TRUE);

Deuda = as.data.frame(cbind(archivo_opt$Empresas,archivo_opt$D_mean,archivo_opt$D_lineal_opt_3.0_3.0,archivo_opt$D_CVAR_Total_3.0_3.0));
colnames(Deuda) <- c("Empresas","D_mean","D_lineal_3.0_3.0","D_CVAR_3.0_3.0");


barplot(height=as.matrix(Deuda[,2:4]),names=Deuda[,1]);
