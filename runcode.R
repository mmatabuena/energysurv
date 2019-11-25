
library("RcppAlgos")
 

   
library("coin")
   
data("GTSG")
   

Rcpp::sourceCpp('Algorithms.cpp')
    

         
         
         x=GTSG$time[GTSG$group=="Chemotherapy+Radiation"]
         y= GTSG$time[GTSG$group=="Chemotherapy"]
         c1= GTSG$event[GTSG$group=="Chemotherapy+Radiation"]
         c2= GTSG$event[GTSG$group=="Chemotherapy"]
         n1= length(x)
         n2= length(y)
         permutaciones=100000
         mpermut=permuteSample(0:(n1+n2-1), n1+n2, n = permutaciones,  repetition = FALSE, Parallel = TRUE) #
         agregada= c(x[c1==1],y[c2==1])
         agregada2= distanciarma(agregada,agregada)[lower.tri(distanciarma(agregada,agregada))]
         agregada2= agregada2^2
         agregada2= agregada2[agregada2>0]
         agregada= sqrt(median(agregada2))
         a1=kernelgausianogeneral(x,y,as.integer(c1),as.integer(c2),mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)),sigma=agregada)[1]
         a4=kernellaplacioanogeneral(x,y,as.integer(c1),as.integer(c2),mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)),sigma=agregada)[1]
         a5=energiageneral(x,y,c1,c2,mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)))[1
           muestra1= x
           muestra2= y
           indicador1= c1
           indicador2= c2
           data= data.frame(time=  c(muestra1,muestra2),trt= c(indicador1,indicador2),grupo= c(rep(0,length(indicador1)),rep(1,length(indicador2))))
           data[,3]= as.factor(data[,3])
           logrank= pvalue(logrank_test(Surv(time,trt)~grupo, data = data,type="logrank",distribution = "approximate", ties.method = "average-scores"))
           gehan= pvalue(logrank_test(Surv(time,trt)~grupo, data = data,type="Gehan-Breslow",distribution = "approximate", ties.method = "average-scores"))
           Tarone= pvalue(logrank_test(Surv(time,trt)~grupo, data = data,type="Tarone-Ware",distribution = "approximate", ties.method = "average-scores"))
           peto= pvalue(logrank_test(Surv(time,trt)~grupo, data = data,type="Fleming-Harrington",rho=1,gamma=0,distribution = "approximate", ties.method = "average-scores"))
           Fleming=  pvalue(logrank_test(Surv(time,trt)~grupo, data = data,type="Fleming-Harrington",rho=1,gamma=1,distribution = "approximate", ties.method = "average-scores"))
           
           a1
           a4
           a5
           logrank
           gehan
           Tarone
           peto
           Fleming
        


