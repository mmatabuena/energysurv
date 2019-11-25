#include <RcppArmadillo.h>
// [[Rcpp::depends("RcppArmadillo")]]
#include <omp.h>
// [[Rcpp::plugins(openmp)]]
#include <cmath>
#include <R.h>
#include <Rmath.h>
# include <cmath>
#include <cstdlib>
using namespace Rcpp;
using namespace std;


 // [[Rcpp::export]]
arma::vec pesos( arma::vec  censura){
   double pesosaux;
   int n= censura.size();

   double auxiliar;
   double auxiliar2;
   double auxiliar3;
   double auxiliar4;
   double auxiliar5;
   double auxiliar6;
   arma::vec pesos(n);
   arma::vec pesos2(n);

   pesosaux=1;


   for(int i=0;i<n;i++){
     if(i==0){
       auxiliar= n-i-1;
       auxiliar2= n-i;
       auxiliar3= auxiliar/auxiliar2;
       if(censura(i)==0){
         auxiliar6= 1;

       }else{
         auxiliar6= auxiliar3;

       }


       pesosaux= pesosaux*auxiliar6;
       pesos(i)=  pesosaux;
       pesos2(i)= censura(i)/n;

     }else{

       auxiliar= n-i-1;
       auxiliar2= n-i;
       auxiliar3= auxiliar/auxiliar2;
       if(censura(i)==0){
         auxiliar6= 1;

       }else{
         auxiliar6= auxiliar3;

       }

       pesosaux= pesosaux*auxiliar6;
       auxiliar4= n-i;
       auxiliar5= censura(i)/auxiliar4;
       pesos(i)=  pesosaux;
       pesos2(i)= auxiliar5*pesos(i-1);
     }
   }
   float suma= sum(pesos2);
   return pesos2/suma;

 }


 // [[Rcpp::export]]
arma::mat distanciarma(arma::vec x, arma::vec y){
   int n1= x.size();
   int n2= y.size();
   arma::mat distancias(n1,n2);

   for(int i=0;i<n1;i++){
     for(int j=0;j<n2;j++){
       distancias(i,j)= x(i)-y(j);

     }
   }

   return abs(distancias);
 }





// [[Rcpp::export]]
arma::mat kernelgausiano2(arma::vec x, arma::vec y, float sigma){
  int n1= x.size();
  int n2= y.size();
  arma::mat distancias(n1,n2);
  for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++){
      distancias(i,j)= exp(-pow(((x(i)-y(j))/sigma),2));
        
        

    }
      }
  
  return distancias;
}

// [[Rcpp::export]]
arma::mat kernellaplaciano2
  (arma::vec x, arma::vec y, float sigma){
  int n1= x.size();
  int n2= y.size();
  arma::mat distancias(n1,n2);
  for(int i=0;i<n1;i++){
    for(int j=0;j<n2;j++){
      distancias(i,j)= exp(-(1/sigma)*abs(x(i)-y(j)));
      
    }
  }
  
  return distancias;
}

// [[Rcpp::export]]
arma::mat kernellaplaciano(const arma::mat x, const float sigma) {
  arma::mat z;
  z= -sigma*x;
  z= expmat(z);
  
  return z;
}

// [[Rcpp::export]]
arma::mat kernelgausiano(const arma::mat x, const float sigma) {
  arma::mat z;
  z= x;
  z= z%z;
  z= -sigma*z;
  z= expmat(z);
  return z;
}




 // [[Rcpp::export]]
arma::vec dist_energy(arma::mat  pe1,arma::mat pe2, arma::mat dist1,
                       arma::mat dist2,arma::mat dist3){
   float m=pe1.size();
   float n=pe2.size();
   arma::vec e=((m*n)/(m+n))*(2*trans(pe1)*dist3*pe2-trans(pe1)*dist1*pe1-trans(pe2)*dist2*pe2);
   return e;
 }


// [[Rcpp::export]]
arma::vec kernel(arma::mat  pe1,arma::mat pe2, arma::mat dist1,
                      arma::mat dist2,arma::mat dist3){
  float m=pe1.size();
  float n=pe2.size();
  arma::vec e=((m*n)/(m+n))*(-2*trans(pe1)*dist3*pe2+trans(pe1)*dist1*pe1+trans(pe2)*dist2*pe2);
  return e;
}





 // [[Rcpp::export]]
arma::vec energiaauxiliar(arma::vec x, arma::vec y, arma::vec censura1, arma:: vec censura2){

   arma::uvec  x1indicesordenado= sort_index(x);
   arma::uvec  y1indicesordenado= sort_index(y);
   arma::vec  x1ordenado= x(x1indicesordenado);
   arma::vec  y1ordenado= y(y1indicesordenado);
   arma::vec censura1ordenado= censura1(x1indicesordenado);
   arma::vec censura2ordenado= censura2(y1indicesordenado);
   arma::vec pesos1= pesos(censura1ordenado);
   arma::vec pesos2= pesos(censura2ordenado);
   arma::mat distancia1= distanciarma(x1ordenado,x1ordenado);
   arma::mat distancia2= distanciarma(y1ordenado,y1ordenado);
   arma::mat distancia3= distanciarma(x1ordenado,y1ordenado);

   return dist_energy(pesos1,pesos2,distancia1,distancia2,distancia3);

 }




// [[Rcpp::export]]
arma::vec kernelauxiliargausiano(arma::vec x, arma::vec y, arma::vec censura1, arma:: vec censura2, float sigma=1){
  
  arma::uvec  x1indicesordenado= sort_index(x);
  arma::uvec  y1indicesordenado= sort_index(y);
  arma::vec  x1ordenado= x(x1indicesordenado);
  arma::vec  y1ordenado= y(y1indicesordenado);
  arma::vec censura1ordenado= censura1(x1indicesordenado);
  arma::vec censura2ordenado= censura2(y1indicesordenado);
  arma::vec pesos1= pesos(censura1ordenado);
  arma::vec pesos2= pesos(censura2ordenado);
  arma::mat kernel1= kernelgausiano2(x1ordenado,x1ordenado,sigma);
  arma::mat kernel2= kernelgausiano2(y1ordenado,y1ordenado,sigma);
  arma::mat kernel3= kernelgausiano2(x1ordenado,y1ordenado,sigma);
  
  return  kernel(pesos1,pesos2,kernel1,kernel2,kernel3);

  
}




// [[Rcpp::export]]
arma::vec kernelauxiliarlaplaciano(arma::vec x, arma::vec y, arma::vec censura1, arma:: vec censura2, float sigma=1){
  
  arma::uvec  x1indicesordenado= sort_index(x);
  arma::uvec  y1indicesordenado= sort_index(y);
  arma::vec  x1ordenado= x(x1indicesordenado);
  arma::vec  y1ordenado= y(y1indicesordenado);
  arma::vec censura1ordenado= censura1(x1indicesordenado);
  arma::vec censura2ordenado= censura2(y1indicesordenado);
  arma::vec pesos1= pesos(censura1ordenado);
  arma::vec pesos2= pesos(censura2ordenado);
  arma::mat kernel1= kernellaplaciano2(x1ordenado,x1ordenado,sigma);
  arma::mat kernel2= kernellaplaciano2(y1ordenado,y1ordenado,sigma);
  arma::mat kernel3= kernellaplaciano2(x1ordenado,y1ordenado,sigma);
  
  return  kernel(pesos1,pesos2,kernel1,kernel2,kernel3);
  
  
}




 // [[Rcpp::export]]
 float pvalor(arma::vec permutaciones, arma::vec energia){
   int n= permutaciones.size();
   int contar=0;
   for(int i=0; i<n;i++){
     if(permutaciones(i)>=energia(0)){
       contar= contar+1;
     }
   }
   float n2=n;
   float contar2= contar;

return contar2/(n2+1);

 }


// [[Rcpp::export]]
arma::uvec indices(arma::uvec x, arma:: uvec primeros){
  arma::uvec salida(primeros.size());
  for(unsigned int i=0;i<primeros.size();i++){
    salida(i)= x(primeros(i));
  }
  return salida;
}


// [[Rcpp::export]]
arma::vec indices2(arma::vec x, arma:: uvec primeros){
  arma::vec salida(primeros.size());
  for(unsigned int i=0;i<primeros.size();i++){
    salida(i)= x(primeros(i));
  }
  return salida;
}

 // [[Rcpp::export]]
float energiageneral(const arma::vec x, const arma::vec y, const arma::vec censura1, const arma:: vec censura2, const arma::umat permutaciones,  arma::uvec primeros,  arma::uvec segundos){

   const arma::umat permutaciones2= permutaciones.t();

   int nperm= permutaciones.n_rows;
   arma::uvec  x1indicesordenado= sort_index(x);
   arma::uvec  y1indicesordenado= sort_index(y);
   arma::vec  x1ordenado= x(x1indicesordenado);
   arma::vec  y1ordenado= y(y1indicesordenado);
   arma::vec censura1ordenado= censura1(x1indicesordenado);
   arma::vec censura2ordenado= censura2(y1indicesordenado);
   arma::vec pesos1= pesos(censura1ordenado);
   arma::vec pesos2= pesos(censura2ordenado);
   arma::mat distancia1= distanciarma(x1ordenado,x1ordenado);
   arma::mat distancia2= distanciarma(y1ordenado,y1ordenado);
   arma::mat distancia3= distanciarma(x1ordenado,y1ordenado);

   arma::vec energia;
   energia=dist_energy(pesos1,pesos2,distancia1,distancia2,distancia3);


   arma::uvec peraux;
   arma::uvec peraux2;
   arma::uvec peraux3;
   arma::vec datosconjuntos;
   arma::vec censuraconjuntos;
   arma::vec x1;
   arma::vec x2;
   arma::vec c1;
   arma::vec c2;
   arma::vec resultadospermut(nperm);




   for(int i=0;i<nperm;i++){
      datosconjuntos= join_cols(x,y);
      censuraconjuntos= join_cols(censura1,censura2);
      peraux= permutaciones2.col(i);
      try{
        peraux2= peraux(primeros);
        peraux3= peraux(segundos);

        x1= datosconjuntos(peraux2);
        x2= datosconjuntos(peraux3);
        c1= censuraconjuntos(peraux2);
        c2= censuraconjuntos(peraux3);


      resultadospermut(i)= energiaauxiliar(x1,x2,c1,c2)(0);

      }catch(...){
        resultadospermut(i)= 1000;

      }

   }


   return pvalor(resultadospermut,energia);




    }



// [[Rcpp::export]]
float kernelgausianogeneral(const arma::vec x, const arma::vec y, const arma::vec censura1, const arma:: vec censura2, const arma::umat permutaciones,  arma::uvec primeros,  arma::uvec segundos, float sigma=1){
  
  const arma::umat permutaciones2= permutaciones.t();
  
  int nperm= permutaciones.n_rows;
  arma::uvec  x1indicesordenado= sort_index(x);
  arma::uvec  y1indicesordenado= sort_index(y);
  arma::vec  x1ordenado= x(x1indicesordenado);
  arma::vec  y1ordenado= y(y1indicesordenado);
  arma::vec censura1ordenado= censura1(x1indicesordenado);
  arma::vec censura2ordenado= censura2(y1indicesordenado);
  arma::vec pesos1= pesos(censura1ordenado);
  arma::vec pesos2= pesos(censura2ordenado);
  arma::mat kernel1= kernelgausiano2(x1ordenado,x1ordenado,sigma);
  arma::mat kernel2= kernelgausiano2(y1ordenado,y1ordenado,sigma);
  arma::mat kernel3= kernelgausiano2(x1ordenado,y1ordenado,sigma);
  
  
  
  
  arma::vec energia;
  energia=kernel(pesos1,pesos2,kernel1,kernel2,kernel3);
  
  
  arma::uvec peraux;
  arma::uvec peraux2;
  arma::uvec peraux3;
  arma::vec datosconjuntos;
  arma::vec censuraconjuntos;
  arma::vec x1;
  arma::vec x2;
  arma::vec c1;
  arma::vec c2;
  arma::vec resultadospermut(nperm);
  
  

  
  for(int i=0;i<nperm;i++){

    try{
      
      datosconjuntos= join_cols(x,y);
      censuraconjuntos= join_cols(censura1,censura2);
      peraux= permutaciones2.col(i);
      
      peraux2= peraux(primeros);
      peraux3= peraux(segundos);
      
      x1= datosconjuntos(peraux2);
      x2= datosconjuntos(peraux3);
      c1= censuraconjuntos(peraux2);
      c2= censuraconjuntos(peraux3);
      
      
      resultadospermut(i)= kernelauxiliargausiano(x1,x2,c1,c2,sigma)(0);
      
    }catch(...){
      resultadospermut(i)= 1000;
      
    }
    
  }
  
  
  return pvalor(resultadospermut,energia);
  
}


// [[Rcpp::export]]
float kernellaplacioanogeneral(const arma::vec x, const arma::vec y, const arma::vec censura1, const arma:: vec censura2, const arma::umat permutaciones,  arma::uvec primeros,  arma::uvec segundos, float sigma=1){
  
  const arma::umat permutaciones2= permutaciones.t();
  
  int nperm= permutaciones.n_rows;
  arma::uvec  x1indicesordenado= sort_index(x);
  arma::uvec  y1indicesordenado= sort_index(y);
  arma::vec  x1ordenado= x(x1indicesordenado);
  arma::vec  y1ordenado= y(y1indicesordenado);
  arma::vec censura1ordenado= censura1(x1indicesordenado);
  arma::vec censura2ordenado= censura2(y1indicesordenado);
  arma::vec pesos1= pesos(censura1ordenado);
  arma::vec pesos2= pesos(censura2ordenado);
  
  arma::mat kernel1= kernellaplaciano2(x1ordenado,x1ordenado,sigma);
  arma::mat kernel2= kernellaplaciano2(y1ordenado,y1ordenado,sigma);
  arma::mat kernel3= kernellaplaciano2(x1ordenado,y1ordenado,sigma);
  
  
  
  arma::vec energia;
  energia=kernel(pesos1,pesos2,kernel1,kernel2,kernel3);
  
  
  arma::uvec peraux;
  arma::uvec peraux2;
  arma::uvec peraux3;
  arma::vec datosconjuntos;
  arma::vec censuraconjuntos;
  arma::vec x1;
  arma::vec x2;
  arma::vec c1;
  arma::vec c2;
  arma::vec resultadospermut(nperm);
  
  

  
  for(int i=0;i<nperm;i++){

    try{
      datosconjuntos= join_cols(x,y);
      censuraconjuntos= join_cols(censura1,censura2);
      peraux= permutaciones2.col(i);
      peraux2= peraux(primeros);
      peraux3= peraux(segundos);
      
      x1= datosconjuntos(peraux2);
      x2= datosconjuntos(peraux3);
      c1= censuraconjuntos(peraux2);
      c2= censuraconjuntos(peraux3);
      
      
      resultadospermut(i)= kernelauxiliarlaplaciano(x1,x2,c1,c2,sigma)(0);
      
    }catch(...){
      resultadospermut(i)= 1000;
      
    }
    
  }
  
  
  return pvalor(resultadospermut,energia);
  
  
  
  
}





/*** R

  library("RcppAlgos")
 

   
   library("coin")
   
     data("GTSG")
     
    

         
         
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
         agregada2= agregada^2
         agregada3= agregada2^2
         a1=kernelgausianogeneral(x,y,as.integer(c1),as.integer(c2),mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)),sigma=agregada)[1]
         a4=kernellaplacioanogeneral(x,y,as.integer(c1),as.integer(c2),mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)),sigma=agregada)[1]
         a5=energiageneral(x,y,c1,c2,mpermut,as.integer(0:(n1-1)),as.integer((n1):(n1+n2-1)))[1]
     
       
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
        

*/

