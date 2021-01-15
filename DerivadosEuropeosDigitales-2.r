rm(list=ls())

#Integrantes del equipo
#Hernández García Yesenia Inés
#Munguía Landín Luis
#Ramírez Guízar Brenda Jazmín
#Reyes Martínez Samuel Joshua

#Código para valuar derivados europeos y digitales

'Acomodar matriz'
mover <- function(A){
  B= A
  n=sqrt(length(A))
  B[1,1]=A[1,1]
  for (i in 2:n){
    for (j in 1:i){
      B[i,j]= A[i,(i+1-j)]
    } 
  }
  
  return(B)
}

'Arbol precio'
construyearbolpreciossub <- function(S1,u,d,N) {
  tree = matrix(0, nrow=N+1, ncol=N+1)
  for (i in 1:(N+1)) {
    for (j in 1:i) {
      tree[i, j] = S1 * u^(j-1) * d^((i-1)-(j-1))
    }  }
  return(tree)
}

'Proba q'
q_prob <- function(r,u, d,N,TTT) {
  delta_t=TTT/N
  return((exp(r*delta_t) - d)/(u-d))
}


'Arbol derivado'
valuarderivado <- function(tree,S1,r,u,d,N,TTT, K, type,valordigital) {
  q = q_prob(r,u,d,N,TTT)
  delta_t=TTT/N 
  option_tree = matrix(0, nrow=nrow(tree), ncol=ncol(tree))
  if(type == 'put europea plain vanilla') {
    option_tree[nrow(option_tree),] = pmax(K - tree[nrow(tree),], 0)
  } else if (type=='call europea plain vanilla'){
    option_tree[nrow(option_tree),] = pmax(tree[nrow(tree),] - K, 0)
  } else if(type=='forward'){
    option_tree[nrow(option_tree),] = tree[nrow(tree),] - K
  }else if(type=='call digital'){
    for(k in 1:nrow(tree)){
      if(tree[nrow(tree),k]< valordigital){
        option_tree[nrow(tree),k]<- valordigital
      }else{
        option_tree[nrow(tree),k]<-0
      }
    }
    
  } else if(type=='put digital'){
    for(k in 1:nrow(tree)){
      if(tree[nrow(tree),k]>valordigital){
        option_tree[nrow(tree),k]<- valordigital
      }else{
        option_tree[nrow(tree),k]<-0
      }
    } 
    
  } else { 
    print("Ingresa un tipo valido de derivado")
  }
  for (i in (nrow(tree)-1):1) {
    for(j in 1:i) {
      option_tree[i,j]=((1-q)*option_tree[i+1,j] + q*option_tree[i+1,j+1])/exp(r*delta_t)
    }
  }
  return(option_tree)
}

'Arbol alpha'
subyacente <-function(arbolito,payoff){
  subyacente_tree = matrix(0, nrow=nrow(arbolito)-1, ncol=ncol(arbolito)-1) 
  for (i in (nrow(arbolito)-1):1) {
    for(j in 1:i) {
      subyacente_tree[i,j]=(payoff[i+1,j+1]-payoff[i+1,j])/(arbolito[i+1,j+1]-arbolito[i+1,j])
    }
  }
  return(subyacente_tree)
}

'Arbol beta'
cuentademercado <-function(arbolito,payoff,r,N,TTT){
  delta_t=TTT/N
  cuentademercado_tree = matrix(0, nrow=nrow(arbolito)-1, ncol=ncol(arbolito)-1) 
  for (i in (nrow(arbolito)-1):1) {
    for(j in 1:i) {
      cuentademercado_tree[i,j]=(payoff[i+1,j+1]-((payoff[i+1,j+1]-payoff[i+1,j])/(arbolito[i+1,j+1]-arbolito[i+1,j]))*arbolito[i+1,j+1])/exp(r*delta_t)
    }
  }
  return(cuentademercado_tree)
}

'Código para graficar el arbol'
my_BinomialTreePlot<-function (BinomialTreeValues,xlab,ylab,dx = -0.00025, dy = 0.004, cex = 1, 
                               digits = 2, ...) 
{
  #Redondeo de los valores
  Tree = round(BinomialTreeValues, digits = digits)
  #Tree=BinomialTreeValues
  depth = ncol(Tree)
  plot(x = c(0, depth-1), y = c(-depth + 1, depth - 1),axes=FALSE,xlab=xlab,ylab=ylab,type = "n", 
       col = 0, ...)
  #Coloque el primer valor
  text(0 + dx, 0 + dy, deparse(Tree[1, 1]), cex = cex)
  for (i in 1:(depth - 1)) {
    y = seq(from = -i, by = 2, length = i + 1)
    x = rep(i, times = length(y)) + 0
    for (j in 1:length(x)) text(x[j] + dx, y[j] + dy, deparse(Tree[length(x) + 
                                                                     1 - j, i + 1]), cex = cex)
    #Dibujar las líneas
    y = (-i):i
    x = rep(c(i, i-1), times = 2 * i)[1:length(y)]
    lines(x, y, col = 2)
  }
  invisible()
}

#'Imprimir'
arboldelprecio <- function(x) {
  arboldelpreciosubyacente<-my_BinomialTreePlot(t(mover(x)),ylab="Precio del subyacente",xlab="Paso")
  
  return(arboldelpreciosubyacente)
}

arboldelder <- function(y) {
  arboldelder <-my_BinomialTreePlot(t(mover(y)),ylab="Precio del derivado",xlab="Paso")
  
  return(arboldelder)
}

arbolalpha <- function(z) {
  arbolalpha <- my_BinomialTreePlot(t(mover(z)),ylab="Cantidad del bien subyacente",xlab="Paso")
  
  return(arbolalpha)
}

arbolbeta <- function(v) {
  arbolbeta <- my_BinomialTreePlot(t(mover(v)),ylab="Cantidad en la cuenta de mercado de dinero",xlab="Paso")  
  return(arbolbeta)
}

'Función que devuelve los precios y los arboles'
valuacionderivado <- function(type,S1,r,u,d,N,TTT,K,valordigital=0) {
  q <- q_prob(r,u, d,N,TTT)
  x <- construyearbolpreciossub(S1,u,d,N)
  y <- valuarderivado(x,S1,r,u,d,N,TTT, K, type,valordigital)
  z <- subyacente(x,y)
  v <-cuentademercado(x,y,r,N,TTT)
  return(list(q=q, arboldelprecio(x), arboldelder(y),arbolalpha(z),arbolbeta(v), prima=y[1,1]))
}

#### S1 es el precio del subyacente a tiempo 0
#### r es la tasa libre de riesgo, tiene que ponerse en decimales
#### u factor por el que sube el precio
#### d factor por el que baja el precio
#### k es el precio strike
#### r es la tasa libre de riesgo, tiene que ponerse en decimales
#### N es el número de pasos del árbol
#### TTT es el vencimiento de la opción en años
#### valor digital (caso de opciones digitales)
#### tipos: 'put europea plain vanilla'
#### 'call europea plain vanilla'
#### 'forward'
#### 'call digital'
#### 'put digital'

valuacionderivado(type='call europea plain vanilla',S1=100,r=0.045,u=1.2,d=.8,N=5,TTT=2,K=100,valordigital = 0)
#valuacionderivado(type='call digital',S1=100,r=0.045,u=1.2,d=.8,N=5,TTT=2,K=100,valordigital=100)
#valuacionderivado(type='forward',S1=100,r=0.045,u=1.2,d=.8,N=5,TTT=2,K=100,valordigital=0)
