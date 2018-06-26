
##--- Fonction de RA total et fonction de renegociation

fct_rat <-function(age=0,cf=NULL){
  #-----age en année 
  #-----Evaluation proba Remb Anticip structurel (sans RA partiel)
  #-----calibrage Avril 2017----
  
  #on reduit l'échelle des variables sur [0,1] 
  #- après les avoir bornées 
  age.max<-24 #la fonction prendra la valeur 0 en age.max
  age.int1<-4/age.max
  age.int2<-12.5/age.max
  dt.tmp<-data.table(x=age/age.max)
  dt.tmp[x>1,x:=1]
  dt.tmp[x<0,x:=0]
  dt.tmp[x<=age.int1,pb.rat:=0.07*x/age.int1]
  dt.tmp[x>age.int1 & x<=age.int2,pb.rat:=0.07]
  dt.tmp[x>age.int2,pb.rat:=-0.07*(x-age.int2)/(1-age.int2)+0.07]
  dt.tmp[pb.rat<0,pb.rat:=0]
  return(dt.tmp[,.(pb.rat)])
}

poly4.2var<-function(x,y,val.max=1,cf.x=c(0,0,0),cf.y=c(0,0),cf.var=c(0,0,0)){
  #----Fonction polynomiale evaluant la probabilité de renégociation
  #---- calibrage Avril 2017
  #---- appelée dans fct.ren.dg5
  res<- x*y*( cf.var[1]*(1 - x)*(1 - y^2) +
                cf.var[2]*(1 - x)*(1 - y) +
                cf.var[3]*(1 - y)*(1 - x^2) +
                val.max -
                cf.x[1]*(1 - x) -
                cf.x[2]*(1 - x^2) -
                cf.x[3]*(1 - x^3)-
                cf.y[1]*(1 - y) -
                cf.y[2]*(1 - y^2) )
  return(res)
}  

fct_renego<-function(matres=0,spread=0,slope=NULL,lst.cf.ren=NULL){ 
  #----Evaluation proba de renégocioation
  #----calibrage Avril 2017----
  #matres en année, spread en taux: 0.01=1%
  #- lst.cf.ren contient la liste des coefficients à appliquer dans la fonction poly4.2var 
  xmax=21
  ymax=0.0375
  
  #on reduit l'échelle des variables sur [0,1] 
  #- après les avoir bornées 
  dt.tmp<-data.table(x=matres/xmax,y=spread/ymax)
  
  dt.tmp[x>1,x:=1]
  dt.tmp[y<0,y:=0]
  dt.tmp[y>1,y:=1]
  dt.tmp[,pb.ren:= poly4.2var(x,y,lst.cf.ren[[4]],lst.cf.ren[[1]],lst.cf.ren[[2]],lst.cf.ren[[3]])]
  dt.tmp[pb.ren<0,pb.ren:=0]
  return(dt.tmp[,.(1-(1-pb.ren)^12)])
  #- on retourne le taux annualisé
}  #- Evaluation proba renégociation
