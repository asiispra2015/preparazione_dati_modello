#aggiornamento 8 ottobre 2018
#Input: 
# file con dati pm10 estratti dal database postgis con il programma python di interrogazione al database
# Output:
# - variabili standardizzate
# - le distanze invertite e standardizzate
# - il parametro data_record_value rinominato in pm10
# - il parametro data_record rinominato yymmdd
# - aggiunta la variabile season rispetto al file di input
# Il file di output si chiama: pm10_analisi.csv

#ATTTENZIONE: le variabili temporali sono standardizzate stazione per stazione, le variabili
# spaziali sono invece standardizzate spazialmente. Per queste variabili sarebbe più corretto
#standardizzare in base ad area climatica? Per il momento questo problema non è affrontato
rm(list=objects())
library("readr")
library("dplyr")
library("purrr")
library("lubridate")
library("imputeTS") #per riempire gli NA dell'ndvi
options(warn=2,error=recover)

#file estratto mediante python da postgis contiene tutte le variabili associate al centroide
#rappresentativo delle centraline
nomeFile<-"pm10_25luglio2018.csv"

if(!file.exists(nomeFile)) stop(sprintf("File %s di input non trovato",nomeFile))



# Funzione per standardizzazione delle variabili quantitative --------------------------


#closure utilizzata per restituire scala o scala.inv quando si standardizzano le variabili spazio-temporali

scala<-function(x){if(all(is.na(x))) return(x); scale(x)->xx; imputeTS::na.interpolation(xx) }
scala.inv<-function(x){if(all(is.na(x))) return(x); scala(1/(x+1))->xx; imputeTS::na.interpolation(xx)}

creaFun<-function(nome){
  
  if(grepl("^pbl..$",nome)){
    print("PLANET BOUNDARY LAYER: scalo 1/x invece di x")
    scala.inv
  }else{
    print("scalo x")
    scala
  }
  
}#fine creaFun


# INIZIO PROGRAMMA --------------------------------------------------------

# Lettura file dati come estratto dal database postgis --------------------
read_delim(nomeFile,delim=";",col_names=TRUE) %>% mutate(data_record=as.character(data_record))->dati


###########################################################################################
SOGLIA<-1
###########################################################################################

###########################################################################################
#arrotondamento commerciale
###########################################################################################
arrotonda=function(x)
{
  
  x %% 1 -> resto
  ceiling(x[resto>=0.5 & !is.na(resto)])->x[resto>=0.5 & !is.na(resto)]
  floor(x[resto<0.5 & !is.na(resto)])->x[resto<0.5 & !is.na(resto)]
  
  x[is.nan(x)]<-NA
  
  x
  
}  
###########################################################################################
# fine arrotondamento commerciale
###########################################################################################



#rinomino data_record_value in pm10
if(!any(grepl("^pm10$",names(dati)))) names(dati)[grepl("^data_record_value$",names(dati))]<-"pm10"
#rinomino data_record in yymmdd
if(!any(grepl("^yymmdd$",names(dati)))) names(dati)[grepl("^data_record$",names(dati))]<-"yymmdd"

#creo la variabile rpm10, i valori di pm10 con arrotondamento commerciale.
#I valori al di sotto di SOGLIA li pongo uguali a SOGLIA.
dati %>% 
  mutate(rpm10=arrotonda(pm10)) %>%
  mutate(rpm10=ifelse(rpm10 < SOGLIA  & !is.na(rpm10),SOGLIA,rpm10))->dati




###########################################################################################
#Per utilizzare crossed random effects in nlme::lme dobbiamo introdurre due variabili dummy che servono per
#la creazione delle matrici identità con pdIdent
#crossed random effects con nlme:unconditional model
###########################################################################################
dati<-within(dati,one1<-one2<-1L)



####################################
#salviamo i valori di pm10 originali e creiamo una variabile di pm10 partendo da rpm10 (valore di pm10 arrotondato) + una componente random
#Quindi: pm10.orig contiene i dati di pm10 originali (alcuni valori interi, altri con cifre decimali); rpm10 contiene i valori di pm10 arrotondati; pm10 contiene i valori
#di rpm10 più una componente random (in modo di avere tutti i valori con cifre decimali)
####################################
dati$pm10->dati$pm10.orig
#aggiungiamo un valore random ai dati interi di rpm10
dati$pm10<-round(runif(nrow(dati),min=-0.5,max=0.5)+dati$rpm10,2)
dati[dati$pm10<SOGLIA & !is.na(dati$pm10),]$pm10<-SOGLIA

####################################
#distinguiamo tra i valori di pm10 >=THRESHOLD e quelli minori di THRESHOLD: serve solo a scopo visivio al momento della validazione
####################################
THRESHOLD<-5
dati$validi<-1
dati[dati$pm10< THRESHOLD & !is.na(dati$pm10), ]$validi<-0


####################################
# Creazione della stagione ------------------------------------------------
####################################
lubridate::month(as.Date(dati$yymmdd,"%Y-%m-%d"))->mesi
mesi->season

season[season %in% c(1,2,12)]<-1
season[season %in% c(3,4,5)]<-2
season[season %in% c(6,7,8)]<-3
season[season %in% c(9,10,11)]<-4
dati$season<-season
rm(mesi)
rm(season)


####################################
# Creazione del logaritmo dell'aod550, pbl00 e pbl12
####################################
dati$log.aod550<-log(dati$aod550)
dati$log.pbl00<-log(dati$pbl00)
dati$log.pbl12<-log(dati$pbl12)


########################### 
#Elenco delle variabili spaziali
########################### 
emissioni<-c("pm10_diff","nh3_diff","co_punt")
distanze<-c("d_a1","d_a2","d_costa","d_aero","d_impianti")
clc<-c("cl_agri","cl_arbl","cl_crop","cl_dcds","cl_evgr","cl_hidv","cl_lwdv","cl_pstr","cl_shrb")
osm<-c("av_buf_a1","av_buf_a23","av_buf_oth")
altre<-c("i_surface","q_dem","p_istat")

########################### 
#Elenco delle variabili spazio-temporali
########################### 
pbl<-c("pbl00","pbl12","log.pbl00","log.pbl12")
era5<-c("t2m","sp","u10","v10","tp")
ndvi<-c("ndvi")
aod<-c("aod550","log.aod550")

########################### 
#Elenco delle variabili qualitative
########################### 
varQualitative<-c("dust","climate_zone","season","cod_reg")

########################### 
#Variabili restanti nel file dati
########################### 
varAna<-c("yymmdd","id_centralina","banda","idcell","x","y")

####################################
#scalo variabili spazio-temporali, scalo le variabili per stazione
####################################
c(era5,ndvi,pbl)->pred1

purrr::map(pred1,.f=function(nomeVar){

  #myfun è la funzione per scalare la variabile x, sarà scala o scala.inv
  creaFun(nomeVar)->myfun
  
  dati[,c("yymmdd","id_centralina",nomeVar)]->mydf
  names(mydf)[3]<-"temp"

  mydf %>%
    tidyr::spread(key=id_centralina,value=temp)->ris0


    ris0 %>%
      dplyr::select(-yymmdd) %>%
        purrr::map_df(.,.f=myfun)->zz
      
    #attenzione id_centralina contiene ":" ma spread li elimina, vanno reinseriti pena il fallimento del join
    #anche per anno mese e giorno il passaggio dal wide format al long format cambia i caratteri
    data.frame(list(ris0[,c("yymmdd")],zz)) %>%
              tidyr::gather(id_centralina,temp,-yymmdd) %>%
                mutate(id_centralina=stringr::str_replace_all(id_centralina,"([0-9]{2})\\.([0-9]{2})\\.([0-9]{2})$","\\1:\\2:\\3")) %>%      
                  mutate(id_centralina=stringr::str_replace_all(id_centralina,"\\.","-"))->ris
    names(ris)[3]<-nomeVar

  ris
  
}) %>%
  purrr::reduce(left_join,by=c("yymmdd"="yymmdd","id_centralina"="id_centralina"))->daUnire0

names(daUnire0)<-c("yymmdd","id_centralina",paste0(era5,".s"),paste0(ndvi,".s"),paste0(pbl,".inv.s"))


#variabili da standardizzare ma non per stazione
purrr::map_df(dati[,c(emissioni,clc,osm,altre)],.f=scala)->pred2.s
names(pred2.s)<-paste0(c(emissioni,clc,osm,altre),".s")

#variabili da standardizzare invertire ma non per stazione
purrr::map_df(dati[,c(distanze)],.f=scala.inv)->pred3.s
names(pred3.s)<-paste0(c(distanze),".inv.s")

data.frame(list(dati[,c(varAna,varQualitative,aod)],pred2.s,pred3.s))->daUnire1
rm(list=objects(pattern="^pred[23]$"))
#rm(dati)


#finale data.frame con tutte le variabili
inner_join(daUnire0,daUnire1,by=c("id_centralina"="id_centralina","yymmdd"="yymmdd"))->finale

stopifnot(nrow(finale)==nrow(daUnire0))
stopifnot(nrow(finale)==nrow(daUnire1))

rm(daUnire0)
rm(daUnire1)

dplyr::left_join(dati %>% select(-climate_zone,-cod_reg,-dust,-season,-x,-y,-aod550,-log.aod550),finale,by=c("idcell"="idcell","banda"="banda","id_centralina"="id_centralina","yymmdd"="yymmdd")) %>%
  dplyr::select(id_centralina,yymmdd,banda,season,idcell,x,y,climate_zone,cod_reg,geom,pm10,rpm10,pm10.orig,everything())->daScrivere


#attenzione: ndvi contiene NA mentre ndvi.s (dove possibile) gli NA sono stati riempiti mediante na.interpolate.
#Quindi: per filtrare le righe con NA in ndvi utilizzare "ndvi.s" e non "ndvi" (altrimenti togliamo anche le righe dove è stato possibile fare l'imputing).
daScrivere[,grepl("ndvi.s",names(daScrivere))]->tmp
daScrivere[!is.na(tmp$ndvi.s),]->daScrivere
rm(tmp)
rm(dati)
rm(finale)
#l'ndvi (ndvi.s) contiene tutti NA in tre stazioni (in Veneto, Calabria e Sardegna)
#Per queste tre stazioni non è stato possibile fare l'imputing dei dati mancanti 
# IT0448A_5_BETA_2003-02-22_00:00:00 
# IT1940A_5_BETA_2009-03-01_00:00:00 
# IT2040A_5_BETA_2011-05-19_00:00:00
#

#scrittura file di output
write_delim(daScrivere,path="pm10_analisi.csv",delim=";",col_names = TRUE)
