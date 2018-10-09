#Input: file con dati pm10 estratti dal database postgis con il programma python di interrogazione al database
# Output:
# - variabili standardizzate
# - le distanze invertite e standardizzate
# - il parametro data_record_value rinominato in pm10
# - il parametro data_record rinominato yymmdd
# - aggiunta la variabile season rispetto al file di input

#Il file di output si chiama: pm10_analisi.csv

#aggiornamento 8 ottobre 2018

rm(list=objects())
library("readr")
library("dplyr")
library("purrr")
library("lme4")
library("lubridate")
library("lattice")
options(warn=2,error=recover)

nomeFile<-"pm10_25luglio2018.csv"
if(!file.exists(nomeFile)) stop(sprintf("File %s di input non trovato",nomeFile))
# Lettura file dati come estratto dal database postgis --------------------
read_delim(nomeFile,delim=";",col_names=TRUE)->dati

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


#imputing NDVI
#scala<-function(x){if(all(is.na(x))) return(x); scale(x)->xx; imputeTS::na.interpolation(xx) } #qui nn lo posso fare, altrimenti mi riempie una serie
#con gli estremi di due serie adiacenti. Devo fare imputing serie per serie
purrr::map(unique(dati$id_centralina),.f=function(idC){
  
  dati %>%
    filter(id_centralina==idC)->subDati
  stopifnot(nrow(subDati)<=365)
  
  if(all(is.na(subDati$ndvi))) return(subDati[,c("yymmdd","id_centralina","ndvi")])
  na.interpolation(subDati$ndvi)->new.ndvi
  
  data.frame(yymmdd=subDati$yymmdd,id_centralina=subDati$id_centralina,ndvi=new.ndvi)
  
  
}) %>% reduce(rbind) %>% as.data.frame->NDVI
NDVI$id_centralina<-as.character(NDVI$id_centralina)
left_join(dati %>% select(-ndvi) ,NDVI,by=c("yymmdd"="yymmdd","id_centralina"="id_centralina"))->dati
rm(NDVI)

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



########################### variabili quantitative (senza aod)
emissioni<-c("pm10_diff","nh3_diff","co_punt")
pbl<-c("pbl00","pbl12","log.pbl00","log.pbl12")
distanze<-c("d_a1","d_a2","d_costa","d_aero","d_impianti")
era5<-c("t2m","sp","u10","v10","tp")
clc<-c("cl_agri","cl_arbl","cl_crop","cl_dcds","cl_evgr","cl_hidv","cl_lwdv","cl_pstr","cl_shrb")
osm<-c("av_buf_a1","av_buf_a23","av_buf_oth")
altre<-c("ndvi","i_surface","q_dem","p_istat")

pred1<-c(era5,altre,clc,osm,emissioni)

#di queste variabili standardizziamo l'inverso 
pred2<-c(pbl,distanze)

#variabili restanti
pred4<-c("yymmdd","id_centralina","banda","idcell")

# Funzione per standardizzazione delle variabili quantitative --------------------------
scala<-function(x){scale(x)}




#scalo variabili in pred1
purrr::map_df(dati[,pred1],.f=scala)->pred1.s
names(pred1.s)<-paste0(names(pred1.s),".s")

#scalo variabili in pred2
purrr::map_df(dati[,pred2],.f=~(scala(1/(.+1))))->pred2.inv.s
names(pred2.inv.s)<-paste0(names(pred2.inv.s),".inv.s")

#finale data.frame con tutte le variabili
finale<-purrr::reduce(list(pred1.s,pred2.inv.s,dati[,pred4]),cbind)
rm(list=objects(pattern="^pred.+$"))
#rm(dati)

dplyr::left_join(dati,finale,by=c("idcell"="idcell","banda"="banda","id_centralina"="id_centralina","yymmdd"="yymmdd")) %>%
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
