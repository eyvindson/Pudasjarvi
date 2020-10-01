# Calculates the the annual carbon balance of drained organic soils as is (approximately)
# done in the Greenhouse Gas Inventory, that is,
# balance = emission - below-ground litter input
#
# Emission is calculated with the aid of site specific emission factors and areas of site types 
# Litter input consists litter from living trees, harvest residues, natural mortality and
# understorey vegetation
#
# Litter from living trees, harvests and natural mortality is given separately for Southern and Northern
# Finland.

#Peatland shares are now (23.10.2018) here BUT will be moved to MELAToC_definitions.R

sh <- read.table("../ParameterData/Peatland-type-shares24.1.2019.txt",header=TRUE)

shares <- function(year) {
    if(year < 1990) {
        ret <- c(sh[sh$year == 1990,2:6][[1]],sh[sh$year == 1990,2:6][[2]],
                 sh[sh$year == 1990,2:6][[3]],sh[sh$year == 1990,2:6][[4]],
                 sh[sh$year == 1990,2:6][[5]])
    }
    else { if(year > 2016) {
        ret <- c(sh[sh$year == 2016,2:6][[1]],sh[sh$year == 2016,2:6][[2]],
               sh[sh$year == 2016,2:6][[3]],sh[sh$year == 2016,2:6][[4]],
               sh[sh$year == 2016,2:6][[5]])
       
    } else {
      ret <- c(sh[sh$year == year,2:6][[1]],sh[sh$year == year,2:6][[2]],
             sh[sh$year == year,2:6][[3]],sh[sh$year == year,2:6][[4]],
             sh[sh$year == year,2:6][[5]])
        
    } }
   ret
}



## LITTER INPUT

## UNDERSTOREY

#on understoreyt und.nwl.NF.min.csv  und.nwl.SF.min.csv  und.org.abv.csv  und.org.bel.csv
#und.org.abv <-  (read.csv(paste(undpath,"und.org.abv.csv", sep=""), header=TRUE))
#undBel <-  (read.csv(paste(undpath,"und.org.bel.csv", sep=""), header=TRUE))

undBel <-  (read.csv(UnderstoreyOrganicBelow, header=TRUE))

## LITTER FROM TREES

treelitter.sf<-  read.table(OrganicSoilBelowSouth, header=TRUE)
treelitter.nf<-  read.table(OrganicSoilBelowNorth, header=TRUE)


## LITTER FROM LOGGING RESIDUES

#melaWaste<- read.table(MelaWaste, header=TRUE)
melaWaste <- MelaWaste


melaTimes <- MelaYears
melaN <- length(melaTimes)
melaNlitter <- melaN - 1 

fol.rat.sf <-  FoliageRatioSouth
fol.rat.nf <-  FoliageRatioNorth

#wastelitter for output
wastelitt<-matrix(data=0,nrow=melaNlitter,ncol=2)

sumlitterWasteBel.sf <-matrix(data=0,nrow=melaNlitter,ncol=1)
sumlitterWasteBel.nf <-matrix(data=0,nrow=melaNlitter,ncol=1)

for(iregion in 1:2) {

  sumlitter <-matrix(data=0,nrow=melaNlitter,ncol=1)
    for(ispecies in 1:3) {
      for(ipart in 1:5) {
        for(itimes in 1:melaNlitter) {
           if(ipart == 1) { #only fine roots here
             if(ispecies == 1) {
               if(iregion==1) flr <- fol.rat.sf[1] else flr <- fol.rat.nf[1]}
               else {
               if(ispecies == 2) {
               if(iregion==1) flr <- fol.rat.sf[2] else flr <- fol.rat.nf[2]}
               else{
               flr = 1}
             }
             sumlitter[itimes] = sumlitter[itimes] +
             carbon(fineroots(flr*1000.0*melaWaste[melaWaste$region==iregion&
             melaWaste$soil==2&melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes],
             ispecies))
             wastelitt[itimes,1] = wastelitt[itimes,1] +
             carbon(fineroots(flr*1000.0*melaWaste[melaWaste$region==iregion&
             melaWaste$soil==2&melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes],
             ispecies))
           }
            if(ipart == 5) {
              sumlitter[itimes] = sumlitter[itimes] +
              carbon(1000.0*melaWaste[melaWaste$region==iregion&melaWaste$soil==2
              &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes])
              wastelitt[itimes,1] = wastelitt[itimes,1] +
              carbon(1000.0*melaWaste[melaWaste$region==iregion&melaWaste$soil==2
              &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes])                        
            }
            if(ipart == 2 | ipart == 3 | ipart == 4 ) {
            wastelitt[itimes,2] = wastelitt[itimes,2] +
            carbon(1000.0*melaWaste[melaWaste$region==iregion&melaWaste$soil==2
            &melaWaste$species==ispecies&melaWaste$part==2,4+itimes]) +
            carbon(1000.0*melaWaste[melaWaste$region==iregion&melaWaste$soil==2
            &melaWaste$species==ispecies&melaWaste$part==3,4+itimes]) +
            carbon(1000.0*melaWaste[melaWaste$region==iregion&melaWaste$soil==2
            &melaWaste$species==ispecies&melaWaste$part==4,4+itimes])
            }
           
        }
      } #for(ipart in
    } #for(ispecies in
         
 if(iregion == 1) {
    sumlitterWasteBel.sf <- sumlitter
  } else {
  sumlitterWasteBel.nf <- sumlitter
  }
}


## LITTER FROM NATURAL MORTALITY

#melaNatMort = read.table(MelaNatMortality, header=TRUE)
melaNatMort = MelaNatMortality
## Natural mortality is specified in terms of stand volume in melaNatMort

sumlitterNMortBel.sf <-matrix(data=0,nrow=melaNlitter,ncol=1)
sumlitterNMortBel.nf <-matrix(data=0,nrow=melaNlitter,ncol=1)

#mortlitter for output
mortlitt<-matrix(data=0,nrow=melaNlitter,ncol=2)


for(iregion in 1:2) {
    sumlitter <-matrix(data=0,nrow=melaNlitter,ncol=1)
    
  for(ispecies in 1:3) {
     for(itimes in 1:melaNlitter) {
         stemVol <- 1000.0*melaNatMort[melaNatMort$region==iregion&melaNatMort$soil==2
                    &melaNatMort$species==ispecies,3+itimes]
         for(ipart in 1:5) {
             if(ipart == 1)  {#only fine roots from foliage
                if(ispecies == 1) {
                   flr <- fol.rat.sf[1]} else {
                   if(ispecies == 2) {
                     flr = fol.rat.sf[1]} else{flr = 1}
                 }
                 sumlitter[itimes] = sumlitter[itimes] +
                 carbon(fineroots(flr*functionBEF(iregion,2,ispecies,ipart)*stemVol,
                                 ispecies))
                 mortlitt[itimes,1] = mortlitt[itimes,1] +
                 carbon(fineroots(flr*functionBEF(iregion,2,ispecies,ipart)*stemVol,
                                 ispecies))
             } #if(ipart ==1
             if(ipart == 5) {  #Roots are bel
                 sumlitter[itimes] = sumlitter[itimes] +
                 carbon(functionBEF(iregion,2,ispecies,ipart)*stemVol)
                 mortlitt[itimes,1] = mortlitt[itimes,1] +
                 carbon(functionBEF(iregion,2,ispecies,ipart)*stemVol)
             }
             if(ipart == 2 | ipart == 3 | ipart == 4) { #all rest is above-ground
                 mortlitt[itimes,2] = mortlitt[itimes,2] +
                 carbon(functionBEF(iregion,2,ispecies,2)*stemVol) +
                 carbon(functionBEF(iregion,2,ispecies,3)*stemVol) +
                 carbon(functionBEF(iregion,2,ispecies,4)*stemVol)
                         }
          }#for(ipart
       }
   }
  if(iregion == 1) {
     sumlitterNMortBel.sf <- sumlitter
  } else {
     sumlitterNMortBel.nf <- sumlitter
  }
}


## CALCULATION OF THE BALANCE

melaYears = seq(from=melaTimes[1], to=melaTimes[length(melaTimes)], by = 1)
approxTimes <- as.integer(melaTimes[1:(melaN-1)] + diff(melaTimes)/2)

#Areas of drained organic soils in Southern and Northern Finland
MelaADraindedOrg.sf = Areas$area[Areas$region==1&Areas$soil==2]
MelaADraindedOrg.nf = Areas$area[Areas$region==2&Areas$soil==2]

# Litter
litterT.sf <- MelaADraindedOrg.sf*rowSums(treelitter.sf[,2:9])[1:length(melaYears)]
litterT.nf <- MelaADraindedOrg.nf*rowSums(treelitter.nf[,2:9])[1:length(melaYears)]
litterW.sf <- approx(approxTimes,sumlitterWasteBel.sf,melaYears,rule=2)$y
litterW.nf <- approx(approxTimes,sumlitterWasteBel.nf,melaYears,rule=2)$y
litterM.sf <- approx(approxTimes,sumlitterNMortBel.sf,melaYears,rule=2)$y
litterM.nf <- approx(approxTimes,sumlitterNMortBel.nf,melaYears,rule=2)$y

#LitterSouth and litterNorth depend on input from trees
litterSouth <- litterT.sf + litterW.sf + litterM.sf
litterNorth <- litterT.nf + litterW.nf + litterM.nf
if(!calculationRegions[1]) {
    litterT.sf <- 0
    litterW.sf <- 0
    litterM.sf <- 0
}
if(!calculationRegions[2]) {
    litterT.nf <- 0
    litterW.nf <- 0
    litterM.nf <- 0
}


litterT <- litterT.sf + litterT.nf
litterW <- litterW.sf + litterW.nf
litterM <- litterM.sf + litterM.nf


#Litters for output
wasteLitter <- data.frame(cbind(melaYears,
               approx(approxTimes,wastelitt[,1],melaYears,rule=2)$y,
               approx(approxTimes,wastelitt[,2],melaYears,rule=2)$y))
mortLitter <- data.frame(cbind(melaYears,
               approx(approxTimes,mortlitt[,1],melaYears,rule=2)$y,
               approx(approxTimes,mortlitt[,2],melaYears,rule=2)$y))
    

#### Annual emissions from drained peat sites  g C m-2 a-1 to ton C /ha (note CARBON)
Rhtkg <- 425.7*0.01
Mtkg <- 312.1*0.01
Ptkg <- 242.3*0.01
Vatkg <- 218.9*0.01
Jatkg <- 185.2*0.01

#Areas of types of peatland  year 2008 (NIR 2010) 1000 ha
# Rhtkg        Mtkg       Ptkg      Vatkg       Jatkg
# 688         1150       1628        860         41
#relMireTypeAreas = c(688, 1150, 1628, 860, 41) / sum(688, 1150, 1628, 860, 41)

#Both emission and understorey litter input per unit area are constant over
#the whole county. If no litter input is coming from trees either south or north,
# emission and understorey litter input are equal to zero in the respective area.
emissionSouth <- 1:length(melaYears)
if(sum(litterSouth) > 1e-6) {
    for(i in 1:length(melaYears)) {   
#    emissionSouth <- MelaADraindedOrg.sf * (relMireTypeAreas[1]*Rhtkg +
#                  relMireTypeAreas[2]*Mtkg + relMireTypeAreas[3]*Ptkg +
#                  relMireTypeAreas[4]*Vatkg + relMireTypeAreas[5]*Jatkg)
    rMTA <- shares(melaYears[i])
    emissionSouth[i] <- MelaADraindedOrg.sf * (rMTA[1]*Rhtkg +
                      rMTA[2]*Mtkg + rMTA[3]*Ptkg +
                      rMTA[4]*Vatkg + rMTA[5]*Jatkg)
    }
    litterUSouth <- MelaADraindedOrg.sf * sum(undBel)
} else {
    emissionSouth  <- 0
    litterUSouth <- 0
}
emissionNorth <- 1:length(melaYears)
if(sum(litterNorth) > 1e-6) {
    for(i in 1:length(melaYears)) {        
#    emissionNorth <- MelaADraindedOrg.nf * (relMireTypeAreas[1]*Rhtkg +
#                      relMireTypeAreas[2]*Mtkg + relMireTypeAreas[3]*Ptkg +
#                      relMireTypeAreas[4]*Vatkg + relMireTypeAreas[5]*Jatkg)
   rMTA <- shares(melaYears[i])
   emissionNorth[i] <- MelaADraindedOrg.nf * (rMTA[1]*Rhtkg +
                    rMTA[2]*Mtkg + rMTA[3]*Ptkg +
                    rMTA[4]*Vatkg + rMTA[5]*Jatkg)
    }
    litterUNorth <- MelaADraindedOrg.nf * sum(undBel)
    
} else {
    emissionNorth  <- 0
    litterUNorth <- 0
}

#emission = emissionSouth + emissionNorth
if(!calculationRegions[1]) {
        litterUSouth <- 0
}
if(!calculationRegions[2]) {
         litterUNorth <- 0
}

litterU <- litterUSouth + litterUNorth


litter <- litterU + litterT + litterW + litterM

##########################

talitt <- as.data.frame(cbind(melaYears,litterT.sf/MelaADraindedOrg.sf,litterW.sf/MelaADraindedOrg.sf,
     litterM.sf/MelaADraindedOrg.sf,litterUSouth/MelaADraindedOrg.sf,litterT.nf/MelaADraindedOrg.nf,
     litterW.nf/MelaADraindedOrg.nf,litterM.nf/MelaADraindedOrg.nf,litterUNorth/MelaADraindedOrg.nf))
colnames(talitt) <- c("year", "lTreesSo", "lWasteSo", "lMortSo", "lUndSo",
                      "lTreesNo", "lWasteNo", "lMortNo", "lUndNo")


######################################


if(EarlierStart) {
    ghgLitt <- read.table(GHGILitterOrganicFile,header=TRUE)
    cond <- ghgLitt$year >= SimulationStartYear & ghgLitt$year < melaTimes[1]
    ghg <- MelaADraindedOrg.sf*ghgLitt$South[cond] +
        MelaADraindedOrg.nf*ghgLitt$North[cond]
    ghgEmissionSouth <-1:sum(cond)          #for those years for which cond = TRUE
    ghgEmissionNorth <-1:sum(cond)          #for those years for which cond = TRUE
    j = 0
    for(i in 1:length(ghgLitt$year)) {
        if(cond[i]) {
            j <- j + 1
           rMTA <- shares(ghgLitt$year[i])
           ghgEmissionSouth[j] <- MelaADraindedOrg.sf * (rMTA[1]*Rhtkg +
           rMTA[2]*Mtkg + rMTA[3]*Ptkg +
           rMTA[4]*Vatkg + rMTA[5]*Jatkg)
           ghgEmissionNorth[j] <- MelaADraindedOrg.nf * (rMTA[1]*Rhtkg +
           rMTA[2]*Mtkg + rMTA[3]*Ptkg +
           rMTA[4]*Vatkg + rMTA[5]*Jatkg)
        }
    }
# GHGI litter are for years 1990 - 2016. If start is earlier than 1990
# ghg needs to augmented with zeros before 1990
    if(SimulationStartYear < 1990) {
	nothing <- rep(0,times = (1990 - SimulationStartYear))
	ghg <- c(nothing,ghg)
        ghgEmissionSouth <- c(nothing,ghgEmissionSouth)
        ghgEmissionNorth <- c(nothing,ghgEmissionNorth)
    }
    litter <- c(ghg,litter)
    emissionSouth <- c(ghgEmissionSouth,emissionSouth)
    emissionNorth <- c(ghgEmissionNorth,emissionNorth)
}

if(!calculationRegions[1]) {
emissionSouth <- 0
}
if(!calculationRegions[2]) {
    emissionNorth <- 0
}

emission = emissionSouth + emissionNorth
cBalance <- litter - emission

# Effect of (above-ground) dead wood pool - constant sink after 1997 or 2001
# This is from GHGI (Paula Ollila) 11.10.2018
#dom.sf <- c(0,0,0,0,0,0,0,0,rep(0.00754106277747484,(yr-1997)))
#dom.nf <- c(0,0,0,0,0,0,0,0,0,0,0,0, rep(0.00720002212248889,(yr-2001)))
#unit is t C / ha (/year)

if(!OnlyDOMSOM) {   #If OnlyDOMSOM deadwood excluded from calculation
        if(!EarlierStart) {
                    timeSpan <- melaYears
        } else {
                    timeSpan <- seq(SimulationStartYear, to=max(melaTimes), by = 1)
        }

            cond.sf <- timeSpan > 1997
            cond.nf <- timeSpan > 2001
            deadWoodChange.sf <- MelaADraindedOrg.sf * cond.sf * 0.00754106277747484
            deadWoodChange.nf <- MelaADraindedOrg.nf * cond.nf * 0.00720002212248889

            if(calculationRegions[1]) {
                        cBalance <- cBalance + deadWoodChange.sf
            }
            if(calculationRegions[2]) {
                        cBalance <- cBalance + deadWoodChange.nf
            }
}


## AN OPTION TO MAKE EMISSION FACTORS DEPENDENT ON STEM VOLUME

# CMON-mallin mukaan myös puuston määrä vaikuttaa emissiokertoimiin:
# 0,833 g CO2/m2/vuosi, eli puuston tilavuuden kasvu 1 m3/ha lisää päästöä tuon verran.
#voladd<-approxfun(MealYears,MelaMeanV - MelaMeanV[1],rule=2)

#cBalance <- cBalance - A*voladd(melaYears)*0.01*0.833*12/44


## POSSIBLY THE EFFECT CLIMATE CHANGE

#Artikkelin summamallissa, jossa PPA:n tilalla oli V, touko-lokakuun keskilämpötilan kerroin oli:
#135 g CO2/m2/vuosi, eli keskilämpötilan nousu 1 asteella lisää päästöä tuon verran.
#Koska kyse on lineaarisesta summamallista, homma hoituu tuolla kertoimella suoraan.
#Aleksi teki mallin PPA:lla, ehkä lämpötilan kerroin on siinä hiukan eri?
#t. Paavo

#HUOM LCF Ilmastonmuutos on vaihdettu vaikuttamaan vuodesta 2010 alkaen
#siten, etta vahennetaan funktioista v. 2010 tilanne  == siis oletetaan,
#etta kaikki kuvautuu oikein v. 2010 ja muuttuva ilmasto vaikuttaa siita eteenpain

#tyrs<-c(2000,2010,2020,2030,2040,2050,2060,2070)
#tadd.sf<-c(0.0,0.55,0.97,1.33,1.62,1.85, 2.01, 2.11)          #SFin
#tadd.sf<-c(0.0,0.56, 1.00, 1.38, 1.68, 1.92, 2.09, 2.19)     #CFin
#tadd.nf<-c(0.0,0.62, 1.10, 1.51, 1.85, 2.11, 2.30, 2.41)
#dt.sf<-approxfun(tyrs,tadd.sf,rule=2)
#dt.nf<-approxfun(tyrs,tadd.nf,rule=2)
#tassa LCF muutos
#dt.sf<-approxfun(tyrs,tadd.sf-dt.sf(2010),rule=2)
#dt.nf<-approxfun(tyrs,tadd.nf-dt.nf(2010),rule=2)
#CliPan muutos 2018
#dt.sf<-approxfun(tyrs,tadd.sf-dt.sf(2011),rule=2)
#dt.nf<-approxfun(tyrs,tadd.nf-dt.nf(2011),rule=2)
dt.sf <- DeltaTSouth
dt.nf <- DeltaTNorth


if(ClimateChange) cc = 1.0 else cc = 0.0

#OBS  * 0.7 * == annual mean temperature --> mean temperature may - oct 
#cBalance <- cBalance - cc * A * 135.0*0.01*(12/44)*0.7*(dt.sf(melaYears)+dt.nf(melaYears))/2.0
#cBalance <- cBalance - cc *(MelaADraindedOrg.sf+MelaADraindedOrg.nf) *
#            135.0*0.01*(12/44)*0.7*(dt.sf(melaYears)+dt.nf(melaYears))/2.0

if(!EarlierStart) {
    cBalance <- cBalance - cc * MelaADraindedOrg.sf *
                135.0*0.01*(12/44)*0.7 * dt.sf(melaYears)
    cBalance <- cBalance - cc * MelaADraindedOrg.nf *
                135.0*0.01*(12/44)*0.7 * dt.nf(melaYears)
} else {   #CC does not affect the balance prior MELA period (melaYears) 
    nolla <- rep(0,times=melaYears[1]-SimulationStartYear)
    cBalance <- cBalance - cc * MelaADraindedOrg.sf *
                 135.0*0.01*(12/44)*0.7 * c(nolla,dt.sf(melaYears))
    cBalance <- cBalance - cc * MelaADraindedOrg.nf *
                 135.0*0.01*(12/44)*0.7 * c(nolla,dt.nf(melaYears))
}


if(!EarlierStart) {
    organic <- cbind(melaYears,litter,cBalance,wasteLitter[,2:3],mortLitter[,2:3],litterT,litterU)
} else {
    nolla <- rep(0,times=melaYears[1]-SimulationStartYear)
    organic <- as.data.frame(cbind(c(SimulationStartYear:(melaYears[1]-1),melaYears),
                     litter,cBalance,       #these already augmented above with prior to MELA 
                     rbind(as.matrix(cbind(nolla,nolla)),as.matrix(wasteLitter[,2:3])),
                     rbind(as.matrix(cbind(nolla,nolla)),as.matrix(mortLitter[,2:3])),
                     c(nolla,litterT),c(nolla,rep(litterU,times=length(melaYears)))))
}
    
    colnames(organic)=c("year","orgLitter","C_organic","owaste_bel","owaste_abv","omort_bel","omort_abv","otree_bel","ound")

