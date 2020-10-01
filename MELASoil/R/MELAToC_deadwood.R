source(YassoFile)
dyn.load(YassoSharedLibrary)

## MELAToC_deadwood.R #########################################################
## This R program (source("MELAToC_deadwood.R")) calculates C stock change of
## deadwood by Yasso (= input consists of only large tree parts)
## Output in data frame deadwoodresult


################################################################################
# MELA data
# Waste contains biomass of all material left in harvesting (whether normal or
# energy) in the forest.
# Each line contains: region soil species part region data2007 ... data 2047
# region = 1,2 (southern, northern Finland)
# soil = 1,2,3 (upland, ditched and non-ditched peatland)
# species = 1,2,3 (pine, spruce, birch & other decidious)
# part = 1,2,3,4,5 (foliage, branches, stem, stump, roots (coarse))

#NatMort is as Waste but only stem volume for region soil species

# sf.min is litter from standing trees southern finland produced by Clean_volToCInput.R
# columns are "year" "nwl_A" "nwl_W" "nwl_E" "nwl_N" "fwl_A" "fwl_W" "fwl_E" "fwl_N"
#( _A is Acid soluble etc)
# rows are years 2007 ... 2057
# unit is t C / ha

#nf.min is as sf_min but for northern Finland

# Initial states for upland soils (southern Finland, northern Finland) are
# read from a file.
# Columns are: nwl_south" "fwl_south" "cwl_south" "nwl_north" "fwl_north" "cwl_north"
# Rows are A W E N H
###############################################################################

melaWaste <- MelaWaste
melaNatMort <- MelaNatMortality

#MELA areas
sf.kankaat = Areas$area[Areas$region==1&Areas$soil==1]
nf.kankaat = Areas$area[Areas$region==2&Areas$soil==1]
sf.turvemaat = Areas$area[Areas$region==1&Areas$soil==2]
nf.turvemaat = Areas$area[Areas$region==2&Areas$soil==2]

#Sum litter inputs melaWaste, melaNatMort

melaTimes <- MelaYears

melaN = length(melaTimes)
melaNlitter = melaN - 1


sumlitter.nwl.sf <- matrix(data=0,nrow=melaNlitter,ncol=4)
sumlitter.fwl.sf <- matrix(data=0,nrow=melaNlitter,ncol=4)
sumlitter.cwl.sf <- matrix(data=0,nrow=melaNlitter,ncol=4)
sumlitter.nwl.nf <- matrix(data=0,nrow=melaNlitter,ncol=4)
sumlitter.fwl.nf <- matrix(data=0,nrow=melaNlitter,ncol=4)
sumlitter.cwl.nf <- matrix(data=0,nrow=melaNlitter,ncol=4)

## First melaWaste, that is cutting residues

#wastelitter for output
wastelitt<-matrix(data=0,nrow=melaNlitter,ncol=6)


## Then litter from natural mortality - all of it is litter. Natural mortality
## is specified in terms of stand volume in melaNatMort

#mortlitter for output
mortlitt<-matrix(data=0,nrow=melaNlitter,ncol=6)
for(iregion in 1:2) {
    sumlitter.nwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
    sumlitter.fwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
    sumlitter.cwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
  for(ispecies in 1:3) {
        for(itimes in 1:melaNlitter) {
          stemVol <- melaNatMort[melaNatMort$region==iregion&melaNatMort$soil==DeadwoodSoil
                     &melaNatMort$species==ispecies,3+itimes]
          ipart <- 3
             sumlitter.cwl[itimes,] = sumlitter.cwl[itimes,] +
             stem.AWEN(carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol)
             sumlitter.cwl[itimes,] = sumlitter.cwl[itimes,] +
             stem.AWEN(carbon(functionBEF(iregion,1,ispecies,6)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,6)*stemVol)   #Bark
       }
  }
  if(iregion == 1) {
     sumlitter.nwl.sf <- sumlitter.nwl.sf + 1.0*sumlitter.nwl
     sumlitter.fwl.sf <- sumlitter.fwl.sf + 1.0*sumlitter.fwl
     sumlitter.cwl.sf <- sumlitter.cwl.sf + 1.0*sumlitter.cwl
    } else {
     sumlitter.nwl.nf <- sumlitter.nwl.nf + 1.0*sumlitter.nwl
     sumlitter.fwl.nf <- sumlitter.fwl.nf + 1.0*sumlitter.fwl
     sumlitter.cwl.nf <- sumlitter.cwl.nf + 1.0*sumlitter.cwl
    }
}


#This set again below, depending whther EarlierStart is active. Here for MELA litter
simulationPeriod <- melaTimes[1]:max(melaTimes)

#Approximation of wastelitter and mortlitter that are a step function = constant during 10 year
#periods: the value is at the midle of the period except the first one.

#approxTimes = value time of litter input. Best (must) be between times when MELA stock values
#have been given. The litter from harvest (waste) and natural mortality is constant for periods
# between MELA stock values - most natural
#is thus in the middle of the period. rule = 2 gives constant value outside approxTimes, thus
#litter input from harvest and natural mortality = constant in the first half of the first MELA period.
#Litter from standing trees (sf.min and nf.min) change from the start because it is calculated with
#the aid of stock values that change between melaTimes.

sf.min <- read.table(MineralSoilSouth, header=TRUE)
nf.min <- read.table(MineralSoilNorth, header=TRUE)


approxTimes = as.integer(melaTimes[1:(melaN-1)] + diff(melaTimes)/2)

inputLitter.sf <- data.frame(cbind(simulationPeriod,
    sf.min[,2] + 1000*approx(approxTimes,sumlitter.nwl.sf[,1],simulationPeriod,rule=2)$y,
    sf.min[,3] + 1000*approx(approxTimes,sumlitter.nwl.sf[,2],simulationPeriod,rule=2)$y,
    sf.min[,4] + 1000*approx(approxTimes,sumlitter.nwl.sf[,3],simulationPeriod,rule=2)$y,
    sf.min[,5] + 1000*approx(approxTimes,sumlitter.nwl.sf[,4],simulationPeriod,rule=2)$y,
    sf.min[,6] + 1000*approx(approxTimes,sumlitter.fwl.sf[,1],simulationPeriod,rule=2)$y,
    sf.min[,7] + 1000*approx(approxTimes,sumlitter.fwl.sf[,2],simulationPeriod,rule=2)$y,
    sf.min[,8] + 1000*approx(approxTimes,sumlitter.fwl.sf[,3],simulationPeriod,rule=2)$y,
    sf.min[,9] + 1000*approx(approxTimes,sumlitter.fwl.sf[,4],simulationPeriod,rule=2)$y,
    1000*approx(approxTimes,sumlitter.cwl.sf[,1],simulationPeriod,rule=2)$y,
    1000*approx(approxTimes,sumlitter.cwl.sf[,2],simulationPeriod,rule=2)$y,    #no cwl litter from
    1000*approx(approxTimes,sumlitter.cwl.sf[,3],simulationPeriod,rule=2)$y,    #standing trees
    1000*approx(approxTimes,sumlitter.cwl.sf[,4],simulationPeriod,rule=2)$y))

 inputLitter.nf <- data.frame(cbind(simulationPeriod,
    nf.min[,2] + 1000*approx(approxTimes,sumlitter.nwl.nf[,1],simulationPeriod,rule=2)$y,
    nf.min[,3] + 1000*approx(approxTimes,sumlitter.nwl.nf[,2],simulationPeriod,rule=2)$y,
    nf.min[,4] + 1000*approx(approxTimes,sumlitter.nwl.nf[,3],simulationPeriod,rule=2)$y,
    nf.min[,5] + 1000*approx(approxTimes,sumlitter.nwl.nf[,4],simulationPeriod,rule=2)$y,
    nf.min[,6] + 1000*approx(approxTimes,sumlitter.fwl.nf[,1],simulationPeriod,rule=2)$y,
    nf.min[,7] + 1000*approx(approxTimes,sumlitter.fwl.nf[,2],simulationPeriod,rule=2)$y,
    nf.min[,8] + 1000*approx(approxTimes,sumlitter.fwl.nf[,3],simulationPeriod,rule=2)$y,
    nf.min[,9] + 1000*approx(approxTimes,sumlitter.fwl.nf[,4],simulationPeriod,rule=2)$y,
    1000*approx(approxTimes,sumlitter.cwl.nf[,1],simulationPeriod,rule=2)$y,
    1000*approx(approxTimes,sumlitter.cwl.nf[,2],simulationPeriod,rule=2)$y,    #no cwl litter from
    1000*approx(approxTimes,sumlitter.cwl.nf[,3],simulationPeriod,rule=2)$y,    #standing trees
    1000*approx(approxTimes,sumlitter.cwl.nf[,4],simulationPeriod,rule=2)$y))


    colnames(inputLitter.sf) <- c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W",
                                  "fwl_E","fwl_N","cwl_A","cwl_W","cwl_E","cwl_N")
    colnames(inputLitter.nf) <- c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W",
                                  "fwl_E","fwl_N","cwl_A","cwl_W","cwl_E","cwl_N")

#litters for output
wasteLitter <- data.frame(cbind(simulationPeriod,
        1000.0*approx(approxTimes,wastelitt[,1],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,wastelitt[,2],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,wastelitt[,3],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,wastelitt[,4],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,wastelitt[,5],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,wastelitt[,6],simulationPeriod,rule=2)$y))

mortLitter <- data.frame(cbind(simulationPeriod,
        1000.0*approx(approxTimes,mortlitt[,1],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,mortlitt[,2],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,mortlitt[,3],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,mortlitt[,4],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,mortlitt[,5],simulationPeriod,rule=2)$y,
        1000.0*approx(approxTimes,mortlitt[,6],simulationPeriod,rule=2)$y))
colnames(wasteLitter) <- c("year","wafol","wabra","wastem","wastump","waroots","wafr")
colnames(mortLitter) <- c("year","mofol","mobra","mostem","mostump","moroots","mofr")


## Check if input litter from trees == 0 in south or north
## If treeinput == 0, Yasso inital value needs to be equal to 0.
if(sum(inputLitter.sf[,2:13]) > 1e-6) {
    treeInputCheckSouth <- 1
} else {
    treeInputCheckSouth <- 0
}
if(sum(inputLitter.nf[,2:13]) > 1e-6) {
    treeInputCheckNorth <- 1
} else {
    treeInputCheckNorth <- 0
}


## SIMULATION MAY START EARLIER THAN AT FIRST YEAR OF melaTimes (=MelaYears), LITTER INPUT
## IS TAKEN THEN FROM LITTER USED IN GREENHOUSE GAS INVENTORY IN MINERAL SOIL ((YASSO)

if(EarlierStart) {
	simulationPeriod <- SimulationStartYear:max(melaTimes)
	
	ghgLitt <- read.table(GHGILitterMineralFile,header=TRUE)
        cond <- ghgLitt$region == 1 & ghgLitt$year >= SimulationStartYear & ghgLitt$year < melaTimes[1]
        ghg <- cbind(ghgLitt[cond,1],sf.kankaat*ghgLitt[cond,3:14])
        colnames(ghg)[1]="year"
	inputLitter.sf <- rbind(ghg,inputLitter.sf)  #South

        cond <- ghgLitt$region == 2 & ghgLitt$year >= SimulationStartYear & ghgLitt$year < melaTimes[1]
         ghg <- cbind(ghgLitt[cond,1],nf.kankaat*ghgLitt[cond,3:14])
         colnames(ghg)[1]="year"
         inputLitter.nf <- rbind(ghg,inputLitter.nf)  #North
}


## inputLitter.sf and inputLitter.nf are litter input to Yasso calculations,
## inputLitterFinland is for output

inputLitterFinland <- cbind(simulationPeriod,inputLitter.sf[,2:13]+
                                       inputLitter.nf[,2:13])
if(!calculationRegions[1]) {
inputLitterFinland[,2:13] <- inputLitterFinland[,2:13] - inputLitter.sf[,2:13]
}
if(!calculationRegions[2]) {
inputLitterFinland[,2:13] <- inputLitterFinland[,2:13] - inputLitter.nf[,2:13]
}

colnames(inputLitterFinland) <- c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W",
                                  "fwl_E","fwl_N","cwl_A","cwl_W","nwl_E","cwl_N")


## YASSO CALCULATIONS

## INITIAL STATE

ist <- read.table(YassoInitialStateFile,header=TRUE)   #ist contains Yasso soil C stocks
                             #for years of greenhouse gas inventory

initialState <- data.frame(cbind(t(ist[ist$year==simulationPeriod[1] & ist$region == 1,3:7]),
                 t(ist[ist$year==simulationPeriod[1] & ist$region == 1,8:12]),
                 t(ist[ist$year==simulationPeriod[1] & ist$region == 1,13:17]),
                 t(ist[ist$year==simulationPeriod[1] & ist$region == 2,3:7]),
                 t(ist[ist$year==simulationPeriod[1] & ist$region == 2,8:12]),
                 t(ist[ist$year==simulationPeriod[1] & ist$region == 2,13:17])
                 ))
colnames(initialState) = c("nwl_south", "fwl_south", "cwl_south","nwl_north","fwl_north", "cwl_north")
rownames(initialState) = c("a","w","e","n","h")

# Intial state is per ha , thus totals by multiplication with area
for(i in 1:3) {
  initialState[,i] <- treeInputCheckSouth *
      Areas$area[Areas$region==1&Areas$soil==DeadwoodSoil] * initialState[,i]
  initialState[,i+3] <- treeInputCheckNorth *
      Areas$area[Areas$region==2&Areas$soil==DeadwoodSoil] * initialState[,i+3]
}

# The weather contitions
# mean temperature (C), precipitation (mm), temperature amplitude (C)
# were earlier constant and given in a two-row table:
# row 1 southern Finland, row 2 northern Finland, named in this code
# as e0, and assigned to the value that was defined in MELAToC_defintions.R
#e0 <- YassoEnvironment
# NOW the weather conditions are given through the function (defined in MELAToC_defintions.R)
# weatherFunction(year, region, variable)
   #region 1 = south, 2 = north
   #variable 1 = mean temp, 2 = temp amplitude, 3 = rainfall
# (weather conditions can be either changing or constant)

# In the code, e0 has been replaced by call
# weatherFunction(i, region, variable).
# NOTE that in e0 the variables were in order according to columns:
# temperature, precipitation, temp. amplitude, in weatherFunction the
# "order" of variables is mean temp, temp amplitude, precipitation.

# Climate change affects add on weatherFunction() as before if they are present. 

dt.sf <- DeltaTSouth   # CC temperature function
dt.nf <- DeltaTNorth

dpr.sf <- PrecChgSouth
dpr.nf <- PrecChgNorth

if(ClimateChange) cc = 1.0 else cc = 0.0

lout<-data.frame(matrix(0,nrow=length(simulationPeriod),ncol=7))
colnames(lout)=c("year","nwls","fwls","cwls","nwln","fwln","cwln")
lout[,1] = simulationPeriod

step = 1
zDummy <- c(0,0,0,0,0)
init<-c(0,0,0,0,0)
input <-c(0,0,0,0,0)

zeroNyrs <- rep(0,times=length(simulationPeriod))

Cnwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cnwl.sf)=c("year","A","W","E","N","H")
Cnwl.nf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cnwl.nf)=c("year","A","W","E","N","H")
Cfwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cfwl.sf)=c("year","A","W","E","N","H")
Cfwl.nf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cfwl.nf)=c("year","A","W","E","N","H")

#cwl
diameter = 15

Ccwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Ccwl.sf)=c("year","A","W","E","N","H")
if(calculationRegions[1]) {      #calculationRegions= c(South, North)
#   Ccwl.sf <- data.frame(simulationPeriod,simulationPeriod,
#            simulationPeriod,simulationPeriod,simulationPeriod,simulationPeriod)
#   colnames(Ccwl.sf)=c("year","A","W","E","N","H")

   initcwl.sf = DeadwoodIniDownscaling * initialState[,3]
   if(DeadwoodSpinup) {
     init <- initcwl.sf
     inity = simulationPeriod[1]
     e = c(weatherFunction(inity,1,1)+cc*dt.sf(inity),(1+cc*dpr.sf(inity))
        *weatherFunction(inity,1,3),weatherFunction(inity,1,2))
     spininput =c(mean(inputLitter.sf[,10]),mean(inputLitter.sf[,11]),
           mean(inputLitter.sf[,12]),mean(inputLitter.sf[,13]))
     initcwl.sf <- yasso07(YassoParameters, 10000, e, init,
             c(spininput,0), diameter, zDummy)
   }
   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initcwl.sf} else
          {init <- Ccwl.sf[Ccwl.sf$year==i-1,2:6]}
#      e = c(e0[1,2]+cc*dt.sf(i),(1+cc*dpr.sf(i))*e0[1,3],e0[1,4])
      e = c(weatherFunction(i,1,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,1,3),
                        weatherFunction(i,1,2))
      lout[lout$year==i,4]=sum(inputLitter.sf[inputLitter.sf$year==i,10:13])
      Ccwl.sf[Ccwl.sf$year==i,2:6] <-
       yasso07(YassoParameters, step, e, init,
       c(inputLitter.sf[inputLitter.sf$year==i,10:13],0), diameter, zDummy)
      }
}

Ccwl.nf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Ccwl.nf)=c("year","A","W","E","N","H")
if(calculationRegions[2]) {
#   Ccwl.nf <- data.frame(simulationPeriod,simulationPeriod,
#           simulationPeriod,simulationPeriod,simulationPeriod,simulationPeriod)
#   colnames(Ccwl.nf)=c("year","A","W","E","N","H")

   initcwl.nf=DeadwoodIniDownscaling * initialState[,6]
   if(DeadwoodSpinup) {
       init <- initcwl.nf
       inity = simulationPeriod[1]
       e = c(weatherFunction(inity,1,1)+cc*dt.nf(inity),(1+cc*dpr.nf(inity))
             *weatherFunction(inity,1,3),weatherFunction(inity,1,2))
       spininput =c(mean(inputLitter.nf[,10]),mean(inputLitter.nf[,11]),
                  mean(inputLitter.nf[,12]),mean(inputLitter.nf[,13]))
       initcwl.nf <- yasso07(YassoParameters, 10000, e, init,
                             c(spininput,0), diameter, zDummy)
       }
    
   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initcwl.nf} else
        {init <- Ccwl.nf[Ccwl.nf$year==i-1,2:6]}
#      e = c(e0[2,2]+cc*dt.nf(i),(1+cc*dpr.nf(i))*e0[2,3],e0[2,4])
      e = c(weatherFunction(i,2,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,2,3),
                        weatherFunction(i,2,2))
      lout[lout$year==i,7]=sum(inputLitter.nf[inputLitter.nf$year==i,10:13])
      Ccwl.nf[Ccwl.nf$year==i,2:6] <-
       yasso07(YassoParameters, step, e, init,
       c(inputLitter.nf[inputLitter.nf$year==i,10:13],0), diameter, zDummy)
      }

}

CFinlandMineralSoil = cbind(Cnwl.sf[,"year"],Cnwl.sf[,2:6] + Cnwl.nf[,2:6] +
                 Cfwl.sf[,2:6] + Cfwl.nf[,2:6] + Ccwl.sf[,2:6] +
                 Ccwl.nf[,2:6])
colnames(CFinlandMineralSoil)=c("year","A","W","E","N","H")

write.table(round(lout,4), file=YassoOutputLitter,row.names=FALSE )


CStockPerHa.sf = cbind(Cnwl.sf[,"year"],Cnwl.sf[,2:6]/sf.kankaat +
                   Cfwl.sf[,2:6]/sf.kankaat
                   + Ccwl.sf[,2:6]/sf.kankaat)
colnames(CStockPerHa.sf)=c("year","A","W","E","N","H")

write.table(round(CStockPerHa.sf,4), file=YassoOutputStockPerHaSouth,row.names=FALSE )

CStockPerHa.nf = cbind(Cnwl.nf[,"year"],Cnwl.nf[,2:6]/nf.kankaat +
                     Cfwl.nf[,2:6]/nf.kankaat +
                     Ccwl.nf[,2:6]/nf.kankaat)
colnames(CStockPerHa.nf)=c("year","A","W","E","N","H")


write.table(round(CStockPerHa.nf,4), file=YassoOutputStockPerHaSouth,row.names=FALSE )


change <- diff(rowSums(CFinlandMineralSoil[,2:length(CFinlandMineralSoil[1,])]),lag=1)


yassoresult <- cbind(simulationPeriod,rowSums(inputLitterFinland[,2:length(inputLitterFinland[1,])]),
                          c(change,change[length(change)]))

if(EarlierStart) {
    nolla <- rep(0,times=length(wasteLitter[1,]))
    for(i in 1:(melaTimes[1]-simulationPeriod[1])) {
        wasteLitter <- rbind(nolla,wasteLitter)
        mortLitter <- rbind(nolla,mortLitter)
    }
}

 

# colnames(deadwoodresult) are:
# [1] "year"     "Litter"   "C_change" "wafol"    "wabra"    "wastem"
# [7] "wastump"  "waroots"  "wafr"     "mofol"    "mobra"    "mostem"
# [13] "mostump"  "moroots"  "mofr"
# totLitter = total amount of litter
# C_change  = change of C stock of deadwood
# wa = harvest residue litter
# mo = natural mortality litter
# fol = foliage, bra = branches, stem = stem, stump = stump, roots = coarse roots (> 2m mm)
# fr = fine roots

## NOTE: only "Litter"   "C_change" and "mostem" are meaniningful output,
## others are only because the code was edited from MELAToC_yasso.R - should be
## cleaned up!

deadwoodresult <- cbind(yassoresult,wasteLitter[,2:7],mortLitter[,2:7])
                         
colnames(deadwoodresult) <- c("year","totLitter","C_change",colnames(wasteLitter)[2:7],
                   colnames(mortLitter)[2:7])



deadwoodresult

