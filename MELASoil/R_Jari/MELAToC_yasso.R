source(YassoFile)
dyn.load(YassoSharedLibrary)
    
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
# rows are years of simulation
# unit is t C / ha

#nf.min is as sf_min but for northern Finland

# Initial states for upland soils (southern Finland, northern Finland) are
# read from a file.
# Columns are: nwl_south" "fwl_south" "cwl_south" "nwl_north" "fwl_north" "cwl_north"
# Rows are A W E N H
###############################################################################

melaWaste <- MelaWaste
melaNatMort <- MelaNatMortality

#Litter from standing trees is read from (intermediate) files.
# The columns in these files are
# "year" "nwl_A" "nwl_W" "nwl_E" "nwl_N" "fwl_A" "fwl_W" "fwl_E" "fwl_N"
# (nwl = non woody litter, fwl = fine woody litter (no coarse woody litter cwl
# from standing trees) )
sf.min <- read.table(MineralSoilSouth, header=TRUE)
nf.min <- read.table(MineralSoilNorth, header=TRUE)
    
# Litter from standing trees is per ha (to be consistent with Aleksi's calculations)
# Totals by multiplication with area
    
#MELA areas
sf.kankaat = Areas$area[Areas$region==1&Areas$soil==1]
nf.kankaat = Areas$area[Areas$region==2&Areas$soil==1]
    
#There are 8 data colums in sf.min and nf.min
for(i in 2:9) {
sf.min[,i] <- sf.kankaat * sf.min[,i]
nf.min[,i] <- nf.kankaat * nf.min[,i]
}


#Sum all litter inputs (melaWaste, melaNatMort and min (sf na nf)

melaTimes <- MelaYears

melaN = length(melaTimes)
melaNlitter = melaN - 1

## First melaWaste, that is cutting residues

#wastelitter for output
wastelitt<-matrix(data=0,nrow=melaNlitter,ncol=6)


#fol.rat is for correcting estimates of fine roots made on the
#basis of biomass calculated with the aid of Repola et al. functions and
#using BEFs based on Marklund's eqns
fol.rat.sf <-  FoliageRatioSouth
fol.rat.nf <-  FoliageRatioNorth

for(iregion in 1:2) {
  sumlitter.nwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
  sumlitter.fwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
  sumlitter.cwl <-matrix(data=0,nrow=melaNlitter,ncol=4)
  #rows = years of MELA input

  #column = A W E N
  for(ispecies in 1:3) {
      for(ipart in 1:5) {
        for(itimes in 1:melaNlitter) {
   #wastelitt is summed here + addition from fine roots that is calculated
   #with the aid of foliage (ipart == 1)

          wastelitt[itimes,ipart] = wastelitt[itimes,ipart] +
          carbon(melaWaste[melaWaste$region==iregion&melaWaste$soil==1
            &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes])
            
          if(ipart == 1) { #fine roots also here with foliage
            sumlitter.nwl[itimes,] = sumlitter.nwl[itimes,] +
                foliage.AWEN(carbon(melaWaste[melaWaste$region==iregion&melaWaste$soil==1
                &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes]),ispecies)
            if(ispecies == 1) {
                if(iregion==1) flr <- fol.rat.sf[1] else flr <- fol.rat.nf[1]} else {
                    if(ispecies == 2) {
                        if(iregion==1) flr <- fol.rat.sf[2] else flr <- fol.rat.nf[2]} else{
                              flr = 1}
                  }
            sumlitter.nwl[itimes,] = sumlitter.nwl[itimes,] +
                  fineroot.AWEN(carbon(fineroots(flr*melaWaste[melaWaste$region==iregion&
                  melaWaste$soil==1&melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes],
                                                 ispecies)),ispecies)
            wastelitt[itimes,6] = wastelitt[itimes,6] +
             carbon(fineroots(flr*melaWaste[melaWaste$region==iregion&
             melaWaste$soil==1&melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes],
             ispecies))
          }
          if(ipart == 2) {
            sumlitter.fwl[itimes,] = sumlitter.fwl[itimes,] +
                 branches.AWEN(carbon(melaWaste[melaWaste$region==iregion&melaWaste$soil==1
                 &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes]))/
                   sum(branches.AWEN(1))
          }
          if(ipart == 5) { #roots are also fwl
            sumlitter.fwl[itimes,] = sumlitter.fwl[itimes,] +
                 stem.AWEN(carbon(melaWaste[melaWaste$region==iregion&melaWaste$soil==1
                 &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes]),ispecies)/
                   sum(stem.AWEN(1,ispecies))
          }
          if(ipart == 3 | ipart == 4) { #stem and stump are cwl
             sumlitter.cwl[itimes,] = sumlitter.cwl[itimes,] +
                  stem.AWEN(carbon(melaWaste[melaWaste$region==iregion&melaWaste$soil==1
                  &melaWaste$species==ispecies&melaWaste$part==ipart,4+itimes]),ispecies)/
                    sum(stem.AWEN(1,ispecies))
           }
        }
      }
    if(iregion == 1) {
      sumlitter.nwl.sf <- sumlitter.nwl
      sumlitter.fwl.sf <- sumlitter.fwl
      sumlitter.cwl.sf <- sumlitter.cwl
    } else {
      sumlitter.nwl.nf <- sumlitter.nwl
      sumlitter.fwl.nf <- sumlitter.fwl
      sumlitter.cwl.nf <- sumlitter.cwl}
    }
  }


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
          stemVol <- melaNatMort[melaNatMort$region==iregion&melaNatMort$soil==1
                     &melaNatMort$species==ispecies,3+itimes]
          for(ipart in 1:5) {
            if(ipart == 1)  {#with foliage also fine roots
             sumlitter.nwl[itimes,] = sumlitter.nwl[itimes,] +
              foliage.AWEN(carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol)
             
              if(ispecies == 1) {
                flr <- fol.rat.sf[1]} else {
                if(ispecies == 2) {
                  flr = fol.rat.sf[1]} else{flr = 1}
              }
             sumlitter.nwl[itimes,] = sumlitter.nwl[itimes,] +
             fineroot.AWEN(carbon(fineroots(flr*functionBEF(iregion,1,ispecies,ipart)*stemVol,
                                            ispecies)),ispecies)
             mortlitt[itimes,6] = mortlitt[itimes,6] +
             carbon(fineroots(flr*functionBEF(iregion,1,ispecies,ipart)*stemVol,
                                                       ispecies))             

           } #if(ipart ==1
           if(ipart == 2) {
             sumlitter.fwl[itimes,] = sumlitter.fwl[itimes,] +
             branches.AWEN(carbon((functionBEF(iregion,1,ispecies,ipart)+
             functionBEF(iregion,1,ispecies,7))*stemVol))
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
             carbon((functionBEF(iregion,1,ispecies,ipart)+
              functionBEF(iregion,1,ispecies,7))*stemVol)
           }
           if(ipart == 5) {  #Roots are also fwl
             sumlitter.fwl[itimes,] = sumlitter.fwl[itimes,] +
             stem.AWEN(carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol)
           }
           if( ((ipart == 3) & !OnlyDOMSOM) | ipart == 4) {       #stem and stump are cwl
             sumlitter.cwl[itimes,] = sumlitter.cwl[itimes,] + #OnlyDOMSOM = not Deadwood
             stem.AWEN(carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,ipart)*stemVol)
             
             if(ipart ==3) {sumlitter.cwl[itimes,] = sumlitter.cwl[itimes,] +
             stem.AWEN(carbon(functionBEF(iregion,1,ispecies,6)*stemVol),ispecies)
             mortlitt[itimes,ipart] = mortlitt[itimes,ipart] +
               carbon(functionBEF(iregion,1,ispecies,6)*stemVol)       #Bark
             }

           }
        }#for(ipart     
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


#Note change of first value of litter with coefficient YassoSs

for(i in 1:4) {
sumlitter.nwl.sf[1,i] = sumlitter.nwl.sf[1,i]*YassoSs
sumlitter.fwl.sf[1,i] = sumlitter.fwl.sf[1,i]*YassoSs
sumlitter.cwl.sf[1,i] = sumlitter.cwl.sf[1,i]*YassoSs
sumlitter.nwl.nf[1,i] = sumlitter.nwl.nf[1,i]*YassoSs
sumlitter.fwl.nf[1,i] = sumlitter.fwl.nf[1,i]*YassoSs
sumlitter.cwl.nf[1,i] = sumlitter.cwl.nf[1,i]*YassoSs
}



#This is set again below, depending whther EarlierStart is active. Here for MELA litter
simulationPeriod <- melaTimes[1]:max(melaTimes)

#Approximation of wastelitter and mortlitter that are a step function = constant during 10 year
#periods: the value is at the middle of the period

#approxTimes = value time of litter input. Best (must) be between times when MELA stock values
#have been given. The litter from harvest (waste) and natural mortality is constant for periods
# between MELA stock values - most natural
#is thus in the middle of the period. rule = 2 gives constant value outside approxTimes, thus
#litter input from harvest and natural mortality = constant in the first half of the first MELA period.
#Litter from standing trees (sf.min and nf.min) change from the start because it is calculated with
#the aid of stock values that change between melaTimes.

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
colnames(wasteLitter) <- c("year","lwastefol","lwastebra","lwastestem","lwastestump","lwasteroots",
                           "lwastefr")
colnames(mortLitter) <- c("year","lmortfol","lmortbra","lmortstem","lmortstump","lmortroots","lmortfr")



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


## UNDERSTOREY LITTER

undSf <-  read.csv(UnderstoreyMineralSouth, header=TRUE)
undNf <-  read.csv(UnderstoreyMineralNorth, header=TRUE)
# understorey litter (constant value) is already AWEN divided

#NOTE IF NO TREELITTEINPUT (treeInputCheck...), ALSO UNDERSTOREY LITTERINPUT = 0
for(i in 1:length(simulationPeriod)) {
      inputLitter.sf[i,2:5] <- inputLitter.sf[i,2:5] +
          sf.kankaat * treeInputCheckSouth * t(undSf)[1:4]
      inputLitter.nf[i,2:5] <- inputLitter.nf[i,2:5] +
          nf.kankaat * treeInputCheckNorth * t(undNf)[1:4]
    }



## SIMULATION MAY START EARLIER THAN AT FIRST YEAR OF melaTimes (=MelaYears), LITTER INPUT
## IS TAKEN THEN FROM LITTER USED IN GREENHOUSE GAS INVENTORY IN MINERAL SOIL (YASSO)

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

if(OnlyDOMSOM) {
initialState <- OnlyDomSomIniValueReduction * initialState
}

# Intial state is per ha , thus totals by multiplication with area
for(i in 1:3) {
  initialState[,i] <- treeInputCheckSouth * sf.kankaat * initialState[,i]
  initialState[,i+3] <- treeInputCheckNorth * nf.kankaat * initialState[,i+3]
}








# Weather conditions are given through the function (defined in MELAToC_defintions.R)
# weatherFunction(year, region, variable)
   #region 1 = south, 2 = north
   #variable 1 = mean temp, 2 = temp amplitude, 3 = rainfall


# Climate change adds on weatherFunction() values if it is present.






if(ClimateChange) cc = 1.0 else cc = 0.0


dt.sf <- DeltaTSouth      # CC temperature function
dt.nf <- DeltaTNorth





dpr.sf <- PrecChgSouth    # CC precipitation function

















dpr.nf <- PrecChgNorth




lout<-data.frame(matrix(0,nrow=length(simulationPeriod),ncol=7))   #for litter output
colnames(lout)=c("year","nwls","fwls","cwls","nwln","fwln","cwln")
lout[,1] = simulationPeriod

step = 1
zDummy <- c(0,0,0,0,0)
init<-c(0,0,0,0,0)
input <-c(0,0,0,0,0)

zeroNyrs <- rep(0,times=length(simulationPeriod))

#nwl
diameter <- 0

Cnwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cnwl.sf)=c("year","A","W","E","N","H")
if(calculationRegions[1]) {      #calculationRegions= c(South, North)
   initnwl.sf = Yassox0Factor * initialState[,1]




   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initnwl.sf} else
            {init <- Cnwl.sf[Cnwl.sf$year==i-1,2:6]}


      e = c(weatherFunction(i,1,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,1,3),
            weatherFunction(i,1,2))
      
      lout[lout$year==i,2]=sum(inputLitter.sf[inputLitter.sf$year==i,2:5])
      Cnwl.sf[Cnwl.sf$year==i,2:6] <-
      yasso07(YassoParameters, step, e, init,
              c(inputLitter.sf[inputLitter.sf$year==i,2:5],0), diameter, zDummy)
    }

}

Cnwl.nf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cnwl.nf)=c("year","A","W","E","N","H")
if(calculationRegions[2]) {




   initnwl.nf = Yassox0Factor * initialState[,4]

   for(i in simulationPeriod) {
    if (i==simulationPeriod[1]) {init <- initnwl.nf} else
         {init <- Cnwl.nf[Cnwl.nf$year==i-1,2:6]}


    e = c(weatherFunction(i,2,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,2,3),
          weatherFunction(i,2,2))

    lout[lout$year==i,5]=sum(inputLitter.nf[inputLitter.nf$year==i,2:5])
    Cnwl.nf[Cnwl.nf$year==i,2:6] <-
         yasso07(YassoParameters, step, e, init,
         c(inputLitter.nf[inputLitter.nf$year==i,2:5],0), diameter, zDummy)
    }
}
#fwl
diameter = 2

Cfwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cfwl.sf)=c("year","A","W","E","N","H")
if(calculationRegions[1]) {      #calculationRegions= c(South, North)




   initfwl.sf=Yassox0Factor * initialState[,2]

    
   for(i in simulationPeriod) {
        if (i==simulationPeriod[1]) {init <- initfwl.sf} else
            {init <- Cfwl.sf[Cfwl.sf$year==i-1,2:6]}


      e = c(weatherFunction(i,1,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,1,3),
            weatherFunction(i,1,2))

        lout[lout$year==i,3]=sum(inputLitter.sf[inputLitter.sf$year==i,6:9])
         Cfwl.sf[Cfwl.sf$year==i,2:6] <-
         yasso07(YassoParameters, step, e, init,
         c(inputLitter.sf[inputLitter.sf$year==i,6:9],0), diameter, zDummy)
      }
}

Cfwl.nf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Cfwl.nf)=c("year","A","W","E","N","H")
if(calculationRegions[2]) {




    initfwl.nf=Yassox0Factor * initialState[,5]

    
   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initfwl.nf} else
          {init <- Cfwl.nf[Cfwl.nf$year==i-1,2:6]}


      e = c(weatherFunction(i,2,1)+cc*dt.sf(i),(1+cc*dpr.sf(i))*weatherFunction(i,2,3),
            weatherFunction(i,2,2))
      
      lout[lout$year==i,6]=sum(inputLitter.nf[inputLitter.nf$year==i,6:9])
      Cfwl.nf[Cfwl.nf$year==i,2:6] <-
       yasso07(YassoParameters, step, e, init,
       c(inputLitter.nf[inputLitter.nf$year==i,6:9],0), diameter, zDummy)
    }
}

#cwl
diameter = 15

Ccwl.sf <- data.frame(simulationPeriod,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs,zeroNyrs)
colnames(Ccwl.sf)=c("year","A","W","E","N","H")
if(calculationRegions[1]) {      #calculationRegions= c(South, North)




   initcwl.sf = Yassox0Factor * initialState[,3]
    
   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initcwl.sf} else
          {init <- Ccwl.sf[Ccwl.sf$year==i-1,2:6]}


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




   initcwl.nf=Yassox0Factor * initialState[,6]
    
   for(i in simulationPeriod) {
      if (i==simulationPeriod[1]) {init <- initcwl.nf} else
        {init <- Ccwl.nf[Ccwl.nf$year==i-1,2:6]}


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


write.table(round(CStockPerHa.nf,4), file=YassoOutputStockPerHaNorth,row.names=FALSE )


change <- diff(rowSums(CFinlandMineralSoil[,2:length(CFinlandMineralSoil[1,])]),lag=1)

# Here yassoresult = c(years, total litter input, C stock change)

yassoresult <- cbind(simulationPeriod,rowSums(inputLitterFinland[,2:length(inputLitterFinland[1,])]),
                          c(change,change[length(change)]))

#Litter from standing trees for output
littTrees <- apply(sf.min[,2:8],1,sum) + apply(nf.min[,2:8],1,sum)


#In case of earlier start than in MELA input augment output litters with
#zeroes at the beginning
if(EarlierStart) {
    nolla <- rep(0,times=length(wasteLitter[1,]))
    for(i in 1:(melaTimes[1]-simulationPeriod[1])) {
        wasteLitter <- rbind(nolla,wasteLitter)
        mortLitter <- rbind(nolla,mortLitter)
    }
    littTrees <- c(rep(0,times=(melaTimes[1]-simulationPeriod[1])),littTrees)
}

#Adding understorey litter to output litter. Understorey litter input is constant,
#divided into AWEN

littUnd <- rep(sf.kankaat*sum(undSf), times=length(simulationPeriod)) +
               rep(nf.kankaat*sum(undNf), times=length(simulationPeriod))


    
yassoresult <- cbind(yassoresult, littTrees, wasteLitter[,2:7],mortLitter[,2:7],littUnd)
                         
colnames(yassoresult) <- c("year","littInput","CStockChange","littTrees",colnames(wasteLitter)[2:7],
                   colnames(mortLitter)[2:7], "littUnd")

#yassoresult contains the following columns (unit t C, except year)
# year littInput CStockChange littTrees  wastefol wastebra wastestem wastestump
# wasteroots wastefr mortfol mortbra mortstem  mortstump mortroots
# mortfr littUnd

#littInput     = total amount of litter
#CStockChange  = change in C stock
#littTrees     = litter from living trees
#waste         = harvest residue litter
#mort          = litter from natural mortality
# fol = foliage, bra = branches, stem = stem, stump = stump, roots = coarse roots (> 2m mm)
# fr = fine roots




