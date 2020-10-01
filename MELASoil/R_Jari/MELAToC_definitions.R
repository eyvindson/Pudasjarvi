

## 1. ALL INPUT (MELAINPUT, YEARS, AREAS) ARE READ FROM
##    ONE FILE - THIS THE NEW WAY
##    OR IF THE MELAINPUT (stock, nat. mortality, harvest waste) ARE READ FROM
##    SEPARATE FILES AND YEARS AND AREAS ARE GIVEN AS A SEPARATE INPUT (IN AN INPUT
##    FILE) - THIS THE OLD WAY
##    

NewWay = TRUE
#MelaInputDirectory <- args[2]#'../Projections/' ##CHANGE THIS 
#MELAfil = args[1]#"malusepo_melatulokset_yassoon_ske1.txt" #CHANGE THIS
#MELAInputFile = paste(MelaInputDirectory, MELAfil,sep="")  #Melainput, Melayears and Areas from here
MELAInputFile <-MelaFilePudasjarvi

#Sometimes it is necessary to print out (on console) some information
#for tracking an error in reading the file (probably the file had not standard structure)
if(NewWay) {
InputFileDebug = TRUE
}

if(!NewWay) {
ReadMelaInputOldWayFile = "FRL_AgeClass_P3.5_D13T_FRA-luokitus_jakelu_20190530-file.R"
}


## 1.5 DIRECTORIES FOR INPUT OF WEATHER, LITTER AND INITIAL STOCK -- NOTE THAT
##    FURTHER DEFINITIONS OF INPUT OF LITTER ARE IN 8. LITTER FILE DEFINITIONS

#WeatherInputDirectory =args[2]#'../ParameterData/' 
StockLitterInputDirectory <- '../ParameterData/'



## 2. IT IS POSSIBLE TO ASSESS THE DEADWOOD POOL JUST AS IS DONE BY YASSO
##    (NORMALLY) FOR THE COMBINED DEADWOOD, DOM AND SOM POOL.
##    THIS DONE BY A SEPARATE PIECE OF CODE THAT HAS BEEN ADAPTED
##    FROM MELAToC_yasso.R. THE CALCULATION CAN BE DONE FOR BOTH
##    MINERAL AND ORGANIC SOILS.
##    ALL YASSO- AND OTHER DEFINTIONS BELOW AND ABOVE APPLY ALSO
##    IN THIS CALCULATION.
##    INITIAL VALUE FOR THIS YASSO CALCULATION IS READ FROM THE SAME
##    FILE AS IN NORMAL YASSO CALCULATIONS (SEE DEFINITIONS RELATED TO YASSO)
##    THIS NEEDS TO BE DOWNSCALED BY DeadwoodIniDownscaling. THIS VALUE MAY NEED TO
##    BE ITERATED IN ORDER TO GET START OF SIMULATION RIGHT.

DeadwoodCalculation = FALSE
if(DeadwoodCalculation) {
    DeadwoodSoil = 2        #1 = Mineral soil,2 = Organic soil
    DeadwoodSpinup = FALSE
    DeadwoodIniDownscaling = 0.5
}

# It may be desirable to calculate soil as DOM + SOM, i.e. without DEADWOOD
# this affects MELAToC_yasso.R
# This is 
OnlyDOMSOM = FALSE
if(OnlyDOMSOM) {
OnlyDomSomIniValueReduction = 1.045   #should be 0.981
# OnlyDOMSOM reduces litterinput, inital value need to be corrected accordingly
# This value (0.981) was estimated from FRL MELA litter input
}


## 2. WHICH REGIONS THE CALCULATIONS CONCERN
##   calculationRegions = c(South, North): If South  or North == TRUE
##   IT IS IN THE CALCULATION

#calculationRegions = c(TRUE,TRUE)
calculationRegions = c(TRUE, FALSE)
#calculationRegions = c(FALSE, TRUE)


## 3. FILES FOR LITTER INPUT BEFORE MELA PERIOD AND INITIAL SOIL C STOCK
##    THEY ARE (NORMALLY) FROM CALCULATIONS OF GREEHOUSE GAS INVENTORY

GHGILitterOrganicIF <- "GHGI-litter-organic29.1.2019.txt"

GHGILitterOrganicFile <- paste(StockLitterInputDirectory,GHGILitterOrganicIF,sep="")


## 4. START YEAR OF SIMULATION
##    START, END AND INTERMEDIATE YEARS OF MELA DATA ARE READ FROM MELAInputFile (to MelaYears)
##    SIMULATION MAY START EARLIER THAN FIRST MELA YEAR. THEN PRIOR LITTER INPUT
##    IS READ FROM SEPARATE FILES (SEE 2.).

EarlierStart <- FALSE
SimulationStartYear <- 1975


## 5. DEFINITIONS RELATED TO YASSO

YPath <- '../EXTRA_FILES/'
YassoFile <- paste(YPath,"yasso07.R",sep="")
YassoSharedLibrary <- paste(YPath,"yasso07.so",sep="")

# If weather conditions are constant or change in time - this affects also
# the litter input and initial state file from the GHGI

VariableWeather = TRUE

if(VariableWeather) {
    YassoInitialFi <- "GHGI-soil-carbon-mineral-variable-weather23.1.2019.txt"
    GHGILitterMineralIF <- "GHGI-litter-mineral-variable-weather23.1.2019.txt"
} else {
    YassoInitialFi <- "GHGI-soil-carbon-mineral-constant-weather.txt"
    GHGILitterMineralIF <- "GHGI-litter-mineral-constant-weather.txt"
}
YassoInitialStateFile <- paste(StockLitterInputDirectory, YassoInitialFi,sep="")
GHGILitterMineralFile <- paste(StockLitterInputDirectory,GHGILitterMineralIF,sep="")

Yassox0Factor <- 1.0 #Multiplies initial state of Yasso
YassoSs <- 1  #Multiplies first value of litter

# Parameters
#I got this from Aleksi May 2014
### Scandinavian new, added 7.5.2012
YassoParameters <- c(-0.5172509,-3.551512,-0.3458914,-0.2660175,0.044852223,
           0.0029265443,0.9779027,0.6373951,0.3124745,0.018712098,0.022490378,
           0.011738963,0.00099046889,0.3361765,0.041966144,0.089885026,0.089501545,-0.0022709155,0.17,
           -0.0015,0.17,-0.0015,0.17,-0.0015,0,-2.935411,0,101.8253,
           260,-0.080983594,-0.315179,-0.5173524,0, 0,-0.00024180325,
           0.0015341907, 101.8253,260,-0.5391662,1.18574,-0.2632936,0,0,0)

# Output of Yasso results to files
YassoOutputPath <- '../Results/'
YassoOutputLitter <- paste(YassoOutputPath,"litt.dat",sep="")
YassoOutputStockPerHaSouth <- paste(YassoOutputPath,"CStockPerHa-south.dat",sep="")
YassoOutputStockPerHaNorth <- paste(YassoOutputPath,"CStockPerHa-north.dat",sep="")



## 6. WEATHER

# mean temperature (C), precipitation (mm), temperature amplitude (C)
# WEATHER CAN BE CONSTANT (E.G. MEAN OF A PERIOD) OR CHANGING FROM YEAR TO YEAR
# IN BOTH CASES THE WEATHER IS DEFINED HERE AS A FUNCTION = F(YEAR,REGION,WHICH QUANTITY)
# Weather variables can be given either as running average(averagingRange) or as one year value,
# see comment around averagingRange below

#This is needed MELAToC_drainedsoils.R to asses a possible Climate change effect
constWeather <- rbind(c(1,3.808594, 564.741837, 13.404949),c(2,0.5153716, 508.3148014, 14.6240573))

if(VariableWeather) {
    #weather <- read.csv(paste(WeatherInputDirectory,"PUDASJARVI_SAA_1961_2016.csv",sep=""),header=TRUE)
    weather <- read.csv(WeatherFilePudasjarvi,header=TRUE)
	#Annual values 1960-2017 from Jaakko Heikkinen 18.10.2018
    #columns of weather are  "Year" "meanT_north" "meanT_south" "Ampli_north" "Ampli_south"
    #"Rainfall_north" "Rainfall_south"

    #In the following it is assumed that annual weather values are given from 1960 onwards.
    #If it is later than 1960, you'll need to replace 1960 in the code with the appropriate
    #value.

    endWeatherYear = max(weather[,1])
    
    #Greenhouse Gas Inventory uses 30 year running average (this and previous years)
    #Here this is generalized for any span of running average. If averagingRange == 1,
    #weatherFunction returns current value, max averagingRange = 55 years.

    averagingRange = 30
    
     weatherFunction <- function(year,region,variable) {
         #region 1 = south, 2 = north
         #variable 1 = mean temp, 2 = temp amplitude, 3 = rainfall
         col <- variable * 2
         if(region == 1) {
            col <- col + 1
         }

    if(averagingRange < 1) {
        averagingRange = 1
    }
    if(averagingRange > 55) {
        averagingRange = 55
    }

    year1 <- year
    if(year1 < 1960) {
        year1 = 1960
    }

    if(year1 > endWeatherYear) {
        year1 = endWeatherYear
    }

    if(year1 < 1960+averagingRange - 1) {
        span = year1 - 1960 + 1
    } else {
        span = averagingRange
    }


    lastRow <- which(weather$Year == year1)
    firstRow <- lastRow - span + 1
         
    returnValue <- mean(weather[firstRow:lastRow,col])
    
    returnValue
     }
} else {  #constant weather
  # row 1 southern Finland, row 2 northern Finland
  #YassoEnvironment <- rbind(c(1,3.832523,566.382,13.3571),c(2,0.530517,510.4763,14.57347)) #   mean of 1971-2016
  # from Jaakko Heikkinen 9.10.2018
constWeather <- rbind(c(1,3.808594, 564.741837, 13.404949),c(2,0.5153716, 508.3148014, 14.6240573))
     weatherFunction <- function(year,region,variable) {
                  #region 1 = south, 2 = north
                  #variable 1 = mean temp, 2 = temp amplitude, 3 = rainfall
         if(variable == 1) {
             colm = 2
         } else {
             if(variable == 2) {
                 colm = 4
             } else {
                 colm = 3
             }
         }
         returnValue <- constWeather[region,colm]
         returnValue
     }
}


# Climate change effect to temeperature and precipitation. Temperature: addition to annual mean
# Precipitation: in terms of fractional change from the base year

ClimateChange = FALSE

#- this is RCP2.6 for MALULU scenarios

# Temperature
tcY <-c(2017,2025,2030,2040,2050,2060,2070,2080)
tempAdd.sf<-c(0.0000000,0.3119646,0.4565189,0.7383018,0.9755001,1.1437633,1.2333681,1.3029488)
tempAdd.nf<-c(0.0000000,0.4385175,0.5935624,0.8945198,1.1409522,1.3158980,1.4281085,1.5184558)

DeltaTSouth <- approxfun(tcY,tempAdd.sf,rule=2)
DeltaTNorth <- approxfun(tcY,tempAdd.nf,rule=2)

#Precipitation in terms of fractional change from base year
pcY <- tcY     
precChg.sf <- c(0.00000000,0.01402408,0.02227665,0.03877513,0.05549172,0.06789298,
                0.07087733,0.06717061)
precChg.nf <- c(0.00000000,0.03049505,0.03836610,0.05349924,0.06710153,0.07893003,
                0.08851426,0.09580855)

PrecChgSouth <- approxfun(pcY,precChg.sf,rule=2)
PrecChgNorth <- approxfun(pcY,precChg.nf,rule=2)



## 7. DEFINITIONS RELATED TO ORGANIC SOIL (FORESTED PEATLANDS) CALCULATION

#Shares of peatland types in Finland
#The emssion factors are specified for these peatland types

sh <- read.table("../ParameterData/Peatland-type-shares24.1.2019.txt",header=TRUE)

PeatlandShares <- function(year) {
if(year < 1990) {
    ret <- c(sh[sh$year == 1990,2:6][[1]],sh[sh$year == 1990,2:6][[2]],
             sh[sh$year == 1990,2:6][[3]],sh[sh$year == 1990,2:6][[4]],
             sh[sh$year == 1990,2:6][[5]])
    }
else {if(year > 2016) {
        ret <- c(sh[sh$year == 2016,2:6][[1]],sh[sh$year == 2016,2:6][[2]],
                 sh[sh$year == 2016,2:6][[3]],sh[sh$year == 2016,2:6][[4]],
                 sh[sh$year == 2016,2:6][[5]])
        }
     else {
         ret <- c(sh[sh$year == year,2:6][[1]],sh[sh$year == year,2:6][[2]],
                  sh[sh$year == year,2:6][[3]],sh[sh$year == year,2:6][[4]],
                  sh[sh$year == year,2:6][[5]])
         }
    }
ret
}


## 8. LITTER FILE (INTERMEDIATE and INPUT) DEFINITIONS FOR YASSO AND ORGANIC SOIL CALCULATION

IntermediateLitterDir <- '../WorkFiles/'
MineralSoilSouth <- paste(IntermediateLitterDir,"sf_min.dat", sep="")
MineralSoilNorth <- paste(IntermediateLitterDir,"nf_min.dat", sep="")
OrganicSoilAboveSouth  <- paste(IntermediateLitterDir,"sf_org_abv.dat", sep="")
OrganicSoilAboveNorth  <- paste(IntermediateLitterDir,"nf_org_abv.dat", sep="")
OrganicSoilBelowSouth  <- paste(IntermediateLitterDir,"sf_org_bel.dat", sep="")
OrganicSoilBelowNorth  <- paste(IntermediateLitterDir,"nf_org_bel.dat", sep="")

#Understorey litter from GHGI
UnderstoreyPath <- '../EXTRA_FILES/'
UnderstoreyMineralSouth <-  paste(UnderstoreyPath,"und.nwl.SF.min.csv", sep="")
UnderstoreyMineralNorth <-  paste(UnderstoreyPath,"und.nwl.NF.min.csv", sep="")
UnderstoreyOrganicBelow <-  paste(UnderstoreyPath,"und.org.bel.csv", sep="")  #The same for South and North



## 9. BIOMASS FUNCTIONS AND BEF'S

FunctionPath <- '../EXTRA_FILES/'
#These function files contain many function definitions made by Aleksi Lehtonen
source(paste(FunctionPath,"biomass.r", sep=""))              
source(paste(FunctionPath,"litter.r", sep=""))

#BEF
path.bef <- '../EXTRA_FILES/'  # using BEFs from Juha
BEF.sp <- read.table(paste(path.bef,"pineBEFs10.dat", sep=""), header=TRUE)
BEF.ns <- read.table(paste(path.bef,"spruceBEFs10.dat", sep=""), header=TRUE)
BEF.dc <- read.table(paste(path.bef,"decidBEFs10.dat", sep=""), header=TRUE)

#The BEF tables (BEF.sp, BEF.ns, BEF.dc) are packaged to a function (see below) that returns
#BEF as BEF(region, soil, species, part). Part of the code uses this BEF function.


# foliage ratio of fine roots in southern and northern finland
BEF7.sf <-  read.csv('../EXTRA_FILES/Bef.SF.csv', header=TRUE)
BEF7.nf <-  read.csv('../EXTRA_FILES/Bef.NF.csv', header=TRUE)
FoliageRatioSouth <-  BEF7.sf[1:4, 9]/ BEF7.sf[1:4,8]  #These are vectors of length 4
FoliageRatioNorth <-  BEF7.nf[1:4, 9]/ BEF7.nf[1:4,8]


# Packaging of BEF.sp, BEF.ns, AND BEF.dc to a function

#This function reads BEF from BEF tables BEF.sp, BEF.ns  and BEF.dc when input is:
# region, soil, species, part
# region = 1,2 (southern, northern Finland)
# soil = 1,2,3 mineral, ditched, non-ditched organic soil
# species = 1,2,3 (pine, spruce, birch & other decidious)
# part = 1,2,3,4,5,6,7,8
# foliage           = 1
# branch (living)   = 2
#stem               = 3
#stump              = 4
# roots (coarse)    = 5
#bark               = 6
#deadBranch         = 7
#total              = 8

#ispecies deternines the BEF table
#ipart determines the column. The BEF values are in columns 3 - 10 in the tables
#in the order
#"stemwood" "bark" "liveBranch" "foliage" "deadBranch" "stump" "roots" "total"
#thus ipart translates to columns as follows:
#  ipart    1    2   3   4   5   6   7   8
#  col      6    5   3   8   9   4   7   10
#
# mineral soil is on row 1 (south) or row 4 (north)
# ditched organic soil (isoil=2) is on lines 3  and 6 (south, north)
# non-ditched organic soil (isoil = 3) on lines 2 and 5 (south, north) 
column = c(6,5,3,8,9,4,7,10)
soil_addrow = c(0,2,1)
region_row = c(1,4)

functionBEF <- function(iregion,isoil,ispecies,ipart) {
icol <- column[ipart]
irow <- region_row[iregion] + soil_addrow[isoil]
if(ispecies == 1)
  bef <- BEF.sp[irow,icol]
else {
  if(ispecies == 2)
    bef <- BEF.ns[irow,icol]
  else
    bef <- BEF.dc[irow,icol]
  }

    bef
}





