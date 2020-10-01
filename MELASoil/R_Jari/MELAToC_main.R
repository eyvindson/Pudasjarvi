## This  "main program" runs transformation of MELA output to carbon balance of forests for
## the whole country. The 
## MELA output consists of three parts: data on stock (amount of stemwood),
## harvest residues (biomass) and natural mortality (amount of stemwood).
## These data are given for two areas: Southern Finland (Forestry Centers 1 - 10) and
## Northern Finland (Forestry Centres 11-13).
## Areas of mineral soils and drained, forested organic soils are also input from MELA.
##
## Calculation takes place by pieces of R-code through source(" "). As a result, the variables
## in any part of the code are visible in other parts. BE CAREFUL NOT TO REDEFINE VARIABLES
## IF YOU ADD CODE TO THIS SYSTEM!
##
## All necessary input & definition of variables for calculations are given in the definition
## file MELAToC_definitions.R. The variables are kind of global variables for calculations.
## List of defined variables are at the end of file. The naming convention is capitalize
## each word & concatenate, e.g. MelaInputDirectory.
##
## MELAToC_volToCInput.R calculates litter input from standing trees. It writes intermediate
## files that Yasso and organic soil calculation will read. Stores results in data frame litt.
##
## MELAToC_stock.R calculates C balance of trees, results are in data frame cstockchange.
##
## MELATOc_yasso.R: Yasso for mineral soils, results are in data frame yassoresult.
##
## MELAToC_drainedsoils.R: calculations for derained organic soils, results are in
## data frame organic.
##
## Final results are in data frames cFinland (C) and CO2result (CO2, IPCC sign convention)
## Contents of the output data frames are explained at the end of this file.

#For Pudasjarvi, take file name from the command line
 args = commandArgs(trailingOnly=TRUE)
 if (length(args)==0) {
   stop("Missing file name", call.=FALSE)
 } else if (length(args)==1) {
    # print(args[1])
}

MelaFilePudasjarvi <- args[1]
WeatherFilePudasjarvi <- args[2]

##Define global variables for Mela yasso input file and Weather file
##(Cannot give parameters source, only Rscript allows for that
print("Yasso file")
print(MelaFilePudasjarvi)
print("Weather file")
print(WeatherFilePudasjarvi)

source("MELAToC_definitions.R")

if(NewWay) {
    source("MELAToC_file.R")
} else {
    source(ReadMelaInputOldWayFile)
}

if(!DeadwoodCalculation) {
source("MELAToC_volToCInput.R")    #litt
source("MELAToC_stock.R")          #cstockchange
source("MELAToC_yasso.R")          #yassoresult


## organic soils need mean stem volume (for changing emission factors). They are evaluated
## here and given as input
#stockcvs <-  data.frame(read.csv(MelaStock, header=TRUE))
stockcvs <-  MelaStock
MelaMeanV = 1:length(MelaYears)
for(i in 1:length(MelaYears)) {

MelaMeanV[i] = (sum(stockcvs[stockcvs$region==1&stockcvs$soil==2,][3+i])/
          (Areas[Areas$region==1&Areas$soil==2,3]/1000)+
        sum(stockcvs[stockcvs$region==2&stockcvs$soil==2,][3+i])/
       (Areas[Areas$region==2&Areas$soil==2,3]/1000))/2
}

source("MELAToC_drainedsoils.R") #organic

## OUTPUT

treel.org.abv.sf<-  read.table(OrganicSoilAboveSouth, header=TRUE)
treel.org.abv.nf<-  read.table(OrganicSoilAboveNorth, header=TRUE)
                               

tl.org.abv <- Areas$area[Areas$region==1&Areas$soil==2]*rowSums(treel.org.abv.sf[,2:9])+
                  Areas$area[Areas$region==1&Areas$soil==2]*rowSums(treel.org.abv.nf[,2:9])

if(!EarlierStart) {
    cFinland <- cbind(yassoresult[,1],yassoresult[,2]*1.0e-6,organic[,2]*1.0e-6,
           cstockchange[,2]*1.0e-6,yassoresult[,3]*1.0e-6,organic[,3]*1.0e-6,
           (yassoresult[,3]+cstockchange[,2]+organic[,3])*1.0e-6,1.0e-6*yassoresult[,4:15],
           1.0e-6*organic[,4:9],1.0e-6*tl.org.abv)
} else {
    nolla <- rep(0,times=MelaYears[1]-SimulationStartYear)
    total <- (yassoresult[yassoresult$year>=MelaYears[1],3]+cstockchange[,2]
        +organic[organic$year>=MelaYears[1],3])*1.0e-6
    cFinland <- cbind(yassoresult[,1],yassoresult[,2]*1.0e-6,organic[,2]*1.0e-6,
           c(nolla,cstockchange[,2]*1.0e-6),yassoresult[,3]*1.0e-6,organic[,3]*1.0e-6,
           c(nolla,total),
           1.0e-6*yassoresult[,4:15], 1.0e-6*organic[,4:9],c(nolla,1.0e-6*tl.org.abv))
}

colnames(cFinland) <- c(colnames(yassoresult)[1:2],colnames(organic)[2],colnames(cstockchange)[2],
                        colnames(yassoresult)[3], colnames(organic)[3],"total",
                        colnames(yassoresult)[4:15], colnames(organic)[4:9], "otree_abv")

CO2result <- as.data.frame(cbind(cFinland[,1],-44*cFinland[,4:6]/12,
                              -44*(cFinland[,5]+cFinland[,6])/12,-44*cFinland[,7]/12))
colnames(CO2result) <- c("year","C_trees", "C_mineral", "C_organic","C_soil","Total")
CO2result
write.table(CO2result,file="CO2result.dat",row.names=FALSE,sep=" ")

} else {    #Only deadwood

source("MELAToC_deadwood.R")          #yassoresult

deadwood <- as.data.frame(cbind(yassoresult[,1],1.0e-6*44*yassoresult[,2]/12,
                                -1.0e-6*44*yassoresult[,3]/12))
colnames(deadwood) <- c("year","litter", "C_stockCh")
deadwood
write.table(deadwood,file="deadwood.dat",row.names=FALSE,sep=" ")
}

#plot(CO2result$year,CO2result$C_mineral,type="l",main="Soil stock change FRL_AgeClass_P3",ylim=c(-10,5),
#     ylab = "C stock change (Mt CO2 / year)",xlab="",xlim=c(2010,2070))
#abline(0,0)
