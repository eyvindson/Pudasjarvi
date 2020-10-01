### DATA FROM MELA

#Yksikko naissa on 1000 m3, jaottelu region,soil,species
#region   1 = south,    2 = north
#soil     1 = mineral   2 = organic
#species  1 = pine      2 = spruce    3 = birch & other decidious
#stock <-  read.csv(paste(path.data,stock.input, sep=""), header=TRUE)
#stock <-  read.csv(MelaStock, header=TRUE)
stock <- MelaStock    

fol.rat.sf <-  FoliageRatioSouth
fol.rat.nf <-  FoliageRatioNorth
    
melaTimes <- MelaYears
melaN = length(melaTimes)

stockC.sf = rep(0,melaN)
stockC.nf = rep(0,melaN)
    
V = rep(0,melaN)
    
for(itimes in 1:melaN) {
  for(iregion in 1:2) {
     sum <- 0.0
     for(isoil in 1:2) {
        for(ispecies in 1:3) {
            # in file stock is per 1000 ha
            stemVol <- 1000.0*stock[stock$region==iregion&stock$soil==isoil
                    &stock$species==ispecies,3+itimes]
            V[itimes] = V[itimes]  + stemVol/1.e6
            #print(stemVol/1000)
            sum <- sum + carbon(functionBEF(iregion,isoil,ispecies,8)*stemVol)
            ## 8 = total biomass
            if(ispecies == 1) {flr <- fol.rat.sf[1]} else {
            if(ispecies == 2) {flr = fol.rat.sf[1]} else{flr = 1}}
            sum = sum +
            carbon(fineroots(flr*functionBEF(iregion,isoil,ispecies,1)*stemVol,
            ispecies))
          }
      }
    if(iregion == 1) stockC.sf[itimes] <-stockC.sf[itimes]+ sum
      else stockC.nf[itimes] = stockC.nf[itimes] + sum
   }
}

# Now changes between the stock years
# 
stockChange.sf = rep(0,(melaN-1))
stockChange.nf = rep(0,(melaN-1))
    

for(itimes in 1:(melaN-1)) {
  stockChange.sf[itimes] = (stockC.sf[itimes+1]-stockC.sf[itimes])
  stockChange.nf[itimes] = (stockC.nf[itimes+1]-stockC.nf[itimes])
}

# apu specifies 1/length of period for melaTimes
apu <- diff(melaTimes)
apu <- 1.0 / apu
    
stockChange.sf <- apu*stockChange.sf
stockChange.nf <-apu*stockChange.nf

changeTimes = as.integer(melaTimes[1:(melaN-1)] + diff(melaTimes)/2)

#simulationPeriod <- min(melaTimes):max(melaTimes)
simulationPeriod <- melaTimes[1]:max(melaTimes)
csf<-approx(changeTimes,stockChange.sf,simulationPeriod,rule=2)$y
cnf<-approx(changeTimes,stockChange.nf,simulationPeriod,rule=2)$y

cstockchange <-cbind(simulationPeriod,csf+cnf)
colnames(cstockchange)=c("year","C_trees")

