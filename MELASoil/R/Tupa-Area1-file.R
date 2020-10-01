##List of defined variables is at the end of this file.

## 1. MELA FILES

MelaInputDirectory <- '~/projections/2011/soil/mela/data/projections/'

MelaStock <- paste(MelaInputDirectory, "Tupa-Area1-stock.csv",sep="")
MelaNatMortality <- paste(MelaInputDirectory,"Tupa-Area1-natmort.txt",sep="")
MelaWaste <- paste(MelaInputDirectory, "Tupa-Area1-waste.txt",sep="")


## 2. YEARS OF STOCK OUPUT FROM MELA

# First and last of them define the time span of simulation
MelaYears <- c(2016,2026,2036,2046,2056,2066)

#Natural mortality and harvest waste are per/year for intervals between stock years,
#Those files have thus one column less the stock file.


## 3. AREAS IN HA AS USED IN MELA CALCULATIONS AS A FUNCTION OF REGION AND SOIL

#region = 1, 2 (South, North)
#soil = 1,2,3 (Upland, drained organic, natural organic)
Areas <- data.frame(rbind(c(1,1,481.7*1000),c(1,2,46.1*1000),c(1,3,1000*17.7),
         c(2,1,481.7*1000),c(2,2,46.1*1000),c(2,3,17.7*1000)))

colnames(Areas)<-c("region","soil","area")


MelaStock1 <- read.csv(MelaStock, header = TRUE);
MelaStock  <- MelaStock1

MelaNatMortality1 <- read.table(MelaNatMortality, header = TRUE)
MelaNatMortality <- MelaNatMortality1

MelaWaste1 <- read.table(MelaWaste, header = TRUE)
MelaWaste <- MelaWaste1
