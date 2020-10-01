##List of defined variables is at the end of this file.

## 1. MELA FILES

MelaInputDirectory <- '../Projections/'

MelaStock <- paste(MelaInputDirectory, "FRL_AgeClass_P3.5_D13T_FRA-luokitus_jakelu_20190530-stock.csv",sep="")
MelaNatMortality <- paste(MelaInputDirectory, "FRL_AgeClass_P3.5_D13T_FRA-luokitus_jakelu_20190530-natmort.txt",sep="")
MelaWaste <- paste(MelaInputDirectory, "FRL_AgeClass_P3.5_D13T_FRA-luokitus_jakelu_20190530-waste.txt",sep="")


## 2. YEARS OF STOCK OUPUT FROM MELA

# First and last of them define the time span of simulation
MelaYears <- c(2011,2016,2021,2026,2031,2041,2051,2061)

#Natural mortality and harvest waste are per/year for intervals between stock years,
#Those files have thus one column less the stock file.


## 3. AREAS IN HA AS USED IN MELA CALCULATIONS AS A FUNCTION OF REGION AND SOIL

#region = 1, 2 (South, North)
#soil = 1,2,3 (Upland, drained organic, natural organic)
Areas <- data.frame(rbind(c(1,1,8808.4*1000),c(1,2,2204.6*1000),c(1,3,1000*373.4),
         c(2,1,7128.5*1000),c(2,2,2239.2*1000),c(2,3,1289.2*1000)))

colnames(Areas)<-c("region","soil","area")


MelaStock1 <- read.csv(MelaStock, header = TRUE);
MelaStock  <- MelaStock1

MelaNatMortality1 <- read.table(MelaNatMortality, header = TRUE)
MelaNatMortality <- MelaNatMortality1

MelaWaste1 <- read.table(MelaWaste, header = TRUE)
MelaWaste <- MelaWaste1
