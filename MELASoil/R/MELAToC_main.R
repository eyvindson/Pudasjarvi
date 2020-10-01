## This  "main program" runs transformation of MELA output to carbon balance of forests.
## MELA output consists of three parts (input files): data on stock (amount of stemwood),
## harvest residues (biomass) and natural mortality (amount of stemwood).
##
## Calculation takes place by pieces of R-code through source(" "). As a result, the variables
## in any part of the code are visible in other parts. BE CAREFUL NOT TO REDEFINE VARIABLES
## IF YOU ADD CODE TO THIS SYSTEM!
##
## All necessary input & definition of variables for calculations are given in the definition
## file MELAToC_definitions.R. The variables are thus kind of global variables for calculations.
## List of defined variables are at the end of that file. The naming convention is capitalize
## each word & concatenate, e.g. MelaInputDirectory. -- list of defined variables is outdated at the moment!
##
## Two versions of MELA input can be handled:
## NewWay: all data is in one file, MELAToC_file.R reads these data.
## !NewWay (i.e. old way): data for tree stock, natural mortality and harvest waste are
## in separate files, they and area information and calculation times are read
## from file specified by the variable ReadMelaInputOldWayFile.
##
## The components of calculation are:
## MELAToC_volToCInput.R calculates litter input from standing trees. It writes intermediate
## files that Yasso and organic soil calculation will read. It stores results in data frame litt.
##
## MELATOc_yasso.R: Yasso for mineral soils, results are in data frame yassoresult.
##
## MELAToC_drainedsoils.R: calculations for forested drained organic soils (peatlands),
## results are in data frame organic.
##
## Final results are in data frames cFinland (C) and CO2result (CO2, IPCC sign convention)

args = commandArgs(trailingOnly=TRUE)
# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  #stop("At least one argument must be supplied (input file).n", call.=FALSE)
  args[1] ="832907_PUDASJARVI_P30061_SKE4_4.txt"
  args[2] = '/scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/'
  args[3] = '/scratch/project_2003225/GIT/Yasso/Output_files_Yasso/'
} else if (length(args)==1) {
  # default output file
  args[2] = '/scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/'
  args[3] = '/scratch/project_2003225/GIT/Yasso/Output_files_Yasso/'

}
#args[1] ="832907_PUDASJARVI_P30061_SKE4_4.txt"
#args[2] = '/scratch/project_2003225/GIT/Yasso/Input_files_Yasso/PUDASJARVI_P30061_2010_2060/'
	
	
source("MELAToC_definitions.R")

if(NewWay) {
    source("MELAToC_file.R")
} else {
    source(ReadMelaInputOldWayFile)
}

                                   #Names of dataframes that contain the results
                                   #See R files for documentation of their contents
if(!DeadwoodCalculation) {
source("MELAToC_volToCInput.R")    #litt

#litt contains the following columns:
#year, mlitfol,mlitbra,mlitstem,mlitstump,mlitroots,mlitfr,olitfol,olitbra,
#olitstem,olitstump,olitroots,olitfr
#unit = t C
#mlit = mineral soil
#olit = organic soil
# fol = foliage, bra = branches, stem = stem, stump = stump, roots = coarse roots (> 2m mm)
# fr = fine roots


source("MELAToC_yasso.R")          #yassoresult

#yassoresult contains the following columns (unit t C, except year) 
# year littInput CStockChange littTrees  wastefol wastebra wastestem wastestump
# wasteroots wastefr mortfol mortbra mortstem  mortstump mortroots
# mortfr littUnd

# littInput = total amount of litter
# CStockChange = change in C stock of soil (sink = +)
# littTrees = amount of litter from standing trees
# Rest of the colums are amounts of litter from waste: lwaste* and
# natural mortality: lmort*
# fol = foliage, bra = branches, stem = stems, stump = stump
# stump = stump, roots = roots, fr = fine roots
# AND littUnd = amount of litter produced by understorey vegetation


source("MELAToC_drainedsoils.R") #organic

#organic contains the following columns ((unit t C, except year):
# year orgLitter CStockChange  lwaste_bel lwaste_abvl mort_bel
# lmort_abv  ltree_bel  lund_bel

# orgLitter = total amount of litter that was used in the calculation of
# CStockChange = change in C stock of soil  (sink = +)
# lwaste_bel = amount of litter from harvest residues, below-ground
# lwaste_abv = amount of litter from harvest residues, above-ground
# lmort_bel = amount of litter from natural mortality, below-ground
# lmort_abv = amount of litter from natural mortality, above-ground
# tree_bel = amount of litter from standing trees, below-ground
# lund_bel = amount of litter from understorey vegetation, below-ground


## OUTPUT

treel.org.abv.sf<-  read.table(OrganicSoilAboveSouth, header=TRUE)
treel.org.abv.nf<-  read.table(OrganicSoilAboveNorth, header=TRUE)
                               

tl.org.abv <- Areas$area[Areas$region==1&Areas$soil==2]*rowSums(treel.org.abv.sf[,2:9])+
                  Areas$area[Areas$region==1&Areas$soil==2]*rowSums(treel.org.abv.nf[,2:9])

if(!EarlierStart) {
    cFinland <- cbind(yassoresult[,1],yassoresult[,2]*1.0e-6,organic[,2]*1.0e-6,
           yassoresult[,3]*1.0e-6,organic[,3]*1.0e-6,
           1.0e-6*yassoresult[,4:15],
           1.0e-6*organic[,4:9],1.0e-6*tl.org.abv)
} else {
    nolla <- rep(0,times=MelaYears[1]-SimulationStartYear)
    cFinland <- cbind(yassoresult[,1],yassoresult[,2]*1.0e-6,organic[,2]*1.0e-6,
           yassoresult[,3]*1.0e-6,organic[,3]*1.0e-6,
           1.0e-6*yassoresult[,4:15], 1.0e-6*organic[,4:9],c(nolla,1.0e-6*tl.org.abv))
}

colnames(cFinland) <- c(colnames(yassoresult)[1],"minLitter","orgLitter",
                        "minCstockchange", "orgCstockchange",
                        colnames(yassoresult)[4:15], colnames(organic)[4:9], "otree_abv")

#cFinland  unit = Tg C, columns:
#year
#minLitter   = amount of litter on mineral soils
#orgLitter   = amount of litter on organic soils
#minCstockchange  = C stock change on mineral soils (sink = +)
#orgCstockchange  = C stock change on organic soils (sink = +)
#littTrees   = amount of litter from trees

#lwaste      = litter from harvest residues
#lmort       = litter from natural mortality
#bel         = below ground
#abv         = above ground
# fol = foliage, bra = branches, stem = stems, stump = stump
# stump = stump, roots = roots, fr = fine roots

#ltree_bel   = litter from trees on organic soil below ground
#lund_bel    = litter from understorey plant on organic soil below ground
#otree_abv   = litter from trees on organic soil above ground



CO2result <- as.data.frame(cbind(cFinland[,1],-44*cFinland[,4:5]/12,
                              -44*(cFinland[,4]+cFinland[,5])/12))
colnames(CO2result) <- c("year","C_mineral", "C_organic","C_soil")


output_file =  paste(args[3], substr(args[1],start=1,stop = nchar(args[1])-4),sep="") 
write.table(CO2result, paste(output_file,"_out.csv",sep=""), sep = ";",)
output_file_litter = gsub("Output_files_Yasso", "Litter_output_files_yasso", output_file)
write.csv(sf.min,paste(output_file_litter,"_litter.csv",sep=""), row.names = TRUE)

#columns of CO2result, IPCC sign convention
# year
# C_mineral    = change of C stock on mineral soils Tg CO2 /year
# C_organic    = change of C stock on organic soils Tg CO2 /year
# C_soil       = C_mineral + C_organic

} else {    #Only deadwood

source("MELAToC_deadwood.R")   

#unit = t C    
#litter =     amount of litter
#C_stockCh =  change of C stock of deadwood
deadwood <- as.data.frame(cbind(yassoresult[,1],1.0e-6*44*yassoresult[,2]/12,
                                -1.0e-6*44*yassoresult[,3]/12))
colnames(deadwood) <- c("year","litter", "C_stockCh")

deadwood
}
