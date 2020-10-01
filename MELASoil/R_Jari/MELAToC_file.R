countLines = function(filepath) {
  con = file(filepath, "r")
  n = 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    n =  n + 1
  }
  close(con)
  if(n > 0)
   n - 1
  	else
  	n 
}


readAllLines = function(filepath, lines) {
  con = file(filepath, "r")
  i = 0
  while ( TRUE ) {
    line = readLines(con, n = 1)
    if ( length(line) == 0 ) {
      break
    }
    i = i + 1
    lines[i] = line
  }

  close(con)
  lines
}






#filepath="~/Riston-D/E/Hankkeet/Luke/Skenaariot/2018/Vertailutasolaskenta-2018/Karikkeet-KHKinventaarista/YASSO_TST.dat"
#filepath="~/Riston-D/E/Hankkeet/Luke/Skenaariot/2018/Vertailutasolaskenta-2018/Karikkeet-KHKinventaarista/YASSO_TST.txt"
filepath = MELAInputFile

n_lines = countLines(filepath)
if(InputFileDebug) {
    print(paste("The file has ", as.character(n_lines+1), " lines.", sep = ""))
}

lines <- 1:n_lines
lines <- readAllLines(filepath, lines)

findFirstAndSecondLine <- function(lines, n_lines) {
	nf = 0
	ns = 0
	first = FALSE
	for(i in 1:n_lines) {
		if(!first) {
			if(substr(lines[i],1,7) == "Kankaat") {
				nf = i
				first = TRUE
				}
			} else {
				if(substr(lines[i],1,7) == "Kankaat") {
				ns = i
				break;
				}
			}
		}
		c(nf, ns)
}


n_first = 0
n_second = 0

alku <- findFirstAndSecondLine(lines,n_lines)
if(InputFileDebug) {
print("Line numbers of beginning of Area information")
print(alku)
}



## MelaYears

MYind <- gregexpr(";",lines[1])[[1]]
nYears <- length(MYind)
MelaYears <- 1:nYears
for(i in 1:nYears) {
	if(i == nYears) 
		last = MYind[nYears] + 20
	else
		last = MYind[i+1] - 1
	MelaYears[i] = as.integer(substr(lines[1],MYind[i]+1,last))
}

if(InputFileDebug) {
    print("MelaYears")
    print(MelaYears)
}

## AREAS

l <- gsub(",",".",lines[alku[1]])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
sfUpland <- as.numeric(substr(l,p[1]+1,p[2]-1))
l <- gsub(",",".",lines[alku[1]+2])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
sfUndrained <- as.numeric(substr(l,p[1]+1,p[2]-1))
l <- gsub(",",".",lines[alku[1]+3])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
sfDrained <- as.numeric(substr(l,p[1]+1,p[2]-1))
#
l <- gsub(",",".",lines[alku[2]])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
nfUpland <- as.numeric(substr(l,p[1]+1,p[2]-1))
l <- gsub(",",".",lines[alku[2]+2])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
nfUndrained <- as.numeric(substr(l,p[1]+1,p[2]-1))
l <- gsub(",",".",lines[alku[2]+3])    #replace , -> . and retu
p <- gregexpr(";",l)[[1]]
nfDrained <- as.numeric(substr(l,p[1]+1,p[2]-1))

Areas <- data.frame(rbind(c(1,1,sfUpland*1000),c(1,2,sfDrained*1000),c(1,3,1000*sfUndrained),
                     c(2,1,nfUpland*1000),c(2,2,nfDrained*1000),c(2,3,nfUndrained*1000)))           
colnames(Areas)<-c("region","soil","area")

if(InputFileDebug) {
        print("AREAS")
        print(Areas)
}


#nYears <- length(MelaYears)


## STOCK
# There are 9 tree - site combinations for Southern and Northern Finland

cNames <- c("region", "soil", "species")
for(i in 1:nYears) {
	cNames <- c(cNames, paste("X",as.character(MelaYears[i]),sep=""))
}

st <- data.frame(t(rep(0,times=length(MelaYears)+3))) 
colnames(st) <- cNames
st <- st[-1,]       #empty

stockRow = 6

#Southern Finland == alku[1]

for(i in 0:8) {
	stock <- 1:nYears
	l <- gsub(",",".",lines[alku[1]+stockRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:nYears) {
		if( j == nYears)  last = ind[nYears] + 20 else last = ind[j+1] - 1
	stock[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(st[,1])
	st[lo+1,] = c(1,soil, species, stock)	
}

#Northern Finland == alku[2]

for(i in 0:8) {
	stock <- 1:nYears
	l <- gsub(",",".",lines[alku[2]+stockRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:nYears) {
		if( j == nYears)  last = ind[nYears] + 20 else last = ind[j+1] - 1
	stock[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(st[,1])
	st[lo+1,] = c(2,soil, species, stock)	
}

if(InputFileDebug) {
     print("First stock line ")
     print(st[1,])
}


## Stem volume of natural mortality (1000 m3 / year)
# There are 9 tree - site combinations for Southern and Northern Finland
# Mortality is given as m3/year between two points of time. There are thus one less
# values of mortality than stock 

cNames <- c("region", "soil", "species")
for(i in 1:(nYears-1)) {
	cNames <- c(cNames, paste("X",as.character(MelaYears[i]),sep=""))
}

mort <- data.frame(t(rep(0,times=length(MelaYears)-1+3))) 
colnames(mort) <- cNames
mort <- mort[-1,]       #empty

mortRow = 17

#Southern Finland == alku[1]

for(i in 0:8) {
	mstock <- 1:(nYears-1)
	l <- gsub(",",".",lines[alku[1]+mortRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:(nYears-1)) {
		if( j == (nYears-1))  last = ind[nYears-1] + 20 else last = ind[j+1] - 1
	mstock[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(mort[,1])
	mort[lo+1,] = c(1,soil, species, mstock)	
}

#Northern Finland == alku[2]

for(i in 0:8) {
	mstock <- 1:(nYears-1)
	l <- gsub(",",".",lines[alku[2]+mortRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:(nYears-1)) {
		if( j == (nYears-1))  last = ind[nYears-1] + 20 else last = ind[j+1] - 1
	mstock[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(mort[,1])
	mort[lo+1,] = c(2,soil, species, mstock)	
}

if(InputFileDebug) {
     print("First Natural mortality line ")
     print(mort[1,])
}


## Cutting waste (1000 t DM / year)
# There are 45 	soil, species, part combinations for Southern and Northern Finland
# Mortality is given as t DM / year between two points of time. There are thus one less
# values of mortality than stock 

cNames <- c("region", "soil", "species", "part")
for(i in 1:(nYears-1)) {
	cNames <- c(cNames, paste("X",as.character(MelaYears[i]),sep=""))
}

waste <- data.frame(t(rep(0,times=length(MelaYears)-1+4))) 
colnames(waste) <- cNames
waste <- waste[-1,]       #empty

wasteRow = 28

#Southern Finland == alku[1]

for(i in 0:44) {
	wb <- 1:(nYears-1)
	l <- gsub(",",".",lines[alku[1]+wasteRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	part <- as.integer(substr(l,5,5))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:(nYears-1)) {
		if( j == (nYears-1))  last = ind[nYears-1] + 20 else last = ind[j+1] - 1
	wb[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(waste[,1])
	waste[lo+1,] = c(1,soil, species, part, wb)	
}

#Northern Finland == alku[2]

for(i in 0:44) {
	wb <- 1:(nYears-1)
	l <- gsub(",",".",lines[alku[2]+wasteRow+i])    #replace , -> . and retu
	soil <- as.integer(substr(l,1,1))
	species <- as.integer(substr(l,3,3))
	part <- as.integer(substr(l,5,5))
	ind <- gregexpr(";",l)[[1]]
	for(j in 1:(nYears-1)) {
		if( j == (nYears-1))  last = ind[nYears-1] + 20 else last = ind[j+1] - 1
	wb[j] = as.numeric(substr(l,ind[j]+1,last))			
	}
	lo = length(waste[,1])
	waste[lo+1,] = c(2,soil, species, part, wb)	
}

if(InputFileDebug) {
 print("First harvest waste line ")
              print(waste[1,])
}




MelaStock <- st
MelaNatMortality <- mort
MelaWaste <- waste

