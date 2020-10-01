### to estimate litter of trees from MELA projections

### Projection data 

#in 1000 m3 stemwood, breakdown by region,soil,species
#region   1 = south,    2 = north
#soil     1 = mineral   2 = organic
#species  1 = pine      2 = spruce    3 = birch & other decidious

#####stock <-  read.csv(paste(path.data,stock.input, sep=""), header=TRUE)
#stock <-  read.csv(MelaStock, header=TRUE)

stock = MelaStock

##! Dividing by area --> m3/ha

mela.area <- Areas

    
v.sf.min.sp <- stock[stock$region==1&stock$soil==1&stock$species==1,4:
                     length(stock[stock$region==1&stock$soil==1&stock$species==1,])]*1000/
                     mela.area$area[mela.area$region==1&mela.area$soil==1] 
v.sf.min.ns <- stock[stock$region==1&stock$soil==1&stock$species==2,4:
                     length(stock[stock$region==1&stock$soil==1&stock$species==2,])]*1000/
                 mela.area$area[mela.area$region==1&mela.area$soil==1]
v.sf.min.dc <- stock[stock$region==1&stock$soil==1&stock$species==3,4:
                     length(stock[stock$region==1&stock$soil==1&stock$species==3,])]*1000/
  mela.area$area[mela.area$region==1&mela.area$soil==1]

v.nf.min.sp <- stock[stock$region==2&stock$soil==1&stock$species==1,4:
                     length(stock[stock$region==2&stock$soil==1&stock$species==1,])]*1000/
             mela.area$area[mela.area$region==2&mela.area$soil==1] 
v.nf.min.ns <- stock[stock$region==2&stock$soil==1&stock$species==2,4:
                     length(stock[stock$region==2&stock$soil==1&stock$species==2,])]*1000/
               mela.area$area[mela.area$region==2&mela.area$soil==1] 
v.nf.min.dc <- stock[stock$region==2&stock$soil==1&stock$species==3,4:
                     length(stock[stock$region==2&stock$soil==1&stock$species==3,])]*1000/
                     mela.area$area[mela.area$region==2&mela.area$soil==1] 


v.sf.org.sp <- stock[stock$region==1&stock$soil==2&stock$species==1,4:
                     length(stock[stock$region==1&stock$soil==2&stock$species==1,])]*1000/
                     mela.area$area[mela.area$region==1&mela.area$soil==2]
v.sf.org.ns <- stock[stock$region==1&stock$soil==2&stock$species==2,4:
                     length(stock[stock$region==1&stock$soil==2&stock$species==2,])]*1000/
                     mela.area$area[mela.area$region==1&mela.area$soil==2]
v.sf.org.dc <- stock[stock$region==1&stock$soil==2&stock$species==3,4:
                     length(stock[stock$region==1&stock$soil==2&stock$species==3,])]*1000/
                     mela.area$area[mela.area$region==1&mela.area$soil==2]

v.nf.org.sp <- stock[stock$region==2&stock$soil==2&stock$species==1,4:
                     length(stock[stock$region==2&stock$soil==2&stock$species==1,])]*1000/
                     mela.area$area[mela.area$region==2&mela.area$soil==2]
v.nf.org.ns <- stock[stock$region==2&stock$soil==2&stock$species==2,4:
                     length(stock[stock$region==2&stock$soil==2&stock$species==2,])]*1000/
                     mela.area$area[mela.area$region==2&mela.area$soil==2]
v.nf.org.dc <- stock[stock$region==2&stock$soil==2&stock$species==3,4:
                     length(stock[stock$region==2&stock$soil==2&stock$species==3,])]*1000/
                     mela.area$area[mela.area$region==2&mela.area$soil==2]

#############
### 1. biomass then litter for SOUTH
#############

##! Lf = litter foliage, Lb = litter branches, Ls = litter stem bark,
##! Lst = litter stump, Lr = litter roots, Lfr = litter fine roots

## Unit t C / ha

## Obs: litter is divided in AWEN -compartment
## Litter (e.g. Lf.sp) is a matrix, columns are A W E N and
## rows are the years on MELA input- note transpose t()

#foliage - mineral soil
Mf.sf.min.sp <- (v.sf.min.sp)* BEF.sp[1,6]
Mf.sf.min.ns <- (v.sf.min.ns)* BEF.ns[1,6]
Mf.sf.min.dc <- (v.sf.min.dc)* BEF.dc[1,6]

Lf.sp <- foliage.AWEN(t(foliage.litter(Mf.sf.min.sp,1,1,1)),1)
Lf.ns <- foliage.AWEN(t(foliage.litter(Mf.sf.min.ns,2,1,1)),2)
Lf.dc <- foliage.AWEN(t(foliage.litter(Mf.sf.min.dc,3,1,1)),3)

Lf.sf.min <- carbon(Lf.sp+Lf.ns+Lf.dc) # here carbon

# foliage - drained org

Mf.sf.org.sp <- (v.sf.org.sp)* BEF.sp[3,6]
Mf.sf.org.ns <- (v.sf.org.ns)* BEF.ns[3,6]
Mf.sf.org.dc <- (v.sf.org.dc)* BEF.dc[3,6]

Lf.sp <- foliage.AWEN(t(foliage.litter(Mf.sf.org.sp,1,1,0)),1)
Lf.ns <- foliage.AWEN(t(foliage.litter(Mf.sf.org.ns,2,1,0)),2)
Lf.dc <- foliage.AWEN(t(foliage.litter(Mf.sf.org.dc,3,1,0)),3)

Lf.sf.org <- carbon(Lf.sp+Lf.ns+Lf.dc) # here carbon 

# branches - mineral soil

Mb.sf.min.sp <- (v.sf.min.sp)* BEF.sp[1,5]
Mb.sf.min.ns <- (v.sf.min.ns)* BEF.ns[1,5]
Mb.sf.min.dc <- (v.sf.min.dc)* BEF.dc[1,5]

Lb.sp <- branches.AWEN(t(branch.litter(Mb.sf.min.sp,1)))/sum(branches.AWEN(1))
Lb.ns <- branches.AWEN(t(branch.litter(Mb.sf.min.ns,2)))/sum(branches.AWEN(1))
Lb.dc <- branches.AWEN(t(branch.litter(Mb.sf.min.dc,3)))/sum(branches.AWEN(1))

Lb.sf.min <- carbon(Lb.sp+Lb.ns+Lb.dc) # here carbon

# branches - drained org

Mb.sf.org.sp <- (v.sf.org.sp)* BEF.sp[3,5]
Mb.sf.org.ns <- (v.sf.org.ns)* BEF.ns[3,5]
Mb.sf.org.dc <- (v.sf.org.dc)* BEF.dc[3,5]

Lb.sp <- branches.AWEN(t(branch.litter(Mb.sf.org.sp,1)))/sum(branches.AWEN(1))
Lb.ns <- branches.AWEN(t(branch.litter(Mb.sf.org.ns,2)))/sum(branches.AWEN(1))
Lb.dc <- branches.AWEN(t(branch.litter(Mb.sf.org.dc,3)))/sum(branches.AWEN(1))

Lb.sf.org <- carbon(Lb.sp+Lb.ns+Lb.dc) # here carbon

#stem bark - mineral soil

Ms.sf.min.sp <- (v.sf.min.sp)* BEF.sp[1,4]
Ms.sf.min.ns <- (v.sf.min.ns)* BEF.ns[1,4]
Ms.sf.min.dc <- (v.sf.min.dc)* BEF.dc[1,4]

#HUOM! stem.AWEN(1,1) does  sum up to 1: sum(stem.AWEN(1,1)) = 0.985
# neither do other (= 0.985, 0.99)

Ls.sp <- stem.AWEN(t(bark.litter(Ms.sf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Ls.ns <- stem.AWEN(t(bark.litter(Ms.sf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Ls.dc <- stem.AWEN(t(bark.litter(Ms.sf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Ls.sf.min <- carbon(Ls.sp+Ls.ns+Ls.dc) # here carbon

# stem bark - drained org

Ms.sf.org.sp <- (v.sf.org.sp)* BEF.sp[3,4]
Ms.sf.org.ns <- (v.sf.org.ns)* BEF.ns[3,4]
Ms.sf.org.dc <- (v.sf.org.dc)* BEF.dc[3,4]

Ls.sp <- stem.AWEN(t(bark.litter(Ms.sf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Ls.ns <- stem.AWEN(t(bark.litter(Ms.sf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Ls.dc <- stem.AWEN(t(bark.litter(Ms.sf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Ls.sf.org <- carbon(Ls.sp+Ls.ns+Ls.dc) # here carbon

# Stump - mineral soil

Mst.sf.min.sp <- (v.sf.min.sp)* BEF.sp[1,8]
Mst.sf.min.ns <- (v.sf.min.ns)* BEF.ns[1,8]
Mst.sf.min.dc <- (v.sf.min.dc)* BEF.dc[1,8]

Lst.sp <- stem.AWEN(t(bark.litter(Mst.sf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Lst.ns <- stem.AWEN(t(bark.litter(Mst.sf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Lst.dc <- stem.AWEN(t(bark.litter(Mst.sf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Lst.sf.min <- carbon(Lst.sp+Lst.ns+Lst.dc) # here carbon

# Stump - drained org

Mst.sf.org.sp <- (v.sf.org.sp)* BEF.sp[3,8]
Mst.sf.org.ns <- (v.sf.org.ns)* BEF.ns[3,8]
Mst.sf.org.dc <- (v.sf.org.dc)* BEF.dc[3,8]

Lst.sp <- stem.AWEN(t(bark.litter(Mst.sf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Lst.ns <- stem.AWEN(t(bark.litter(Mst.sf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Lst.dc <- stem.AWEN(t(bark.litter(Mst.sf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Lst.sf.org <- carbon(Lst.sp+Lst.ns+Lst.dc) # here carbon

## Roots - mineral soil

Mr.sf.min.sp <- (v.sf.min.sp)* BEF.sp[1,9]
Mr.sf.min.ns <- (v.sf.min.ns)* BEF.ns[1,9]
Mr.sf.min.dc <- (v.sf.min.dc)* BEF.dc[1,9]

Lr.sp <- stem.AWEN(t(root.litter(Mr.sf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Lr.ns <- stem.AWEN(t(root.litter(Mr.sf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Lr.dc <- stem.AWEN(t(root.litter(Mr.sf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Lr.sf.min <- carbon(Lr.sp+Lr.ns+Lr.dc) # here carbon

# Roots - drained org

Mr.sf.org.sp <- (v.sf.org.sp)* BEF.sp[3,9]
Mr.sf.org.ns <- (v.sf.org.ns)* BEF.ns[3,9]
Mr.sf.org.dc <- (v.sf.org.dc)* BEF.dc[3,9]

Lr.sp <- stem.AWEN(t(root.litter(Mr.sf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Lr.ns <- stem.AWEN(t(root.litter(Mr.sf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Lr.dc <- stem.AWEN(t(root.litter(Mr.sf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Lr.sf.org <- carbon(Lr.sp+Lr.ns+Lr.dc) # here carbon

## Fineroots - mineral soil

#fol.rat.sf <-  BEF7.sf[1:4, 9]/ BEF7.sf[1:4,8] 
fol.rat.sf <- FoliageRatioSouth

################NOTE Fol Rat file here!!!

Mfr.sf.min.sp <- fineroots(((v.sf.min.sp)* BEF.sp[1,6]*fol.rat.sf[1]),1)
Mfr.sf.min.ns <- fineroots(((v.sf.min.ns)* BEF.ns[1,6]*fol.rat.sf[3]),2)
Mfr.sf.min.dc <- fineroots(((v.sf.min.dc)* BEF.dc[1,6]),3)

Lfr.sp <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.min.sp)),1)
Lfr.ns <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.min.ns)),2)
Lfr.dc <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.min.dc)),3)

Lfr.sf.min <- carbon(Lfr.sp+Lfr.ns+Lfr.dc) # here carbons

# Fineroots - drained org

Mfr.sf.org.sp <- fineroots(((v.sf.org.sp)* BEF.sp[3,6]*fol.rat.sf[2]),1)
Mfr.sf.org.ns <- fineroots(((v.sf.org.ns)* BEF.ns[3,6]*fol.rat.sf[4]),2)
Mfr.sf.org.dc <- fineroots(((v.sf.org.dc)* BEF.dc[3,6]),3)

Lfr.sp <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.org.sp)),1)
Lfr.ns <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.org.ns)),2)
Lfr.dc <- fineroot.AWEN(t(fineroot.litter(Mfr.sf.org.dc)),3)

Lfr.sf.org <- carbon(Lfr.sp+Lfr.ns+Lfr.dc) # here carbon

####LITTER DATAFRAME
df.litter <- data.frame(year = c(2000,2005,2010,2015,2020,2030,2040,2050,2060),
Mfr.sf.min.sp = unlist(Mfr.sf.min.sp,use.names = FALSE),
Mfr.sf.min.ns= unlist(Mfr.sf.min.ns,use.names = FALSE),
Mfr.sf.min.dc=unlist(Mfr.sf.min.dc,use.names = FALSE),
Mf.sf.min.sp = unlist(Mf.sf.min.sp,use.names = FALSE),
Mf.sf.min.ns= unlist(Mf.sf.min.ns,use.names = FALSE),
Mf.sf.min.dc=unlist(Mf.sf.min.dc,use.names = FALSE),
stringsAsFactors=FALSE) 


#############
### 2. biomass then litter for NORTH
#############

#Foliage - mineral soil
Mf.nf.min.sp <- (v.nf.min.sp)* BEF.sp[4,6]
Mf.nf.min.ns <- (v.nf.min.ns)* BEF.ns[4,6]
Mf.nf.min.dc <- (v.nf.min.dc)* BEF.dc[4,6]

Lf.sp <- foliage.AWEN(t(foliage.litter(Mf.nf.min.sp,1,2,1)),1)
Lf.ns <- foliage.AWEN(t(foliage.litter(Mf.nf.min.ns,2,2,1)),2)
Lf.dc <- foliage.AWEN(t(foliage.litter(Mf.nf.min.dc,3,2,1)),3)

Lf.nf.min <- carbon(Lf.sp+Lf.ns+Lf.dc) # here carbon

# Foliage - drained org

Mf.nf.org.sp <- (v.nf.org.sp)* BEF.sp[6,6]
Mf.nf.org.ns <- (v.nf.org.ns)* BEF.ns[6,6]
Mf.nf.org.dc <- (v.nf.org.dc)* BEF.dc[6,6]

Lf.sp <- foliage.AWEN(t(foliage.litter(Mf.nf.org.sp,1,2,0)),1)
Lf.ns <- foliage.AWEN(t(foliage.litter(Mf.nf.org.ns,2,2,0)),2)
Lf.dc <- foliage.AWEN(t(foliage.litter(Mf.nf.org.dc,3,2,0)),3)

Lf.nf.org <- carbon(Lf.sp+Lf.ns+Lf.dc) # here carbon

# Branch - mineral soil

Mb.nf.min.sp <- (v.nf.min.sp)* BEF.sp[4,5]
Mb.nf.min.ns <- (v.nf.min.ns)* BEF.ns[4,5]
Mb.nf.min.dc <- (v.nf.min.dc)* BEF.dc[4,5]

Lb.sp <- branches.AWEN(t(branch.litter(Mb.nf.min.sp,1)))/sum(branches.AWEN(1))
Lb.ns <- branches.AWEN(t(branch.litter(Mb.nf.min.ns,2)))/sum(branches.AWEN(1))
Lb.dc <- branches.AWEN(t(branch.litter(Mb.nf.min.dc,3)))/sum(branches.AWEN(1))

Lb.nf.min <- carbon(Lb.sp+Lb.ns+Lb.dc) # here carbon

# Branch - drained org

Mb.nf.org.sp <- (v.nf.org.sp)* BEF.sp[6,5]
Mb.nf.org.ns <- (v.nf.org.ns)* BEF.ns[6,5]
Mb.nf.org.dc <- (v.nf.org.dc)* BEF.dc[6,5]

Lb.sp <- branches.AWEN(t(branch.litter(Mb.nf.org.sp,1)))/sum(branches.AWEN(1))
Lb.ns <- branches.AWEN(t(branch.litter(Mb.nf.org.ns,2)))/sum(branches.AWEN(1))
Lb.dc <- branches.AWEN(t(branch.litter(Mb.nf.org.dc,3)))/sum(branches.AWEN(1))

Lb.nf.org <- carbon(Lb.sp+Lb.ns+Lb.dc) # here carbon

# Stem bark - mineral soil

Ms.nf.min.sp <- (v.nf.min.sp)* BEF.sp[4,4]
Ms.nf.min.ns <- (v.nf.min.ns)* BEF.ns[4,4]
Ms.nf.min.dc <- (v.nf.min.dc)* BEF.dc[4,4]

####

Ls.sp <- stem.AWEN(t(bark.litter(Ms.nf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Ls.ns <- stem.AWEN(t(bark.litter(Ms.nf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Ls.dc <- stem.AWEN(t(bark.litter(Ms.nf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Ls.nf.min <- carbon(Ls.sp+Ls.ns+Ls.dc) # here carbon

# Stem bark - drained org

Ms.nf.org.sp <- (v.nf.org.sp)* BEF.sp[6,4]
Ms.nf.org.ns <- (v.nf.org.ns)* BEF.ns[6,4]
Ms.nf.org.dc <- (v.nf.org.dc)* BEF.dc[6,4]

Ls.sp <- stem.AWEN(t(bark.litter(Ms.nf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Ls.ns <- stem.AWEN(t(bark.litter(Ms.nf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Ls.dc <- stem.AWEN(t(bark.litter(Ms.nf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Ls.nf.org <- carbon(Ls.sp+Ls.ns+Ls.dc) # here carbon

# Stump - mineral soil

Mst.nf.min.sp <- (v.nf.min.sp)* BEF.sp[4,8]
Mst.nf.min.ns <- (v.nf.min.ns)* BEF.ns[4,8]
Mst.nf.min.dc <- (v.nf.min.dc)* BEF.dc[4,8]

Lst.sp <- stem.AWEN(t(bark.litter(Mst.nf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Lst.ns <- stem.AWEN(t(bark.litter(Mst.nf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Lst.dc <- stem.AWEN(t(bark.litter(Mst.nf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Lst.nf.min <- carbon(Lst.sp+Lst.ns+Lst.dc) # here carbon

# Stump - drained org

Mst.nf.org.sp <- (v.nf.org.sp)* BEF.sp[6,8]
Mst.nf.org.ns <- (v.nf.org.ns)* BEF.ns[6,8]
Mst.nf.org.dc <- (v.nf.org.dc)* BEF.dc[6,8]

Lst.sp <- stem.AWEN(t(bark.litter(Mst.nf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Lst.ns <- stem.AWEN(t(bark.litter(Mst.nf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Lst.dc <- stem.AWEN(t(bark.litter(Mst.nf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Lst.nf.org <- carbon(Lst.sp+Lst.ns+Lst.dc) # here carbon

## Roots - mineral soil

Mr.nf.min.sp <- (v.nf.min.sp)* BEF.sp[4,9]
Mr.nf.min.ns <- (v.nf.min.ns)* BEF.ns[4,9]
Mr.nf.min.dc <- (v.nf.min.dc)* BEF.dc[4,9]

Lr.sp <- stem.AWEN(t(root.litter(Mr.nf.min.sp,1)),1)/sum(stem.AWEN(1,1))
Lr.ns <- stem.AWEN(t(root.litter(Mr.nf.min.ns,2)),2)/sum(stem.AWEN(1,2))
Lr.dc <- stem.AWEN(t(root.litter(Mr.nf.min.dc,3)),3)/sum(stem.AWEN(1,3))

Lr.nf.min <- carbon(Lr.sp+Lr.ns+Lr.dc) # here carbon

# Roots - drained org

Mr.nf.org.sp <- (v.nf.org.sp)* BEF.sp[6,9]
Mr.nf.org.ns <- (v.nf.org.ns)* BEF.ns[6,9]
Mr.nf.org.dc <- (v.nf.org.dc)* BEF.dc[6,9]

Lr.sp <- stem.AWEN(t(root.litter(Mr.nf.org.sp,1)),1)/sum(stem.AWEN(1,1))
Lr.ns <- stem.AWEN(t(root.litter(Mr.nf.org.ns,2)),2)/sum(stem.AWEN(1,2))
Lr.dc <- stem.AWEN(t(root.litter(Mr.nf.org.dc,3)),3)/sum(stem.AWEN(1,3))

Lr.nf.org <- carbon(Lr.sp+Lr.ns+Lr.dc) # here carbon


## Fineroots - mineral soil

#fol.rat.nf <-  BEF7.nf[1:4, 9]/ BEF7.nf[1:4,8] 
fol.rat.nf <- FoliageRatioNorth

Mfr.nf.min.sp <- fineroots(((v.nf.min.sp)* BEF.sp[4,6]*fol.rat.nf[1]),1)
Mfr.nf.min.ns <- fineroots(((v.nf.min.ns)* BEF.ns[4,6]*fol.rat.nf[3]),2)
Mfr.nf.min.dc <- fineroots(((v.nf.min.dc)* BEF.dc[4,6]),3)

Lfr.sp <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.min.sp)),1)
Lfr.ns <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.min.ns)),2)
Lfr.dc <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.min.dc)),3)

Lfr.nf.min <- carbon(Lfr.sp+Lfr.ns+Lfr.dc) # here carbon

# Fineroots - drained org

Mfr.nf.org.sp <- fineroots(((v.nf.org.sp)* BEF.sp[6,6]*fol.rat.nf[2]),1)
Mfr.nf.org.ns <- fineroots(((v.nf.org.ns)* BEF.ns[6,6]*fol.rat.nf[4]),2)
Mfr.nf.org.dc <- fineroots(((v.nf.org.dc)* BEF.dc[6,6]),3)

Lfr.sp <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.org.sp)),1)
Lfr.ns <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.org.ns)),2)
Lfr.dc <- fineroot.AWEN(t(fineroot.litter(Mfr.nf.org.dc)),3)

Lfr.nf.org <- carbon(Lfr.sp+Lfr.ns+Lfr.dc) # here carbon



##################
### 3. Summing up litter input 
#################

#Mineral soils

# south
nwl.sf.min <- Lf.sf.min + Lfr.sf.min
fwl.sf.min <- Lb.sf.min + Ls.sf.min + Lst.sf.min + Lr.sf.min

# north
nwl.nf.min <- Lf.nf.min + Lfr.nf.min
fwl.nf.min <- Lb.nf.min + Ls.nf.min + Lst.nf.min + Lr.nf.min

# Drained organic soils

# south
nwl.sf.org.abv <- Lf.sf.org 
nwl.sf.org.bel <- Lfr.sf.org

fwl.sf.org.abv <- Lb.sf.org + Ls.sf.org + Lst.sf.org
fwl.sf.org.bel <- Lr.sf.org

# north
nwl.nf.org.abv <- Lf.nf.org 
nwl.nf.org.bel <- Lfr.nf.org

fwl.nf.org.abv <- Lb.nf.org + Ls.nf.org + Lst.nf.org
fwl.nf.org.bel <- Lr.nf.org

########################
###  4. Interpolating time series for Mineral soils 
########################


###    index <- c(2011,2016,2021,2026,2031,2041,2051, 2061)
index <- MelaYears
    
xout<-seq(from=min(index),to=max(index),by=1)

nwl.sf.min.final <- cbind(approx(index,nwl.sf.min[,1],xout,rule=2)$y,approx(index,nwl.sf.min[,2],xout,rule=2)$y,
                         approx(index,nwl.sf.min[,3],xout,rule=2)$y, approx(index,nwl.sf.min[,4],xout,rule=2)$y)

fwl.sf.min.final <- cbind(approx(index,fwl.sf.min[,1],xout,rule=2)$y,approx(index,fwl.sf.min[,2],xout,rule=2)$y,
                         approx(index,fwl.sf.min[,3],xout,rule=2)$y, approx(index,fwl.sf.min[,4],xout,rule=2)$y)
#IF nf == 0:
if(nwl.nf.min[,1][1]  == "NaN" ){
nwl.nf.min[is.na(nwl.nf.min)] <- 0
fwl.nf.min[is.na(fwl.nf.min)] <- 0
}
nwl.nf.min.final <- cbind(approx(index,nwl.nf.min[,1],xout,rule=2)$y,approx(index,nwl.nf.min[,2],xout,rule=2)$y,
                         approx(index,nwl.nf.min[,3],xout,rule=2)$y, approx(index,nwl.nf.min[,4],xout,rule=2)$y)

fwl.nf.min.final <- cbind(approx(index,fwl.nf.min[,1],xout,rule=2)$y,approx(index,fwl.nf.min[,2],xout,rule=2)$y,
                          approx(index,fwl.nf.min[,3],xout,rule=2)$y, approx(index,fwl.nf.min[,4],xout,rule=2)$y)


sf.min <- cbind(xout,nwl.sf.min.final,fwl.sf.min.final)
colnames(sf.min)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")

nf.min <- cbind(xout,nwl.nf.min.final,fwl.nf.min.final)
colnames(nf.min)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")

##!  t C/ha
#These intermdiate files are used in Yasso calculations

write.table(round(sf.min,4), file=MineralSoilSouth,row.names=FALSE )
write.table(round(nf.min,4), file=MineralSoilNorth, row.names=FALSE)
    

########################
###  5. Interpolating time series for drained organic soils 
########################
if(nwl.sf.org.abv[,1][1]  == "NaN" ){
nwl.sf.org.abv[is.na(nwl.sf.org.abv)] <- 0
fwl.sf.org.abv[is.na(fwl.sf.org.abv)] <- 0
}

nwl.sf.org.abv.final <- cbind(approx(index,nwl.sf.org.abv[,1],xout)$y,approx(index,nwl.sf.org.abv[,2],xout)$y,
                              approx(index,nwl.sf.org.abv[,3],xout)$y, approx(index,nwl.sf.org.abv[,4],xout)$y)
fwl.sf.org.abv.final <- cbind(approx(index,fwl.sf.org.abv[,1],xout)$y,approx(index,fwl.sf.org.abv[,2],xout)$y,
                              approx(index,fwl.sf.org.abv[,3],xout)$y, approx(index,fwl.sf.org.abv[,4],xout)$y)
if(nwl.nf.org.abv[,1][1]  == "NaN" ){
nwl.nf.org.abv[is.na(nwl.nf.org.abv)] <- 0
fwl.nf.org.abv[is.na(fwl.nf.org.abv)] <- 0
}

nwl.nf.org.abv.final <- cbind(approx(index,nwl.nf.org.abv[,1],xout)$y,approx(index,nwl.nf.org.abv[,2],xout)$y,
                              approx(index,nwl.nf.org.abv[,3],xout)$y, approx(index,nwl.nf.org.abv[,4],xout)$y)
fwl.nf.org.abv.final <- cbind(approx(index,fwl.nf.org.abv[,1],xout)$y,approx(index,fwl.nf.org.abv[,2],xout)$y,
                              approx(index,fwl.nf.org.abv[,3],xout)$y, approx(index,fwl.nf.org.abv[,4],xout)$y)

if(nwl.sf.org.bel[,1][1]  == "NaN" ){
nwl.sf.org.bel[is.na(nwl.sf.org.bel)] <- 0
fwl.sf.org.bel[is.na(fwl.sf.org.bel)] <- 0
}
nwl.sf.org.bel.final <- cbind(approx(index,nwl.sf.org.bel[,1],xout)$y,approx(index,nwl.sf.org.bel[,2],xout)$y,
                              approx(index,nwl.sf.org.bel[,3],xout)$y, approx(index,nwl.sf.org.bel[,4],xout)$y)
fwl.sf.org.bel.final <- cbind(approx(index,fwl.sf.org.bel[,1],xout)$y,approx(index,fwl.sf.org.bel[,2],xout)$y,
                              approx(index,fwl.sf.org.bel[,3],xout)$y, approx(index,fwl.sf.org.bel[,4],xout)$y)
if(nwl.nf.org.bel[,1][1]  == "NaN" ){
nwl.nf.org.bel[is.na(nwl.nf.org.bel)] <- 0
fwl.nf.org.bel[is.na(fwl.nf.org.bel)] <- 0
}

nwl.nf.org.bel.final <- cbind(approx(index,nwl.nf.org.bel[,1],xout)$y,approx(index,nwl.nf.org.bel[,2],xout)$y,
                              approx(index,nwl.nf.org.bel[,3],xout)$y, approx(index,nwl.nf.org.bel[,4],xout)$y)
fwl.nf.org.bel.final <- cbind(approx(index,fwl.nf.org.bel[,1],xout)$y,approx(index,fwl.nf.org.bel[,2],xout)$y,
                              approx(index,fwl.nf.org.bel[,3],xout)$y, approx(index,fwl.nf.org.bel[,4],xout)$y)

sf.org.abv <- cbind(xout,nwl.sf.org.abv.final,fwl.sf.org.abv.final)
nf.org.abv <- cbind(xout,nwl.nf.org.abv.final,fwl.nf.org.abv.final)

sf.org.bel <- cbind(xout,nwl.sf.org.bel.final,fwl.sf.org.bel.final)
nf.org.bel <- cbind(xout,nwl.nf.org.bel.final,fwl.nf.org.bel.final)

colnames(sf.org.abv)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")
colnames(nf.org.abv)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")
colnames(sf.org.bel)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")
colnames(nf.org.bel)=c("year","nwl_A","nwl_W","nwl_E","nwl_N","fwl_A","fwl_W","fwl_E","fwl_N")

#These intermdiate files are used in organic soils calculations

write.table(round(sf.org.abv,4), file=OrganicSoilAboveSouth, row.names=FALSE )
write.table(round(nf.org.abv,4), file=OrganicSoilAboveNorth, row.names=FALSE )
write.table(round(sf.org.bel,4), file=OrganicSoilBelowSouth ,row.names=FALSE )
write.table(round(nf.org.bel,4), file=OrganicSoilBelowNorth ,row.names=FALSE )


########################
###  5. Litter for output
########################

years <- seq(from = min(index), to = max(index), by = 1)

asfmin = mela.area$area[mela.area$region==1&mela.area$soil==1]
anfmin = mela.area$area[mela.area$region==2&mela.area$soil==1]
asforg = mela.area$area[mela.area$region==1&mela.area$soil==2]
anforg = mela.area$area[mela.area$region==2&mela.area$soil==2]

if(Lf.nf.min[,1][1]  == "NaN" ){
Lf.nf.min[is.na(Lf.nf.min)] <- 0
Lb.nf.min[is.na(Lb.nf.min)] <- 0
Ls.nf.min[is.na(Ls.nf.min)] <- 0
Lst.nf.min[is.na(Lst.nf.min)] <- 0
Lr.nf.min[is.na(Lr.nf.min)] <- 0
Lfr.nf.min[is.na(Lfr.nf.min)] <- 0

Lf.sf.org[is.na(Lf.sf.org)] <- 0
Lb.sf.org[is.na(Lb.sf.org)] <- 0
Ls.sf.org[is.na(Ls.sf.org)] <- 0
Lst.sf.org[is.na(Lst.sf.org)] <- 0
Lr.sf.org[is.na(Lr.sf.org)] <- 0
Lfr.sf.org[is.na(Lfr.sf.org)] <- 0

Lf.nf.org[is.na(Lf.nf.org)] <- 0
Lb.nf.org[is.na(Lb.nf.org)] <- 0
Ls.nf.org[is.na(Ls.nf.org)] <- 0
Lst.nf.org[is.na(Lst.nf.org)] <- 0
Lr.nf.org[is.na(Lr.nf.org)] <- 0
Lfr.nf.org[is.na(Lfr.nf.org)] <- 0
}

litter.mineral<-cbind(years,approx(index,rowSums(asfmin*Lf.sf.min+anfmin*Lf.nf.min),years)$y,
                            approx(index,rowSums(asfmin*Lb.sf.min+anfmin*Lb.nf.min),years)$y,
                            approx(index,rowSums(asfmin*Ls.sf.min+anfmin*Ls.nf.min),years)$y,
                            approx(index,rowSums(asfmin*Lst.sf.min+anfmin*Lst.nf.min),years)$y,
                            approx(index,rowSums(asfmin*Lr.sf.min+anfmin*Lr.nf.min),years)$y,
                            approx(index,rowSums(asfmin*Lfr.sf.min+anfmin*Lfr.nf.min),years)$y)

litter.organic<-cbind(approx(index,rowSums(asforg*Lf.sf.org+anforg*Lf.nf.org),years)$y,
                            approx(index,rowSums(asforg*Lb.sf.org+anforg*Lb.nf.org),years)$y,
                            approx(index,rowSums(asforg*Ls.sf.org+anforg*Ls.nf.org),years)$y,
                            approx(index,rowSums(asforg*Lst.sf.org+anforg*Lst.nf.org),years)$y,
                            approx(index,rowSums(asforg*Lr.sf.org+anforg*Lr.nf.org),years)$y,
                            approx(index,rowSums(asforg*Lfr.sf.org+anforg*Lfr.nf.org),years)$y)

litt <- data.frame(cbind(litter.mineral,litter.organic))
colnames(litt) <- c("year","mlitfol","mlitbra","mlitstem","mlitstump","mlitroots","mlitfr",
                          "olitfol","olitbra","olitstem","olitstump","olitroots","olitfr")

#unit = t C 
#mlit = mineral soil litter
#olit = organic soil litter
# fol = foliage, bra = branches, stem = stem, stump = stump, roots = coarse roots (> 2m mm)
# fr = fine roots

