README Segment files with ditch spacing, relevant fields in csv
 : FID index
Unnamed: 0
gridcode: segment id (not unique among different segments files)
ptype: 1=mineral 4=drained peatland 7=peatland
LENGTH: length of ditches within segment (calculated with a buffer)
m_ha: length of ditches per hectare
part_count: number of parts in segment
Shape_Length: circumference in meters of the segment
Shape_Area: area in m2 of the segment 
m_hacorr: adjusted ditch length per hectare (assigned for segments under 6000m2)
DS: ditch spacing in meters
DSlim: ditch spacing with limits (under 20m DS values replaced with average DS of 36m, over 60m DS values replaced with value 60m)
DSlim: ditch spacing with limits (under 20m DS values replaced with 20m, over 60m DS values replaced with value 60m)

