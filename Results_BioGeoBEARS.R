library(ape)
library(cladoRcpp)
library(rgdal)
library(raster)
library(rgeos)
library(dismo)
library(letsR)
library(BioGeoBEARS)


####################################
# Probabilities of states/ranges at each node
####################################

setwd("C:/Users/carol/Desktop/Biogeo_carol")
res<-readRDS("resBAYAREALIKEj.rds")
#teste<-as.matrix(res$ML_marginal_prob_each_state_at_branch_top_AT_node)

# In this table:
"trfn" = "sigmodontinae_tree"
extdata_dir = np(file.path("extdata", package="BioGeoBEARS"))
getwd()

trfn = np("sigmotree2.nwk")
moref(trfn)
tr = read.tree(trfn)
tr

#  You can see the node numbers in the same APE order with:
trtable = prt(tr, printflag=FALSE)
head(trtable)
tail(trtable)

# You can plot APE node labels with:
plot(tr)
axisPhylo()
nodelabels()
tiplabels(1:length(tr$tip.label))

geogfn = np("MATRIXGEO-3areas.txt") ###ou aqui são os poligonos?
# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=geogfn)
tipranges
# Get your states list (assuming, say, 4-area analysis, with max. rangesize=4)
max_range_size = 3
areas = getareas_from_tipranges_object(tipranges)

# This is the list of states/ranges, where each state/range
# is a list of areas, counting from 0
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=TRUE)

# Make the list of ranges
ranges_list = NULL
for (i in 1:length(states_list_0based)){    
  if ( (length(states_list_0based[[i]]) == 1) && (is.na(states_list_0based[[i]])) )
  {
    tmprange = "_"
  } else {
    tmprange = paste(areas[states_list_0based[[i]]+1], collapse="")
  }
  ranges_list = c(ranges_list, tmprange)
}

# Look at the ranges list
ranges_list

# Make the node numbers the row names
# Make the range_list the column names
range_probabilities = as.data.frame(res$ML_marginal_prob_each_state_at_branch_top_AT_node)
row.names(range_probabilities) = trtable$node
names(range_probabilities) = ranges_list[1:232]
length(ranges_list)
# Look at the table (first six rows)
head(range_probabilities)


# Write the table to a tab-delimited text file (for Excel etc.)
write.table(range_probabilities, file="range_probabilities.txt", quote=FALSE, sep="/t")

# Look at the file
moref("range_probabilities.txt")

rages_absolutos<-range_probabilities[261:519,2:232]
colunas_max <- apply(rages_absolutos, 1, which.max)

#resultados<-c(rages_absolutos,colunas_max)

## Get areas names
col_numbers <- seq_along(range_probabilities)
col_names <- colnames(range_probabilities)
col_info_matrix <- cbind(col_numbers, col_names)
colnames(col_info_matrix) <- c("Column_Number", "Column_Name")
col_info_matrix<-as.data.frame(col_info_matrix)
#join higher probability area names for each node 

indices <- match(colunas_max, col_info_matrix$Column_Number)

# Criar a nova tabela com a lista de números e os nomes correspondentes
new_table <- data.frame(
  col_numbers = colunas_max,
  col_names = col_info_matrix$Column_Name[indices],
  stringsAsFactors = FALSE
)

new_table$nohs <- rownames(new_table)

### pegar shape com as áreas de 1 a 11 que equivalem a A a K
shape.bior <- readOGR(choose.files()) # Folder: Bioregions Sigmodontinae 1x1 - 10 trials
#plot(shape.bior)

# rasterize (turn shape into raster)
r<-raster(ncol=180,nrow=180)
extent(r)<-extent(shape.bior)
raster.bior<-rasterize(shape.bior,r,"BIOREGIO_N")
plot(raster.bior)
str(raster.bior)

celulaeregiao<-shape.bior@data$BIOREGIO_N
length(celulaeregiao)

final_table <- matrix(0, nrow = length(celulaeregiao), ncol = length(261:519))
colnames(final_table) <- as.character(261:519)
rownames(final_table) <- 1:length(celulaeregiao)
final_table<-cbind(celulaeregiao,final_table)

#tabela com letras A:K e valores 1:11
area_values <- data.frame(
  area = LETTERS[1:11],
  value = 1:11,
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(new_table))) {
  node <- new_table$nohs[i]
  areas <- unlist(strsplit(new_table[i,2], ""))
  area_values_subset <- area_values$value[area_values$area %in% areas]
  
  for (value in area_values_subset) {
    final_table[, as.character(node)][celulaeregiao == value] <- 1
  }
}

#saveRDS(final_table,file="novo_noh_por_celula.rds")
#write.csv(new_table,"no_area.csv", row.names = F)

############################################################

final_table<-as.data.frame(final_table)
number_to_letter <- setNames(LETTERS[1:11], 1:11)
final_table$letters <- number_to_letter[as.character(final_table$celulaeregiao)]
final_table

#write.csv(final_table,"C:/Users/carol/Documents/Workshop_meninas/biogeobears-meninas/resultados/final_table.csv", row.names = F)

####### novos grids!!!

r<-raster(ncol=60,nrow=60)
extent(r)<-extent(shape.bior)
raster.bior<-rasterize(shape.bior,r,"BIOREGIO_N")
plot(raster.bior)

values_bioregio <- getValues(raster.bior)
grid <- as(raster.bior, "SpatialGridDataFrame")
centroids <- coordinates(grid)
str(grid)
results <- data.frame(
  BIOREGIO_N = values_bioregio,
  longitude = centroids[, 1],
  latitude = centroids[, 2],
  cel_number = 1:length(values_bioregio)
)
results <- na.omit(results)
nrow(results)

#raster.bior_10<-raster.bior ## 6x6
#raster.bior_20<-raster.bior ## 9x9
#raster.bior_50<-raster.bior ## 13x13
#raster.bior_100<-raster.bior ## 19x19
#raster.bior_500<-raster.bior ## 43x43
#raster.bior_1000<-raster.bior ## 60x60
#gridss<-c(raster.bior_10,raster.bior_20,raster.bior_50,raster.bior_100,raster.bior_500,raster.bior_1000)
#saveRDS(gridss, file="C:/Users/carol/Documents/Workshop_meninas/biogeobears-meninas/resultados/raster_grid.rds")

gridsss<-readRDS("raster_grid.rds")

### Plot grids
par(mfrow=c(2, 3))
my_colors <- colorRampPalette(c("yellow","pink", "green", "blue", "red"))
titles<-c("grid_10","grid_20","grid_50","grid_100","grid_500","grid_1000")
letter_to_index <- setNames(1:12, LETTERS[1:12])
rev_names <- rev(names(letter_to_index))
rev_colors <- my_colors(12)[rev(letter_to_index)]
# Plotar cada raster na lista
for (r in 1:length(gridsss)) {
  plot(gridsss[[r]], col = my_colors(12),main=titles[r])
  legend("bottomleft", legend = rev_names, fill = rev_colors,cex = 0.58,y.intersp = 0.6,bty="n")
}

