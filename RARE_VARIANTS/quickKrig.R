
###############################
## R script for kriging
###############################
##Required packages: "fields"
##Based on function tessplot.R
##with addition and changes
###############################
require(fields)
require(raster)

# color palette easy to change
# (100) refers to the number of colors in the palette. 
colrp = colorRampPalette(c("blue","cyan","orange","red","darkred"))(100)

# Define a function to create a grid from an ascii raster
createGridFromAsciiRaster = function(file)
  
  # return a grid to use to project estimated data on a map
  # the grid is computed from an ascii raster file
{  

  info=read.table(file,nrows=6)
  grid.info=info[,2]
  names(grid.info)=info[,1]
  if(is.na(grid.info["XLLCORNER"])==TRUE){
  lat.pix=seq(from=grid.info["YLLCENTER"],by=grid.info["CELLSIZE"],length=grid.info["NROWS"])
  long.pix=seq(from=grid.info["XLLCENTER"],by=grid.info["CELLSIZE"],length=grid.info["NCOLS"])
  } else {
  lat.pix=seq(from=grid.info["YLLCORNER"]-grid.info["CELLSIZE"]/2,by=grid.info["CELLSIZE"],length=grid.info["NROWS"])
  long.pix=seq(from=grid.info["XLLCORNER"]-grid.info["CELLSIZE"]/2,by=grid.info["CELLSIZE"],length=grid.info["NCOLS"])
  }
  
  grid=make.surface.grid( list( long.pix,lat.pix))
  
  return(grid)
}

getConstraintsFromAsciiRaster = function(file,cell_value_min=NULL,cell_value_max=NULL)
  
  # return a matrix with boolean cells
  # cell is set to TRUE if user wants to project on this map's cell
  # cell is set to FALSE otherwise
  # if needed user is invited to add/remove constraints by modifying the function
  
{
  map=read.table(file,skip=6)
  map=t(map)[,nrow(map):1]
  
  # 3 suggested constraints that you can remove if needed:
  
  map_constraints=!is.na(map) 	
  
  if (!is.null(cell_value_min)) {
    map_constraints[map<cell_value_min]=FALSE
  }
  
  if (!is.null(cell_value_max)) {
    map_constraints[map>cell_value_max]=FALSE
  }
  
  # you can add constraints here, example:
  # map_constraints["condition where you do not want to project clusters"]=FALSE	
  
  return(map_constraints)
}


# Define tessplot function
quickKrig = function( 
  coordinates = mydata.coord, 
  data_to_Krig = outfile.txt,
  datacolumn = 1,
  asciifile = NULL,
  cell_value_min = NULL,
  cell_value_max = NULL,
  colpalette = colrp,
  mapadd=T,
  pts.size = .4,
  pts.shape = 19,
  theta=10,
  ...
){
  
  
  # read spatial coordinates from the tess input file  
  coordinates = coordinates
  
  # read admixture (tess or clummp output)
  cluster = data_to_Krig[, datacolumn]

  # Create a grid on which the data will be evaluated
  if(is.null(asciifile)==FALSE){
  grid <- createGridFromAsciiRaster(asciifile)
  } else {grid <- NULL}
  
  # Get constraints for prediction
  constraints <- NULL
  if(is.null(cell_value_min)==FALSE |is.null(cell_value_max)==FALSE ){
  constraints <- getConstraintsFromAsciiRaster(asciifile,cell_value_min,cell_value_max)
  }
  
  # Make Krigging
  fit = Krig(coordinates,cluster,m = 1,theta = theta)

  if(is.null(grid)){
  surface(fit, col = colpalette, levels = c(-1), extrap = T,...)
  points(coordinates, cex = pts.size, pch = pts.shape)
  if(mapadd) {map(add=T, interior = F, lwd = 1)}
  } else {
    look<- predict(fit,grid) # evaluate on a grid of points
    out<- as.surface( grid, look)
    
    if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
  ncolors=length(colpalette)
  image.plot(out,col=colpalette,...)
  points(coordinates,pch=pts.shape, cex=pts.size)
  }
}

# Define tessplot function
quickKrig2 = function( 
  coordinates = mydata.coord, 
  data_to_Krig = outfile.txt,
  datacolumn = 1,
  asciifile = NULL,
  cell_value_min = NULL,
  cell_value_max = NULL,
  colpalette = colrp,
  mapadd=T,
  pts.size = .4,
  pts.shape = 19,
  theta=10,
  ...
){
  
  
  # read spatial coordinates from the tess input file  
  coordinates = coordinates
  
  # read admixture (tess or clummp output)
  cluster = data_to_Krig[, datacolumn]
  
  # Create a grid on which the data will be evaluated
  if(is.null(asciifile)==FALSE){
    grid <- createGridFromAsciiRaster(asciifile)
  } else {grid <- NULL}
  
  # Get constraints for prediction
  constraints <- NULL
  if(is.null(cell_value_min)==FALSE |is.null(cell_value_max)==FALSE ){
    constraints <- getConstraintsFromAsciiRaster(asciifile,cell_value_min,cell_value_max)
  }
  
  # Make Krigging
  fit = Krig(coordinates,cluster,m = 1,theta = theta)
  
  if(is.null(grid)){
    surface(fit, col = colpalette, levels = c(-1), extrap = T,...)
    points(coordinates, cex = pts.size, pch = pts.shape)
    if(mapadd) {map(add=T, interior = F, lwd = 1)}
  } else {
    look<- predict(fit,grid) # evaluate on a grid of points
    out<- as.surface( grid, look)
    
    if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
    ncolors=length(colpalette)
    plot(raster(out),col=colpalette,...)
    points(coordinates,pch=pts.shape, cex=pts.size)
  }
}

# Define tessplot function
quickKrig3 = function( 
  coordinates = mydata.coord, 
  data_to_Krig = outfile.txt,
  datacolumn = 1,
  asciifile = NULL,
  cell_value_min = NULL,
  cell_value_max = NULL,
  colpalette = colrp,
  mapadd=T,
  pts.size = .4,
  pts.shape = 19,
  theta=10,
  ...
){
  
  
  # read spatial coordinates from the tess input file  
  coordinates = coordinates
  
  # read admixture (tess or clummp output)
  cluster = data_to_Krig[, datacolumn]
  
  # Create a grid on which the data will be evaluated
  if(is.null(asciifile)==FALSE){
    grid <- createGridFromAsciiRaster(asciifile)
  } else {grid <- NULL}
  
  # Get constraints for prediction
  constraints <- NULL
  if(is.null(cell_value_min)==FALSE |is.null(cell_value_max)==FALSE ){
    constraints <- getConstraintsFromAsciiRaster(asciifile,cell_value_min,cell_value_max)
  }
  
  # Make Krigging
  fit = Krig(coordinates,cluster,m = 1,theta = theta)
  
  if(is.null(grid)){
    surface(fit, col = colpalette, levels = c(-1), extrap = T,...)
    points(coordinates, cex = pts.size, pch = pts.shape)
    if(mapadd) {map(add=T, interior = F, lwd = 1)}
  } else {
    look<- predict(fit,grid) # evaluate on a grid of points
    out<- as.surface( grid, look)
    
    if (class(constraints)!= "NULL") { out[[8]][ !constraints ] = NA }
    ncolors=length(colpalette)
    image(raster(out),col=colpalette,...)
    points(coordinates,pch=pts.shape, cex=pts.size)
  }
}