#### Prepare Input Data ####
#' read in dimensions
#'
#' @description: reads dimensions of ERA5 data
#' @param filename name of file to read in
#' @return no return 
#' @export
#' @importFrom ncdf4 nc_open ncvar_get 
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data_dims = readin_dim(myfile)
readin_dim=function(filename) {
	dat=ncdf4::nc_open(filename)
	lon=ncdf4::ncvar_get(dat,'longitude')
	lat=ncdf4::ncvar_get(dat,'latitude')
	lev=ncdf4::ncvar_get(dat,'level')*100
	nc_close(dat)
	return(list('lon'=lon,'lat'=lat,'lev'=lev))
}

#' read in ERA5 data
#'
#' @description: reads ERA5 data
#' @param filename name of file to read in
#' @return no return
#' @export
#' @importFrom ncdf4 nc_open ncvar_get nc_close
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
readin_era5=function(filename) {
	dat=ncdf4::nc_open(filename)
	lon=ncdf4::ncvar_get(dat,'longitude')
	lat=ncdf4::ncvar_get(dat,'latitude')
	lev=ncdf4::ncvar_get(dat,'level')*100
	z=ncdf4::ncvar_get(dat,'z')
	temp=ncdf4::ncvar_get(dat,'t')
	u=ncdf4::ncvar_get(dat,'u')
	v=ncdf4::ncvar_get(dat,'v')
	w=ncdf4::ncvar_get(dat,'w')
	vo=ncdf4::ncvar_get(dat,'vo')
	pv=ncdf4::ncvar_get(dat,'pv')
	nc_close(dat)
	return(list('lon'=lon,'lat'=lat,'lev'=lev,'z'=z,'temp'=temp,'u'=u,'v'=v,'w'=w,'vo'=vo,'pv'=pv))
}
