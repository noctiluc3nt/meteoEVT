#' Thermic Front Parameter (TFP)
#'
#' @description Calculates the thermic front parameter based on the potential temperature
#' @param t_fld temperature field [K]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat only for lonlat mode: vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode the horizontal coordinate system, options are 'lonlat' for a longitude-latitude-grid (default), or 'cartesian' for an equidistant cartesian grid
#'
#' @return thermic front parameter [K/m^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' tfp=calc_tfp(data$temp,data$lev,data$lat)
calc_tfp=function(t_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
	th_fld=calc_theta(t_fld,lev_p)
	thgrad=grad(th_fld,lat,d=2,system='p',rho=NULL,dx,dy,lev_p,mode) #2d gradient
	thgradh_abs=sqrt(thgrad[,,,1]^2+thgrad[,,,2]^2) #absolute value horizontal gradient
	thgrad2=grad(thgradh_abs,lat,d=2,system='p',rho=NULL,dx,dy,lev_p,mode)
	tfp_fld=-(thgrad2[,,,1]*thgrad[,,,1]+thgrad2[,,,2]*thgrad[,,,2])/thgradh_abs
	return(tfp_fld)
}


#' F diagnostic
#'
#' @description Calculates the F diagnostic
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat only for lonlat mode: vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode the horizontal coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#'
#' @return F diagnostic (dimensionless)
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' fdiag=calc_fdiag(data$temp,data$u,data$v,data$w,data$lev,data$lat)
calc_fdiag=function(t_fld,u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
	tgrad=grad(t_fld,lat,d=2,system='p',rho=NULL,dx,dy,lev_p,mode) #2d gradient
	tgrad0=0.45/(100*1000) #constant reference value
	rvort=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=TRUE,zvort_fld=NULL,mode)
	avort=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=FALSE,zvort_fld=NULL,mode)
	coriolis=avort-rvort
	fdiag_fld=rvort[,,,3]*sqrt(tgrad[,,,1]^2+tgrad[,,,2]^2)/(coriolis[,,,3]*tgrad0)
	return(fdiag_fld)
}


#' Front Identification und Statistics
#'
#' @description Calculates frontal zones based on a chosen method (TFP, F diagnostic, DSI) and provides statistics of the distribution of meteorological quantities inside the determined frontak zones.
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param phi_fld geopotential height [gpm]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat only for lonlat mode: vector containing latitude
#' @param method character containing the method, use \code{'tfp'} for TFP method, \code{'f'} for F diagnostic and \code{'dsi'} for DSI method
#' @param threshold scalar containing a suitable threshold (e.g., 2*10^-10 for TFP method, 1 or for F diagnostic, 10^-16 for DSI method)
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param fronts_only if you only want to calculate the frontal regions and not their properties (default FALSE)
#' @param mode the horizontal coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#'
#' @return list containing the used method and used threshold, field with logicals containing the detected frontal zones and numerics of temperature, u-wind, v-wind, w-wind, geopotential, vorticity, PV and DSI inside the determined frontal zones
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#'
#' #front identification using the thermic front parameter (example without front statistic)
#' tfp_fronts=frontid(data$temp,lev_p=data$lev,lat=data$lat,fronts_only=TRUE)
#'
#' #front identification using F diagnostic (example with front statistic)
#' f_fronts=frontid(data$temp,data$u,data$v,data$w,data$z,lev_p=data$lev,lat=data$lat,
#'	method='f',threshold=2,fronts_only=FALSE)
#'
#' #front identification using the dynamic state index (example with statistic)
#' dsi_fronts=frontid(data$temp,data$u,data$v,data$w,data$z,lev_p=data$lev,lat=data$lat,
#'	method='dsi',threshold=4*10^-16,fronts_only=FALSE)
frontid=function(t_fld,u_fld=NULL,v_fld=NULL,w_fld=NULL,phi_fld=NULL,lev_p,lat=NULL,method='tfp',threshold=2*10^-10,dx=0.25,dy=0.25,fronts_only=FALSE,mode='lonlat') {
	# identify frontal zones based on threshold exceedance
	if (method=='tfp') {
		tfp_fld=calc_tfp(t_fld,lev_p,lat,dx,dy,mode)
		fzone=(tfp_fld>threshold)
	} else if (method=='f') {
		fdiag_fld=calc_fdiag(t_fld,u_fld,v_fld,w_fld,lev_p,lat,dx,dy,mode)
		fzone=(fdiag_fld>threshold)
	} else if (method=='dsi') {
		dsi_fld=calc_dsi(t_fld,u_fld,v_fld,w_fld,phi_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=FALSE,pv_fld=NULL)
		fzone=(abs(dsi_fld)>threshold)
	} else {
		message('An unknown method was used. Use either tfp, f, or dsi.')
	}
	# determine meteorological properties of the above detected frontal zones
	if (fronts_only==FALSE) {	
	zeta=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=FALSE,zvort_fld=NULL,mode)
	pv_fld=calc_pv(t_fld,u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=FALSE,zvort_fld=NULL,mode)
	dsi_fld=calc_dsi(t_fld,u_fld,v_fld,w_fld,phi_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative=FALSE,pv_fld=NULL,mode)
	ft = t_fld[fzone]
	fu = u_fld[fzone]
	fv = v_fld[fzone]
	fw = w_fld[fzone]
	fz = phi_fld[fzone]
	fzeta = zeta[,,,3][fzone]
	fpv = pv_fld[fzone]
	fdsi = dsi_fld[fzone]
	return(list(method=method,threshold=threshold,fronts=fzone,temperature=ft,uwind=fu,vwind=fv,wwind=fw,geopotential=fz,vorticity=fzeta,pv=fpv,dsi=fdsi))
	} else {
		return(list(method=method,threshold=threshold,fronts=fzone))
	}
}


#' Petterssen Frontogenesis Function
#'
#' @description Calculates the Petterssen frontogenesis function based on the potential temperature
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @param lat only for lonlat mode: vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#'
#' @return Petterssen Frontogenesis Function
#' @export
calc_frontogenesis=function(t_fld,u_fld,v_fld,w_fld,lev_p,mode='lonlat',lat=NULL,dx=0.25,dy=0.25) {
	th_fld=calc_theta(t_fld,lev_p)
	thgrad=grad(th_fld,d=3,mode='p',rho=NULL,dx=dx,dy=dy,plev=lev_p) #3d gradient
	thgrad_abs=sqrt(thgrad[,,,1]^2+thgrad[,,,2]^2+thgrad[,,,3]^2) #absolute value horizontal gradient
	dth_dx=df_dx(th_fld,lat,dx)
	dth_dy=df_dy(th_fld,dy)
	dth_dp=df_dp(th_fld,lev_p)
	du_dx=df_dx(u_fld,lat,dx)
	du_dy=df_dy(u_fld,dy)
	du_dp=df_dp(u_fld,lev_p)
	dv_dx=df_dx(v_fld,lat,dx)
	dv_dy=df_dy(v_fld,dy)
	dv_dp=df_dp(v_fld,lev_p)
	dw_dx=df_dx(w_fld,lat,dx)
	dw_dy=df_dy(w_fld,dy)
	dw_dp=df_dp(w_fld,lev_p)
	fronto_fld=-(dth_dx*(du_dx*dth_dx+dv_dx*dth_dy+dw_dx*dth_dp)-dth_dy*(du_dy*dth_dx+dv_dy*dth_dy+dw_dy*dth_dp)-dth_dp*(du_dp*dth_dx+dv_dp*dth_dy+dw_dp*dth_dp))/thgrad_abs
	return(fronto_fld)
}

