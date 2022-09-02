### meteoEVT-package ###
#'
#' Introduction
#' @name meteoEVT-package
#' @docType package
#' @description Energy-Vorticity theory (EVT) is the fundamental theory to describe processes in the atmosphere by combining conserved quantities from hydrodynamics and thermodynamics. The package 'meteoEVT' provides functions to calculate many energetic and vortical quantities, like potential vorticity, Bernoulli function and dynamic state index (DSI) (Weber and Nevir, 2008), for given gridded data, like ERA5 reanalyses. These quantities can be studied directly or can be used for many applications in meteorology, e.g., the objective identification of atmospheric fronts. For this purpose, separate function are provided that allow the detection of fronts based on the thermic front parameter (Hewson, 1998), the F diagnostic (Parfitt et al., 2017) and the DSI (Mack et al., 2022).
#' @details Phenomenons in the Earth's atmosphere, like tropical hurricanes or extratropical cyclones, can adequately be chararcterized by a combination of energetic and vortical quantities. These quantities can also be used for a consistent theoretical description of these phenomenons. This package provides functions to calculate Bernoulli function, vorticity, enstrophy, helicity, Lamb vector and potential vorticity based on given gridded data sets. Addiotionally, by using energy-vortex theory an adiabatic, stationary and invisicid basic state of the Earth's atmosphere can be derived, which is itself a solution of the primitive equations. The derivation from this basic state is given by the dynamic state index (DSI), which can be used for the study of, e.g., cyclones and fronts. Recently, the DSI was used to identify atmospheric fronts objectively from reanalysis data and thereby provides an alternative way for front detection. For this purpose, this package provides funtions to calculate the DSI and use it to identify atmospheric fronts. This method can be compared with state-of-the-art front identification methods based on the thermic front parameter or the F diagnostic.
#' @references
#' * Weber, T. and Névir, P. (2008). Storm tracks and cyclone development using the theoretical concept of the Dynamic State Index (DSI). Tellus A, 60(1):1–10, doi:10.1111/j.1600-0870.2007.00272.x.
#' * Parfitt, R., Czaja, A., and Seo, H. (2017). A simple diagnostic for the detection of atmospheric fronts. Geophys. Res. Lett., 44:4351–4358, doi:10.1002/2017GL073662.
#' * Hewson, T. D. (1998). Objective fronts. Meteorol. Appl., 5:37–65, doi:10.1017/S1350482798000553.
#' * Mack, L., Rudolph, A. and Névir, P. (2022). Identifying atmospheric fronts based on diabatic processes using the dynamic state index (DSI), arXiv:2208.11438.
#' @md
NULL

#' Density
#'
#' @description Calculates the density of an ideal fluid
#' @param t_fld temperature field [K]
#' @param lev_p vector containing pressure levels [Pa]
#' @return density [kg/m^3]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' density=calc_density(data$temp,data$lev)
calc_density=function(t_fld,lev_p) {
	if (dim(t_fld)[3]!=length(lev_p)) {
		message('The third dimension of t_fld and the length of lev_p are not equal.')
	}
	rho_fld=array(NA,dim=dim(t_fld))
	for (i in 1:length(lev_p)) {
		rho_fld[,,i]=lev_p[i]/(Rd()*t_fld[,,i])
	}
	return(rho_fld)
}


#' Potential temperature
#'
#' @description Calculates the potential temperature
#' @param t_fld temperature field [K]
#' @param lev_p vector containing pressure levels [Pa]
#' @return density [kg/m^3]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' theta=calc_theta(data$temp,data$lev)
calc_theta=function(t_fld,lev_p) {
	if (dim(t_fld)[3]!=length(lev_p)) {
		message('The third dimension of t_fld and the length of lev_p are not equal.')
	}
	th_fld=array(NA,dim=dim(t_fld))
	p0=100000 #reference pressure [Pa]
	for (i in 1:length(lev_p)) {
		th_fld[,,i]=t_fld[,,i]*(p0/lev_p[i])^(Rd()/cp())
	}
	return(th_fld)
}


#' Bernoulli function
#'
#' @description Calculates the Bernoulli function, i.e. total energy density, as sum of potential, kinetic and thermal energy density
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param phi_fld geopotential height [gpm]
#' @return Bernoulli function field [m^2/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' bernoulli=calc_bernoulli(data$temp,data$u,data$v,data$w,data$z)
calc_bernoulli=function(t_fld,u_fld,v_fld,w_fld,phi_fld) {
	#if (dim(t_fld)!=dim(u_fld) | dim(t_fld)!=dim(v_fld) | dim(t_fld)!=dim(w_fld) | dim(t_fld)!=dim(phi_fld)) {
	#	message('All input fields must have the same dimensions.')
	#}
	return(cp()*t_fld+(u_fld^2+v_fld^2+w_fld^2)/2+phi_fld)
}


#' Vorticity
#'
#' @description Calculates the vorticity
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param zvort_only logical, TRUE: if only the vertical vorticity (zvort) should be calculated, FALSE: for the whole vorticity vector, default: FALSE
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity, default: FALSE
#' @param zvort_fld optional zvort field (if e.g., zvort is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return vorticity field [1/s]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' #3d vorticity
#' xi=calc_vorticity(data$u,data$v,data$w,data$lev,lat=data$lat)
#' #z-vorticity as scalar
#' zeta=calc_vorticity(data$u,data$v,data$w,data$lev,lat=data$lat,zvort_only=TRUE)
calc_vorticity=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,zvort_only=FALSE,relative=FALSE,zvort_fld=NULL,mode='lonlat') {
	dims=dim(u_fld)
	if (zvort_only == TRUE) {
		if (is.null(zvort_fld)) {
			zvort_fld=df_dx(v_fld,lat,dx)-df_dy(u_fld,dy)
			return(zvort_fld)
		} else {
			return(zvort_fld)
		}
	} else {
		du_dx=df_dx(u_fld,lat,dx,mode)
		du_dy=df_dy(u_fld,dy,mode)
		du_dp=df_dp(u_fld,lev_p)
		dv_dx=df_dx(v_fld,lat,dx,mode)
		dv_dy=df_dy(v_fld,dy,mode)
		dv_dp=df_dp(v_fld,lev_p)
		dw_dx=df_dx(w_fld,lat,dx,mode)
		dw_dy=df_dy(w_fld,dy,mode)
		dw_dp=df_dp(w_fld,lev_p)
		vort_fld=array(NA,dim=c(dim(u_fld),3))
		vort_fld[,,,1]=du_dy-dv_dp
		vort_fld[,,,2]=du_dp-dw_dx
		if (is.null(zvort_fld)){
			zvort_fld=df_dx(v_fld,lat,dx,mode)-df_dy(u_fld,dy,mode)
		}
		if (relative==FALSE) { #add planetary vorticity
			for (j in 1:dims[2]) {
				vort_fld[,j,,2]=vort_fld[,j,,2]+2*omega()*cos(pi/180*lat[j])*matrix(1,nrow=dims[1],ncol=dims[3])
				zvort_fld[,j,]=zvort_fld[,j,]+2*omega()*sin(pi/180*lat[j])*matrix(1,nrow=dims[1],ncol=dims[3])		
			}
		}
		vort_fld[,,,3]=zvort_fld
		return(vort_fld)
	}
}


#' Enstrophy density
#'
#' @description Calculates the enstrophy density (vorticity squared) either in 2d or 3d
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param zvort_only logical, TRUE: if only 2d enstrophy (based on z-vorticity) should be calculated, FALSE: for 3d enstrophy (based on 3d vorticity), default: TRUE
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity should be used for calculation of enstrophy, default: TRUE
#' @param zvort_fld optional zvort field (if e.g., zvort is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return enstrophy density field [1/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' #3d enstropy
#' ens3d=calc_enstrophy(data$u,data$v,data$w,data$lev,lat=data$lat)
#' #2d enstropy as scalar
#' ens2d=calc_enstrophy(data$u,data$v,lev_p=data$lev,lat=data$lat,zvort_only=TRUE)
calc_enstrophy=function(u_fld,v_fld,w_fld=NULL,lev_p,lat=NULL,dx=0.25,dy=0.25,zvort_only=TRUE,relative=TRUE,zvort_fld=NULL,mode='lonlat') {
	vort_fld=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only,relative,zvort_fld,mode)
	if (zvort_only==TRUE) { #only 2d enstrophy
		return(vort_fld^2)
	} else { #3d enstrophy
		return(vort_fld[,,,1]^2+vort_fld[,,,2]^2+vort_fld[,,,3]^2)
	}
}

#' Helicity density
#'
#' @description Calculates the helicity density (scalar product of wind vector and vorticity vector) either for the whole vector (3d) or only for the vertical component (updraft helicity)
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param vert_only logical, TRUE: if only the updraft helicity w*zeta (based on z-vorticity) should be calculated, FALSE: for 3d helicity (based on 3d vorticity), default: FALSE
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity should be used for calculation of enstrophy, default: TRUE
#' @param zvort_fld optional zvort field (if e.g., zvort is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return helicity density field [m/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' #3d helicity
#' hel=calc_helicity(data$u,data$v,data$w,data$lev,lat=data$lat)
#' #updraft helicity
#' up_hel=calc_helicity(data$u,data$v,data$w,data$lev,lat=data$lat,vert_only=TRUE)
calc_helicity=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,vert_only=FALSE,relative=TRUE,zvort_fld=NULL,mode='lonlat') {
	vort_fld=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=vert_only,relative,zvort_fld,mode)
	if (vert_only==TRUE) { #only 2d enstrophy
		return(vort_fld*w_fld)
	} else {
		return(vort_fld[,,,1]*u_fld+vort_fld[,,,2]*v_fld+vort_fld[,,,3]*w_fld)
	}
}

#' Lamb vector (sometimes called vortex energy)
#'
#' @description Calculates the Lamb vector (cross product of wind vector and vorticity vector)
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity should be used for calculation of enstrophy, default: TRUE
#' @param zvort_fld optional zvort field (if e.g., zvort is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return lamb vector [m/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' lamb=calc_lamb(data$u,data$v,data$w,data$lev,lat=data$lat)
calc_lamb=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,relative=TRUE,zvort_fld=NULL,mode='lonlat') {
	vort_fld=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only=FALSE,relative,zvort_fld,mode)
	wind_fld=array(NA,dim=c(dim(u_fld),3))
	wind_fld[,,,1]=u_fld
	wind_fld[,,,2]=v_fld
	wind_fld[,,,3]=w_fld
	return(crossprod(wind_fld,vort_fld))
}


#' Potential Vorticity (PV)
#'
#' @description Calculates the potential vorticity
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param zvort_only logical, TRUE: if only the vertical vorticity (zvort) should be calculated, FALSE: for the whole vorticity vector, default: FALSE
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity, default: FALSE
#' @param zvort_fld optional zvort field (if e.g., zvort is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return potential vorticity field [K*m^2/(kg*s)]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' #PV based on all three components
#' pv=calc_pv(data$temp,data$u,data$v,data$w,data$lev,lat=data$lat)
#' #PV only based on vertical component
#' pv_vert=calc_pv(data$temp,data$u,data$v,data$w,lev_p=data$lev,lat=data$lat,zvort_only=TRUE)
calc_pv=function(t_fld,u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,zvort_only=FALSE,relative=FALSE,zvort_fld=NULL,mode='lonlat') {
	th_fld=calc_theta(t_fld,lev_p)
	dth_dx=df_dx(th_fld,lat,dx)
	dth_dy=df_dy(th_fld,dy)
	dth_dp=df_dp(th_fld,lev_p)
	vort=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only,relative,zvort_fld,mode)
	rho=calc_density(t_fld,lev_p)
	diab=-rho*g() #diabatic coefficient, for converting p- to z-system
	if (zvort_only==TRUE) {
		return(dth_dp*vort/rho*diab)
	}else{
		return(1/rho*(dth_dx*vort[,,,1]+dth_dy*vort[,,,2]+dth_dp*vort[,,,3]*diab))
	}
		
}

#' Dynamic State Index (DSI)
#'
#' @description Calculates the dynamic state index DSI
#' @param t_fld temperature field [K]
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param phi_fld geopotential height [gpm]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param zvort_only logical, TRUE: if only the vertical vorticity (zvort) should be calculated, FALSE: for the whole vorticity vector, default: FALSE
#' @param relative logical, TRUE: only relative vorticity, FALSE: whole (absolute) vorticity, default: FALSE
#' @param pv_fld optional pv field (if e.g., PV is directly taken from ERA5 and not calculated separately)
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return dynamic state index [K^2*m^4/(kg^2*s^3)]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' dsi=calc_dsi(data$temp,data$u,data$v,data$w,data$z,lev_p=data$lev,lat=data$lat)
calc_dsi=function(t_fld,u_fld,v_fld,w_fld,phi_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,zvort_only=FALSE,relative=FALSE,pv_fld=NULL,mode='lonlat') {
	th_fld=calc_theta(t_fld,lev_p)
	b_fld=calc_bernoulli(t_fld,u_fld,v_fld,w_fld,phi_fld)
	pv_fld=calc_pv(t_fld,u_fld,v_fld,w_fld,lev_p,lat,dx,dy,zvort_only,relative,zvort_fld=NULL,mode=mode)
	dth_dx=df_dx(th_fld,lat,dx)
	dth_dy=df_dy(th_fld,dy)
	dth_dp=df_dp(th_fld,lev_p)
	db_dx=df_dx(b_fld,lat,dx)
	db_dy=df_dy(b_fld,dy)
	db_dp=df_dp(b_fld,lev_p)
	dpv_dx=df_dx(pv_fld,lat,dx)
	dpv_dy=df_dy(pv_fld,dy)
	dpv_dp=df_dp(pv_fld,lev_p)
	rho=calc_density(t_fld,lev_p)
	diab=-rho*g() #diabatic coefficient, for converting p- to z-system
	#dsi_fld=1/rho*diab*jacobian(th_fld,b_fld,pv_fld,lat,d=3,mode='p',rho=NULL,dy=,dx,plev=lev_p)$jacobian
	dsi_fld=1/rho*(dth_dx*db_dy*dpv_dp+dth_dy*db_dp*dpv_dx+dth_dp*db_dx*dpv_dy-dth_dp*db_dy*dpv_dx-dth_dy*db_dx*dpv_dp-dth_dx*db_dp*dpv_dy)*diab
	return(dsi_fld)
}
