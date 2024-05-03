#' Vorticity tensor
#'
#' @description Calculates the vorticity tensor (Omega matrix)
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return vorticity tensor (3x3 tensor) [1/s]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' Omega=calc_vorticity_tensor(data$u,data$v,data$w,data$lev,lat=data$lat)
#'
calc_vorticity_tensor=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
	du_dx=df_dx(u_fld,lat,dx,mode)
	du_dy=df_dy(u_fld,dy,mode)
	du_dp=df_dp(u_fld,lev_p)
	dv_dx=df_dx(v_fld,lat,dx,mode)
	dv_dy=df_dy(v_fld,dy,mode)
	dv_dp=df_dp(v_fld,lev_p)
	dw_dx=df_dx(w_fld,lat,dx,mode)
	dw_dy=df_dy(w_fld,dy,mode)
	dw_dp=df_dp(w_fld,lev_p)
    a21=dv_dx-du_dy
    a31=dw_dx-du_dp
    a32=dw_dy-dv_dp
    mat=matrix(0,ncol=3,nrow=3)
    mat[2,1]=a21
    mat[1,2]=-a21
    mat[3,1]=a31
    mat[1,3]=-a31
    mat[3,2]=a32
    mat[2,3]=-a32
	return(mat)
}

#' Strain-rate tensor
#'
#' @description Calculates the strain-rate tensor (S matrix)
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return Strain-rate tensor (3x3 tensor) [1/s]
#' @export
#'
#' @examples
#' myfile = system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' S = calc_strain_rate_tensor(data$u,data$v,data$w,data$lev,lat=data$lat)
#'
calc_strain_rate_tensor=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
	du_dx=df_dx(u_fld,lat,dx,mode)
	du_dy=df_dy(u_fld,dy,mode)
	du_dp=df_dp(u_fld,lev_p)
	dv_dx=df_dx(v_fld,lat,dx,mode)
	dv_dy=df_dy(v_fld,dy,mode)
	dv_dp=df_dp(v_fld,lev_p)
	dw_dx=df_dx(w_fld,lat,dx,mode)
	dw_dy=df_dy(w_fld,dy,mode)
	dw_dp=df_dp(w_fld,lev_p)
    a21=du_dy+dv_dx
    a31=dw_dx+du_dp
    a32=dw_dy+dv_dp
    mat=matrix(0,ncol=3,nrow=3)
    mat[2,1]=a21
    mat[1,2]=a21
    mat[3,1]=a31
    mat[1,3]=a31
    mat[3,2]=a32
    mat[2,3]=a32
    mat[1,1]=2*du_dx
    mat[2,2]=2*dv_dy
    mat[3,3]=2*dw_dp
	return(mat)
}

#' Q-invariant
#'
#' @description Calculates the Q-invariant := 1/2 * ( frobenius_norm(Omega)^2 - frobenius_norm(S)^2 )
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return Q-invariant [1/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' Q=calc_Qinvariant(data$u,data$v,data$w,data$lev,lat=data$lat)
#'
calc_Qinvariant = function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
    Omega=calc_vorticity_tensor(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=dx,dy=dy,mode=mode)
    S=calc_strain_rate_tensor(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=dx,dy=dy,mode=mode)
    return(1/2*(frobenius_norm(Omega)^2-frobenius_norm(S)^2))
}

#' Kinematic vorticity number
#'
#' @description Calculates kinematic vorticity number = frobenius_norm(Omega)/frobenius_norm(S)
#' @param u_fld zonal velocity field [m/s]
#' @param v_fld meridional velocity field [m/s]
#' @param w_fld vertical velocity field [m/s]
#' @param lev_p vector containing pressure levels [Pa]
#' @param lat vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode use 'lonlat' if the data is given on a lon-lat-grid or 'cartesian' if the data is given on an equidistant cartesian grid
#' @return kinematic vorticity number [-]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' eta=calc_kinematic_vorticity_number(data$u,data$v,data$w,data$lev,lat=data$lat)
#'
calc_kinematic_vorticity_number=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat') {
    Omega=calc_vorticity_tensor(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat')
    S=calc_strain_rate_tensor(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,mode='lonlat')
    return(frobenius_norm(Omega)/frobenius_norm(S))
}

#' Hydrodynamic charge
#'
#' @description Calculates the hydrodynamic charge (density), i.e. the divergence of the Lamb vector
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
#' @return hydrodynamic charge [1/s^2]
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' qh=calc_hydrodynamic_charge(data$u,data$v,data$w,data$lev,lat=data$lat)
calc_hydrodynamic_charge=function(u_fld,v_fld,w_fld,lev_p,lat=NULL,dx=0.25,dy=0.25,relative=TRUE,zvort_fld=NULL,mode='lonlat') {
	vort_fld=calc_vorticity(u_fld,v_fld,w_fld,lev_p,lat,dx,dy,FALSE,relative,zvort_fld,mode)
	wind_fld=array(NA,dim=c(dim(u_fld),3))
	wind_fld[,,,1]=u_fld
	wind_fld[,,,2]=v_fld
	wind_fld[,,,3]=w_fld
	return(div(crossprod(wind_fld,vort_fld)))
}