#### partial derivatives ####
#' df_dx
#'
#' @description Calculates the x derivative using central differences (for lonlat-grid or cartesian grid)
#' @param fld field with dimensions (lon,lat,p)
#' @param lat only for lonlat mode: vector containing latitude
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return field containing the partial derivative w.r.t. x
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' theta=calc_theta(data$temp,data$lev)
#' dtheta_dx=df_dx(theta,data$lat)
df_dx=function(fld,lat=NULL,dx=0.25,mode='lonlat') {
	dims=dim(fld)
	if (mode=='lonlat') {
		if (dims[2]!=length(lat)) {
			message('The latitude dimensions do not match.')
		}
		fld_dx=array(NA,dim=dims)
		#boundaries
		fld_dx[1,,]=2*(fld[2,,]-fld[1,,])
		fld_dx[dims[1],,]=2*(fld[dims[1],,]-fld[dims[1]-1,,])
		fld_dx[2:(dims[1]-1),,]=fld[3:dims[1],,]-fld[1:(dims[1]-2),,]
		for (i in 1:dims[2]) {
			fld_dx[,i,]=fld_dx[,i,]/(2*pi/180*re()*cos(lat[i]*pi/180)*dx)
		}
		return(fld_dx)
	} else if (mode=='cartesian') {
		fld_dx=array(NA,dim=dims)
		#boundaries
		fld_dx[1,,]=(fld[2,,]-fld[1,,])/dx
		fld_dx[dims[1],,]=(fld[dims[1],,]-fld[dims[1]-1,,])/dx
		fld_dx[2:(dims[1]-1),,]=(fld[3:dims[1],,]-fld[1:(dims[1]-2),,])/(2*dx)
		return(fld_dx)
	} else {
		message('An unkown mode is used. You can use lonlat or cartesian.')
	}
}


#' df_dy
#'
#' @description Calculates the y derivative using central differences
#' @param fld with dimensions (lon,lat,p)
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return field containing the partial derivative w.r.t. y
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' theta=calc_theta(data$temp,data$lev)
#' dtheta_dy=df_dy(theta,dy=0.25)
df_dy=function(fld,dy=0.25,mode='lonlat') {
	dims=dim(fld)
	fld_dy=array(NA,dim=dims)
	dy=dy/360*2*pi*re()
	#boundaries
	fld_dy[,1,]=(fld[,2,]-fld[,1,])/dy
	fld_dy[,dims[2],]=(fld[,dims[2],]-fld[,dims[2]-1,])/dy
	fld_dy[,2:(dims[2]-1),]=(fld[,3:dims[2],]-fld[,1:(dims[2]-2),])/(2*dy)
	return(-fld_dy)
}


#' df_dp
#'
#' @description Calculates the p derivative (pressure system) using central differences
#' @param fld field with dimensions (lon,lat,p)
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @return field containing the partial derivative w.r.t. p
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' theta=calc_theta(data$temp,data$lev)
#' dtheta_dp=df_dp(theta)
df_dp=function(fld,plev=5000) {
	dims=dim(fld)
	fld_dp=array(NA,dim=dims)
	if (length(plev)==1) {
		dp=plev
		#boundaries
		fld_dp[,,1]=(fld[,,2]-fld[,,1])/dp
		fld_dp[,,dims[3]]=(fld[,,dims[3]]-fld[,,dims[3]-1])/dp
		fld_dp[,,2:(dims[3]-1)]=(fld[,,3:dims[3]]-fld[,,1:(dims[3]-2)])/(2*dp)
	} else {
		if (dims[3]!=length(plev)) {
			message('The pressure dimensions do not match.')
		}
		#boundaries
		fld_dp[,,1]=(fld[,,2]-fld[,,1])/(plev[2]-plev[1])
		fld_dp[,,dims[3]]=(fld[,,dims[3]]-fld[,,dims[3]-1])/(plev[dims[3]]-plev[dims[3]-1])
		for (i in 2:(dims[3]-1)) {
			fld_dp[,,i]=(fld[,,i+1]-fld[,,i-1])/(plev[i+1]-plev[i-1])
		}
	}
	return(fld_dp)
}	


#' df_dz
#'
#' @description Calculates the z derivative 
#' @param fld field with dimensions (lon,lat,p)
#' @param rho field with dimensions (lon,lat,p) for density or a scalar rho (for constant density)
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @return field containing the partial derivative w.r.t. z
#' @export
df_dz=function(fld,rho,plev=5000) {
	if (dim(fld)!=dim(rho) & dim(rho)!=1) {
		message('The dimension of the input fields do not match.')
	}
	return(-g()/rho*df_dp(fld,plev))
}


#' gradient of a scalar field
#'
#' @description Calculates the gradient
#' @param fld field with dimensions (lon,lat,p)
#' @param lat vector containing latitude
#' @param d scalar for dimension (use d=2 for horizontal gradient and d=3 for 3d-gradient)
#' @param system for type of coordinate system (use 'p' for pressure system and 'z' for height system)
#' @param rho field with dimensions (lon,lat,p) for density or a scalar rho (for constant density)
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return field containing the gradient with dimension (lon,lat,p,d)
#' @export
#'
#' @examples
#' myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
#' data = readin_era5(myfile)
#' theta=calc_theta(data$temp,data$lev)
#' theta_grad=grad(theta,data$lat)
grad=function(fld,lat=NULL,d=3,system='p',rho=NULL,dx=0.25,dy=0.25,plev=5000,mode='lonlat') {
	fld_grad=array(NA,dim=c(dim(fld),d))
	fld_grad[,,,1]=df_dx(fld,lat,dx)
	fld_grad[,,,2]=df_dy(fld,dy)
	if (d==2) {
		return(fld_grad) 
	} else if (d==3 & mode=='p') {
		fld_grad[,,,3]=df_dp(fld,plev)
		return(fld_grad)
	} else if (d==3 & mode=='z') {
		if (is.null(rho)) {
			message('The input density is still NULL.')
		}
		fld_grad[,,,3]=df_dz(fld,rho,plev)
		return(fld_grad)
	} else {
		message('Mode and dimension are not consistent.')
	}
}


#' divergence 
#'
#' @description Calculates the divergence of a vector field
#' @param fld field with dimensions (lon,lat,p,d)
#' @param lat vector containing latitude (only for \code{mode='lonlat'})
#' @param d scalar for dimension (use d=2 for horizontal gradient and d=3 for 3d-gradient)
#' @param system for type of coordinate system (use 'p' for pressure system and 'z' for height system)
#' @param rho field with dimensions (lon,lat,p) for density or a scalar rho (for constant density)
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return field containing the divergence of fld
#' @export
div=function(fld,lat=NULL,d=3,system='p',rho=NULL,dx=0.25,dy=0.25,plev=5000,mode='lonlat') {
	fld_div=array(0,dim=dim(fld))
	fld_div=df_dx(fld,lat,dx)+df_dy(fld,dy)
	if (d==2) {
		return(fld_div) 
	} else if (d==3 & mode=='p') {
		fld_div=fld_div+df_dp(fld,plev)
		return(fld_div)
	} else if (d==3 & mode=='z') {
		if (is.null(rho)) {
			message('The input density is still NULL.')
		}
		fld_div=fld_div+df_dz(fld,rho,plev)
		return(fld_div)
	} else {
		message('Mode and dimension are not consistent.')
	}
}


#' rotation 
#'
#' @description Calculates the rotation of a vector field
#' @param fld with dimensions (lon,lat,p,d)
#' @param lat vector containing latitude
#' @param d scalar for dimension (use d=2 for horizontal gradient and d=3 for 3d-gradient)
#' @param system for type of coordinate system (use 'p' for pressure system and 'z' for height system)
#' @param rho field with dimensions (lon,lat,p) for density or a scalar rho (for constant density)
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return field containing the divergence of fld
#' @export
rot=function(fld,lat=NULL,d=3,system='p',rho=NULL,dx=0.25,dy=0.25,plev=5000,mode='lonlat') {
	fld_rot=array(NA,dim=c(dim(fld),d))
	if (d==2) { #invertibility condition
		fld_rot[,,,1]=df_dx(fld[,,,3],lat=lat,dx=dx)
		fld_rot[,,,2]=df_dy(fld[,,,1],dy=dy)
		return(fld_rot)
	} else if (d==3 & mode=='p') {
		fld_rot[,,,1]=df_dy(fld[,,,3],dy=dy)-df_dp(fld[,,,2],plev)
		fld_rot[,,,2]=-(df_dx(fld[,,,3],lat=lat,dx=dx)-df_dp(fld[,,,1],plev))
		fld_rot[,,,3]=df_dx(fld[,,,2],lat,dx=dx)-df_dy(fld[,,,1],dy=dy)
		return(fld_rot)
	} else if (d==3 & mode=='z') {
		if (is.null(rho)) {
			message('The input density is still NULL.')
		}
		fld_rot[,,,1]=df_dy(fld[,,,3],dy)-df_dz(fld[,,,2],rho,plev)
		fld_rot[,,,2]=-(df_dx(fld[,,,3],lat,dx)-df_dz(fld[,,,1],rho,plev))
		fld_rot[,,,3]=df_dx(fld[,,,2],lat,dx)-df_dy(fld[,,,1],dy)
		return(fld_rot)
	} else {
		message('Mode and dimension are not consistent.')
	}
}

#' Jacobian matrix and determinant
#'
#' @description Calculates the Jacobian matrix and Jacobian determinant for 2 or 3 given scalar fields
#' @param fld1 field 1 with dimensions (lon,lat,p)
#' @param fld2 field 2 with dimensions (lon,lat,p)
#' @param fld3 field 3 with dimensions (lon,lat,p)
#' @param lat vector containing latitude
#' @param d scalar for dimension (use d=2 for 2 input fields and d=3 for 3 inpt fields)
#' @param system for type of coordinate system (use 'p' for pressure system and 'z' for height system)
#' @param rho field with dimensions (lon,lat,p) for density or a scalar rho (for constant density)
#' @param dx x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param dy y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})
#' @param plev a scalar containing the p resolution (if equidistant) or a vector containing pressure levels in Pa (for non-equidistant)
#' @param mode the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid
#' @return list containing Jacobian matrix and determinant
#' @export
jacobian=function(fld1,fld2,fld3=NULL,lat=NULL,d=3,system='p',rho=NULL,dx=0.25,dy=0.25,plev=5000,mode='lonlat') {
	jmat=array(NA,dim=c(dim(fld1),d,d))
	if (d==2) { #2d matrix
		jmat[,,,1,1]=df_dx(fld1,lat,dx)
		jmat[,,,2,1]=df_dx(fld2,lat,dx)
		jmat[,,,1,2]=df_dy(fld1,dy)
		jmat[,,,2,2]=df_dy(fld2,dy)
		det=jmat[,,,1,1]*jmat[,,,2,2]-jmat[,,,1,2]*jmat[,,,2,1]
		return(list('jacobian'=jmat,'determinant'=det,'dimensions'=d))
	} else if (d==3 & mode=='p') {
		jmat[,,,1,1]=df_dx(fld1,lat,dx)
		jmat[,,,2,1]=df_dx(fld2,lat,dx)
		jmat[,,,3,1]=df_dx(fld3,lat,dx)
		jmat[,,,1,2]=df_dy(fld1,dy)
		jmat[,,,2,2]=df_dy(fld2,dy)
		jmat[,,,3,2]=df_dy(fld3,dy)
		jmat[,,,1,3]=df_dp(fld1,plev)
		jmat[,,,2,3]=df_dp(fld2,plev)
		jmat[,,,3,3]=df_dp(fld3,plev)
		det=jmat[,,,1,1]*jmat[,,,2,2]*jmat[,,,3,3]+jmat[,,,1,2]*jmat[,,,2,3]*jmat[,,,3,1]+jmat[,,,1,3]*jmat[,,,2,1]*jmat[,,,3,2]-jmat[,,,1,3]*jmat[,,,2,2]*jmat[,,,3,1]-jmat[,,,1,2]*jmat[,,,2,1]*jmat[,,,3,3]-jmat[,,,1,1]*jmat[,,,2,3]*jmat[,,,3,2]
		return(list('jacobian'=jmat,'determinant'=det,'dimensions'=d,'mode'=mode))
	} else if (d==3 & mode=='z') {
		if (is.null(rho)) {
			message('The input density is still NULL.')
		}
		jmat[,,,1,1]=df_dx(fld1,lat,dx)
		jmat[,,,2,1]=df_dx(fld2,lat,dx)
		jmat[,,,3,1]=df_dx(fld3,lat,dx)
		jmat[,,,1,2]=df_dy(fld1,dy)
		jmat[,,,2,2]=df_dy(fld2,dy)
		jmat[,,,3,2]=df_dy(fld3,dy)
		jmat[,,,1,3]=df_dz(fld1,rho,plev)
		jmat[,,,2,3]=df_dz(fld2,rho,plev)
		jmat[,,,3,3]=df_dz(fld3,rho,plev)
		det=jmat[,,,1,1]*jmat[,,,2,2]*jmat[,,,3,3]+jmat[,,,1,2]*jmat[,,,2,3]*jmat[,,,3,1]+jmat[,,,1,3]*jmat[,,,2,1]*jmat[,,,3,2]-jmat[,,,1,3]*jmat[,,,2,2]*jmat[,,,3,1]-jmat[,,,1,2]*jmat[,,,2,1]*jmat[,,,3,3]-jmat[,,,1,1]*jmat[,,,2,3]*jmat[,,,3,2]
		return(list('jacobian'=jmat,'determinant'=det,'dimensions'=d,'mode'=mode))
	} else {
		message('Mode and dimension are not consistent.')
	}
}


#' scalar product
#'
#' @description Calculates the scalar product of two given fields
#' @param fld1 field 1 with dimensions (lon,lat,p,d)
#' @param fld2 field 2 with dimensions (lon,lat,p,d)
#' @return field of the scalar product with dimensions (lon,lat,p)
#' @export
scalarprod=function(fld1,fld2) {
	if (dim(fld1)!=dim(fld2)) {
		message('Input fields have different dimensions.')
	}
	dims=dim(fld1)
	if (length(dims)==3) { #i.e. both are scalar fields
		return(fld1*fld2)
	} else {
		fld_sp=array(0,dim=dims[1:3])
		for (i in 1:dims[4]) {
			fld_sp=fld_sp+fld1[,,,i]*fld2[,,,i]
		}
		return(fld_sp)
	}
}

#' cross product
#'
#' @description Calculates the cross product of two given 3d vector fields
#' @param fld1 field 1 with dimensions (lon,lat,p,3)
#' @param fld2 field 2 with dimensions (lon,lat,p,3)
#' @return field containing the cross product
#' @export
crossprod=function(fld1,fld2) {
	dims=dim(fld1)
	if (dims[4]!=3) {
		message('The input fields must have the dimensions (lon,lat,p,3).')
	}	
	fld_cp=array(0,dim=dims)
	fld_cp[,,,1]=fld1[,,,2]*fld2[,,,3]-fld1[,,,3]*fld2[,,,2]
	fld_cp[,,,2]=-(fld1[,,,1]*fld2[,,,3]-fld1[,,,3]*fld2[,,,1])
	fld_cp[,,,3]=fld1[,,,1]*fld2[,,,2]-fld1[,,,2]*fld2[,,,1]
	return(fld_cp)
}


