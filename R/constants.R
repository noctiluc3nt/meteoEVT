#### Constants #####
#' R
#'
#' ideal gas constant [J/(mol*K)]
#' @keywords internal
R=function(){
	return(8.314462)
}

#' Rd
#'
#' specific gas constant for dry air [J/(kg*K)]
#' @keywords internal
Rd=function(){
	return(287.058)
}

#' Rv
#'
#' specific gas constant for water vapor [J/(kg*K)]
#' @keywords internal
Rv=function(){
	return(461.4)
}

#' re
#'
#' Earth's radius [m]
#' @keywords internal
re=function(){
	return(6371008.767)
}

#' cp
#'
#' heat capacity for constant pressure [J/(kg*K)]
#' @keywords internal
cp=function(){
	return(1005.7)
}

#' g
#'
#' gravitional accelaration [m/s^2]
#' @keywords internal
g=function(){
	return(9.8062)
}

#' omega
#'
#' Earth's rotation frequency [1/s]
#' @keywords internal
omega=function(){
	return(2*pi/(24*3600))
}



