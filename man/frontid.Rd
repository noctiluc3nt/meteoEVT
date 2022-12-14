% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_fronts.R
\name{frontid}
\alias{frontid}
\title{Front Identification und Statistics}
\usage{
frontid(
  t_fld,
  u_fld = NULL,
  v_fld = NULL,
  w_fld = NULL,
  phi_fld = NULL,
  lev_p,
  lat = NULL,
  method = "tfp",
  threshold = 2 * 10^-10,
  dx = 0.25,
  dy = 0.25,
  fronts_only = FALSE,
  mode = "lonlat"
)
}
\arguments{
\item{t_fld}{temperature field [K]}

\item{u_fld}{zonal velocity field [m/s]}

\item{v_fld}{meridional velocity field [m/s]}

\item{w_fld}{vertical velocity field [m/s]}

\item{phi_fld}{geopotential height [gpm]}

\item{lev_p}{vector containing pressure levels [Pa]}

\item{lat}{only for lonlat mode: vector containing latitude}

\item{method}{character containing the method, use \code{'tfp'} for TFP method, \code{'f'} for F diagnostic and \code{'dsi'} for DSI method}

\item{threshold}{scalar containing a suitable threshold (e.g., 2*10^-10 for TFP method, 1 or for F diagnostic, 10^-16 for DSI method)}

\item{dx}{x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})}

\item{dy}{y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})}

\item{fronts_only}{if you only want to calculate the frontal regions and not their properties (default FALSE)}

\item{mode}{the horizontal coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid}
}
\value{
list containing the used method and used threshold, field with logicals containing the detected frontal zones and numerics of temperature, u-wind, v-wind, w-wind, geopotential, vorticity, PV and DSI inside the determined frontal zones
}
\description{
Calculates frontal zones based on a chosen method (TFP, F diagnostic, DSI) and provides statistics of the distribution of meteorological quantities inside the determined frontak zones.
}
\examples{
myfile=system.file("extdata", "era5_storm-zeynep.nc", package = "meteoEVT")
data = readin_era5(myfile)

#front identification using the thermic front parameter (example without front statistic)
tfp_fronts=frontid(data$temp,lev_p=data$lev,lat=data$lat,fronts_only=TRUE)

#front identification using F diagnostic (example with front statistic)
f_fronts=frontid(data$temp,data$u,data$v,data$w,data$z,lev_p=data$lev,lat=data$lat,
method='f',threshold=2,fronts_only=FALSE)

#front identification using the dynamic state index (example with statistic)
dsi_fronts=frontid(data$temp,data$u,data$v,data$w,data$z,lev_p=data$lev,lat=data$lat,
method='dsi',threshold=4*10^-16,fronts_only=FALSE)
}
