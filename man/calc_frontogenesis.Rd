% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_fronts.R
\name{calc_frontogenesis}
\alias{calc_frontogenesis}
\title{Petterssen Frontogenesis Function}
\usage{
calc_frontogenesis(
  t_fld,
  u_fld,
  v_fld,
  w_fld,
  lev_p,
  mode = "lonlat",
  lat = NULL,
  dx = 0.25,
  dy = 0.25
)
}
\arguments{
\item{t_fld}{temperature field [K]}

\item{u_fld}{zonal velocity field [m/s]}

\item{v_fld}{meridional velocity field [m/s]}

\item{w_fld}{vertical velocity field [m/s]}

\item{lev_p}{vector containing pressure levels [Pa]}

\item{mode}{the coordinate system, options are lonlat for a longitude-latitude-grid (default), or cartesian for an equidistant cartesian grid}

\item{lat}{only for lonlat mode: vector containing latitude}

\item{dx}{x resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})}

\item{dy}{y resolution in the corresponding unit (e.g. 0.25 degree for ERA5 with \code{mode='lonlat'} or e.g. 1000 m in cartesian coordinates with \code{mode='cartesian'})}
}
\value{
Petterssen Frontogenesis Function
}
\description{
Calculates the Petterssen frontogenesis function based on the potential temperature
}
