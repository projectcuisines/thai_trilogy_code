# -*- coding: utf-8 -*-
"""Utilities for the ROCKE3D output."""
import dask.array as da

import xarray as xr

from grid import reverse_along_dim, roll_da_to_pm180
from model_um import calc_um_rel


__all__ = ("adjust_rocke3d_grid", "calc_rocke3d_rei", "calc_rocke3d_rel")

calc_rocke3d_rel = calc_um_rel


def adjust_rocke3d_grid(darr, lon_name="lon", lat_name="lat"):
    """
    Adjust the grid of a ROCKE3D data array.
    Reverse the latitude dimension and shift the substellar coordinate
    from -180 degrees to 0 degree in longitude.
    """
    out = darr
    if lat_name in out.dims:
        out = reverse_along_dim(out, lat_name)
    if lon_name in out.dims:
        # Shift data along the longitude to center the substellar at (0,0)
        out = roll_da_to_pm180(
            out.assign_coords(**{lon_name: out[lon_name] + 180}), lon_name=lon_name
        )
    return out


def calc_rocke3d_rei(rocke3d_ds):
 
    """     Aggregate parametrization based on effective dimension.                                                                                                      
    In the initial form, the same approach is used for stratiform                                                                                                
    and convective cloud.                                                                                                                                        
    
    The fit provided here is based on Stephan Havemann's fit of                                                                                                  
    Dge with temperature, consistent with David Mitchell's treatment                                                                                             
    of the variation of the size distribution with temperature. The                                                                                              
    parametrization of the optical properties is based on De                                                                                                     
    (=(3/2)volume/projected area), whereas Stephan's fit gives Dge                                                                                               
    (=(2*SQRT(3)/3)*volume/projected area), which explains the                                                                                                   
    conversion factor. The fit to Dge is in two sections, because                                                                                                
    Mitchell's relationship predicts a cusp at 216.208 K. Limits                                                                                                 
    of 8 and 124 microns are imposed on Dge: these are based on this                                                                                             
    relationship and should be reviewed if it is changed. Note also                                                                                              
    that the relationship given here is for polycrystals only. 
      
    Parameters
    ----------
    rocke3d_ds: ROCKE-3D data set
    
    These are the parameters used in the temperature dependent 
    parameterizations for ice cloud particle sizes below.
    Parameters for the aggregate parametrization
    a0_agg_cold = 7.5094588E-04,
    b0_agg_cold = 5.0830326E-07,
    a0_agg_warm = 1.3505403E-04,
    b0_agg_warm = 2.6517429E-05,
    t_switch    = 216.208,
    t0_agg      = 279.5,
    s0_agg      = 0.05, 
    
    Returns
    -------
    rei: xarray.DataArray
        Ice effective radius [um], dimension (time,lev,lat,lon).
    """
    
    a0_agg_cold = 7.5094588E-04
    b0_agg_cold = 5.0830326E-07
    a0_agg_warm = 1.3505403E-04
    b0_agg_warm = 2.6517429E-05
    t_switch    = 216.208
    t0_agg      = 279.5
    s0_agg      = 0.05
    
    # Air temperature in ROCKE-3D
    air_temp = rocke3d_ds.t 

    # Calculate the R_eff
    rei = xr.where(air_temp < t_switch,
    
    a0_agg_cold*da.exp(s0_agg *(air_temp - t0_agg)) + b0_agg_cold,
     
    a0_agg_warm*da.exp(s0_agg *(air_temp - t0_agg)) + b0_agg_warm)
    
    # Limit of the parameterization
    rei = (3/2)*(3/(2 * da.sqrt(3))) * xr.ufuncs.minimum(1.24E-04, xr.ufuncs.maximum(8.0E-06,rei))

    rei.rename("ice_cloud_condensate_effective_radius")
    rei.attrs = {"long_name": "ice_cloud_condensate_effective_radius", "units": "micron"}
    return rei
    
def calc_rocke3d_rel(rocke3d_ds, rho = 1.225):   
	"""     Calculation of the effective radius of liquid water
	following: 
	r_eff = ( 3 * rho_air * qcl / ( 4 * pi * rho_water * 0.8 * 1e8 ))**(1./3.)
    corresponding to Eq 14 
    in https://journals.ametsoc.org/jas/article/51/13/1823/23387/The-Measurement-and-Parameterization-of-Effective
    Parameters
    ----------
    liq: xarray.DataArray
    	mass mixing ratio of liquid [kg/kg].
 	rho: xarray.DataArray,
        Array of air density [kg m-3].
    
     Returns
    -------
    rei: xarray.DataArray
        Ice effective radius [um], dimension (time,lev,lat,lon).
    """  
	liq = rocke3d_ds.qcl
	rel = (3 * rho * liq / (4*np.pi*1000*0.8*1e8))**(1./3.)
	
	rel.rename("liquid_cloud_condensate_effective_radius")
	rel.attrs = {"long_name": "liquid_cloud_condensate_effective_radius", "units": "micron"}
	return rel
