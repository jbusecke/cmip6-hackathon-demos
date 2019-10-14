from xgcm.autogenerate import generate_grid_ds
from xgcm import Grid
import xarray as xr
import pyproj
    
def distance(lon0, lat0, lon1, lat1):
    geod = pyproj.Geod(ellps='WGS84')
    _, _, distance = geod.inv(lon0, lat0, lon1, lat1)
    return distance

# TODO: check if all these hardcoded lines work for each CMIP6 model...
def recreate_full_grid(ds, lon_name="lon", lat_name="lat"):
    ds_full = generate_grid_ds(ds, {"X": "x", "Y": "y"}, position=("center", "right"))
    grid = Grid(ds_full, periodic=['X'])

    # infer dx at eastern bound from tracer points
    lon0, lon1 = grid.axes["X"]._get_neighbor_data_pairs(ds_full.lon.load(), "right")
    lat0 = lat1 = ds_full.lat.load().data
    dx = distance(lon0, lat0, lon1, lat1)
    ds_full.coords["dxe"] = xr.DataArray(
        dx, coords=grid.interp(ds_full.lon, "X").coords
    )

    # infer dy at northern bound from tracer points
    lat0, lat1 = grid.axes["Y"]._get_neighbor_data_pairs(
        ds_full.lat.load(), "right", boundary="extrapolate"
    )
    
    lon0 = lon1 = ds_full.lon.load().data
    dy = distance(lon0, lat0, lon1, lat1)
    ds_full.coords["dyn"] = xr.DataArray(
        dy, coords=grid.interp(ds_full.lat, "Y", boundary="extrapolate").coords
    )

    # now simply interpolate all the other metrics
    ds_full.coords['dxt'] = grid.interp(ds_full.coords['dxe'], 'X')
    ds_full.coords['dxne'] = grid.interp(ds_full.coords['dxe'], 'Y', boundary="extrapolate")
    ds_full.coords['dxn'] = grid.interp(ds_full.coords['dxt'], 'Y', boundary="extrapolate")

    ds_full.coords['dyt'] = grid.interp(ds_full.coords['dyn'], 'Y', boundary="extrapolate")
    ds_full.coords['dyne'] = grid.interp(ds_full.coords['dyn'], 'X')
    ds_full.coords['dye'] = grid.interp(ds_full.coords['dyt'], 'X')

    ds_full.coords['area_t'] = ds_full.coords['dxt'] * ds_full.coords['dyt']
    ds_full.coords['area_e'] = ds_full.coords['dxe'] * ds_full.coords['dye']
    ds_full.coords['area_ne'] = ds_full.coords['dxne'] * ds_full.coords['dyne']
    ds_full.coords['area_n'] = ds_full.coords['dxn'] * ds_full.coords['dyn']

    # should i return the coords to dask?
    return ds_full
    
# TODO: check if all these hardcoded lines work for each CMIP6 model...
def recreate_full_grid_old(ds, lon_bound_name='lon_bounds', lat_bound_name='lat_bounds'):
    ds_full = generate_grid_ds(ds, {'X':'x', 'Y':'y'}, position=('center','right'))
    
    grid = Grid(ds_full)
    
    # Derive distances at u point
    dlon = 0 # ill interpolate that later
    dlat = ds[lat_bound_name].isel(vertex=1) - ds[lat_bound_name].isel(vertex=2)
    # interpolate the centered position
    lat = (ds[lat_bound_name].isel(vertex=1) + ds[lat_bound_name].isel(vertex=2)) / 2
    lon = (ds[lon_bound_name].isel(vertex=1) + ds[lon_bound_name].isel(vertex=2)) / 2
    _, dy = dll_dist(dlon, dlat, lon, lat)
    # strip coords and rename dims
    ds_full.coords['dye'] = (['y', 'x_right'], dy.data)
    
    # Derive distances at v point
    dlon = ds[lon_bound_name].isel(vertex=3) - ds[lon_bound_name].isel(vertex=2)
    
    print(dlon)
    # special check for dlon
    temp = dlon.load().data
    temp[temp < 0] = temp[temp < 0] + 360.0
    dlon.data = temp
    
    dlat = 0  # ill interpolate that later
    # interpolate the centered position
    lat = (ds[lat_bound_name].isel(vertex=3) + ds[lat_bound_name].isel(vertex=2)) / 2
    lon = (ds[lon_bound_name].isel(vertex=3) + ds[lon_bound_name].isel(vertex=2)) / 2
    dx, _ = dll_dist(dlon, dlat, lon, lat)
    # strip coords and rename dims
    ds_full.coords['dxn'] = (['y_right', 'x'], dx.data)
    
    # interpolate the missing metrics
    ds_full.coords['dxt'] = grid.interp(ds_full.coords['dxn'], 'Y')
    ds_full.coords['dxne'] = grid.interp(ds_full.coords['dxn'], 'X')
    
    ds_full.coords['dyt'] = grid.interp(ds_full.coords['dye'], 'X')
    ds_full.coords['dyne'] = grid.interp(ds_full.coords['dye'], 'Y')
    
    return ds_full

def dll_dist_old(dlon, dlat, lon, lat):
        """Converts lat/lon differentials into distances in meters

        PARAMETERS
        ----------
        dlon : xarray.DataArray longitude differentials
        dlat : xarray.DataArray latitude differentials
        lon  : xarray.DataArray longitude values
        lat  : xarray.DataArray latitude values

        RETURNS
        -------
        dx  : xarray.DataArray distance inferred from dlon
        dy  : xarray.DataArray distance inferred from dlat
        """

        distance_1deg_equator = 111000.0
        dx = dlon * xr.ufuncs.cos(xr.ufuncs.deg2rad(lat)) * distance_1deg_equator
        dy = ((lon * 0) + 1) * dlat * distance_1deg_equator
        return dx, dy

