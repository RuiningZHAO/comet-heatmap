import h5py

# NumPy
import numpy as np
# AstroPy
import astropy.units as u
from astropy.table import Table
from astropy.coordinates import SkyCoord
# Astroplan
from astroplan import Observer
# sbplan
from sbplan.ephemerides.utils import siteLocation

from . import conf
from utils import loadCometCatalog, silentRemove


def observability(location=conf.location, horizon=conf.minimum_elevation, 
                  mag_lim=conf.limiting_magnitude, 
                  path_to_library=conf.path_to_library, database=conf.database):
    """
    """

    # Load time tag
    time_tag = Table.read(
        os.path.join(path_to_library, "ephemerides.hdf5"), path="Datetime")["time_tag"]

    # Load comet ephemerides
    catalog = loadCometCatalog(database)
    pdes_arr = catalog["designation"].value
    RA_mat   = np.zeros((time_tag.shape[0], pdes_arr.shape[0]))
    DEC_mat  = np.zeros((time_tag.shape[0], pdes_arr.shape[0]))
    Tmag_mat = np.zeros((time_tag.shape[0], pdes_arr.shape[0]))
    for i, pdes in enumerate(pdes_arr):
        tbl = Table.read(os.path.join(path_to_library, "ephemerides.hdf5"), path=pdes)
        RA_mat[:, i] = tbl["RA"].value
        DEC_mat[:, i] = tbl["DEC"].value
        Tmag_mat[:, i] = tbl["Tmag"].value
    sky_coord_mat = SkyCoord(RA_mat, DEC_mat, unit="deg")

    # Location
    location = siteLocation(location)
    tianma = Observer(location=location, name="Tianma")
    
    # Fraction of lunar illumination (FLI)
    lunar_illum = tianma.moon_illumination(time_tag)
    
    # Twilight
    t_twilight_evening = tianma.twilight_evening_astronomical(
        time_tag, which="previous", n_grid_points=150)
    t_twilight_morning = tianma.twilight_morning_astronomical(
        time_tag, which="next", n_grid_points=150)
    
    t_set  = np.ma.zeros((time_tag.shape[0], pdes_arr.shape[0]))
    t_rise = np.ma.zeros((time_tag.shape[0], pdes_arr.shape[0]))
    is_up  = np.zeros((time_tag.shape[0], pdes_arr.shape[0]), dtype=bool)
    for i, t_even in enumerate(t_twilight_evening):
        t_set[i] = tianma.target_set_time(
            time=t_even, target=sky_coord_mat[i], which="next", 
            horizon=(horizon * u.Unit("deg")), n_grid_points=150).jd
        t_rise[i] = tianma.target_rise_time(
            time=t_even, target=sky_coord_mat[i], which="next", 
            horizon=(horizon * u.Unit("deg")), n_grid_points=150).jd
        is_up[i] = tianma.target_is_up(
            time=t_even, target=sky_coord_mat[i], horizon=(horizon * u.Unit("deg")))
    # Case 0: always up
    is_high_0 = t_set.mask & is_up
    # Case 1: set before morning twilight
    is_high_1 = t_set < t_twilight_morning.jd[:, np.newaxis]
    # Case 2: set after morning twilight & not rise between morning twilight and set
    is_high_2 = (
        (t_set > t_twilight_morning.jd[:, np.newaxis]) 
        & ~((t_twilight_morning.jd[:, np.newaxis] < t_rise) & (t_rise < t_set))
    )
    is_high = is_high_0 | is_high_1.filled(False) | is_high_2.filled(False)
    
    # Filter
    mask = (np.ma.array(Tmag_mat, mask=~is_high).min(axis=0) > mag_lim).filled(True)
    catalog = catalog[~mask]
    colnames = Table.read(
        os.path.join(path_to_library, "ephemerides.hdf5"), 
        path=catalog["designation"].value[0]).colnames
    matrices = {
        "is_high": is_high[:, ~mask],
    }
    for colname in colnames:
        matrices[colname] = np.zeros((time_tag.shape[0], (~mask).sum()))
    for i, pdes in enumerate(catalog["designation"].value):
        tbl = Table.read(
            os.path.join(path_to_library, "ephemerides.hdf5"), path=pdes)
        for colname in colnames:
            matrices[colname][:, i] = tbl[colname].value

    # Write
    silentRemove(os.path.join(path_to_library, "matrices.hdf5"))
    with h5py.File(os.path.join(path_to_library, "matrices.hdf5"), "w") as hf:
        # Catalog
        catalog.write(hf, path="catalog", serialize_meta=True)
        # Fraction of lunar illumination (FLI)
        Table(data={"time_tag": time_tag, "lunar_illum": lunar_illum,}).write(
            hf, path="FLI", append=True, serialize_meta=True)
        # Matrices
        for key, val in matrices.items():
            hf.create_dataset(key, data=val)


def debug():
    """
    """

    observability(horizon=25, mag_lim=20)

    return None


if __name__ == "__main__":
    debug()