# NumPy
import numpy as np
# AstroPy
import astropy.units as u
from astropy.table import Table
from astropy.time import Time, TimeDelta
# sbplan
from sbplan.ephemerides import CometEphemerides
from sbplan.ephemerides.utils import siteLocation

from . import conf
from .utils import loadCometCatalog, silentRemove


def ephemerides(location=conf.location, timezone=conf.timezone, start=conf.start, 
                end=conf.end, path_to_library=conf.path_to_library, 
                database=conf.database):
    """
    """

    utc_offset = TimeDelta(timezone * u.Unit("h"))

    # Convert local mid-night to UTC
    mjd_start = (
        Time(f"{start} 00:00:00.000") + TimeDelta(1 * u.Unit("d")) - utc_offset).mjd
    mjd_end = (
        Time(f"{end} 00:00:00.000") + TimeDelta(1 * u.Unit("d")) - utc_offset).mjd

    # Time tag
    time_tag = Time(np.arange(mjd_start, mjd_end + 1, 1), format='mjd')

    location = siteLocation(location)
    
    catalog = loadCometCatalog(database)
    
    params = ["RA", "DEC", "delta", "r", "alpha", "lunar_elong", "Tmag"]

    comet_ephemerides = CometEphemerides(time_tag=time_tag, location=location)
    _, ephemerides = comet_ephemerides.get(catalog=catalog, params=params)

    # Write
    silentRemove(os.path.join(path_to_library, "ephemerides.hdf5"))
    Table({"time_tag": time_tag}).write(
        os.path.join(path_to_library, "ephemerides.hdf5"), path="Datetime", 
        serialize_meta=True)
    for pdes, tbl in ephemerides.items():
        tbl.write(
            os.path.join(path_to_library, "ephemerides.hdf5"), path=pdes, append=True, 
            serialize_meta=True)


def debug():
    """
    """

    ephemerides()

    return None


if __name__ == "__main__":
    debug()