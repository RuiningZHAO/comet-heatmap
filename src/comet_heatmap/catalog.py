import os

# AstroPy
from astropy.table import Table
# skyfield
from skyfield.api import Loader
from skyfield.data import mpc
# sbplan
from sbplan.catalog import CometCatalog
from sbplan.catalog.utils import toSkyfieldFormat

from . import conf


def catalog(path_to_library=conf.path_to_library, database=conf.database):
    """
    """
    
    if database == "mpc":

        # - Loader
        load = Loader(path_to_library)

        # - Download MPC catalog of observable comets
        with load.open(mpc.COMET_URL, reload=True) as f:
            catalog = mpc.load_comets_dataframe(f)

        catalog = Table.from_pandas(catalog)
        toSkyfieldFormat(catalog, database=database)

    elif database == "jpl":

        catalog = CometCatalog.get(cat_type="obs", update=True)

        catalog = toSkyfieldFormat(catalog, database=database)

    catalog.write(
        os.path.join(path_to_library, f"{database}.csv"), overwrite=True)

    return catalog


def debug():
    """
    """

    catalog()
    
    return None


if __name__ == "__main__":
    debug()