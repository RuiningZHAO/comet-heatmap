import os, errno
# AstroPy
from astropy.table import Table

def loadCometCatalog(database):
    """
    """

    catalog = Table.read(f"lib/{database}.csv")
    
    return catalog

def silentRemove(path):
    """
    """

    try:
        os.remove(path)
    except OSError as e:
        if e.errno != errno.ENOENT: # errno.ENOENT = no such file or directory
            raise # re-raise exception if a different error occurred