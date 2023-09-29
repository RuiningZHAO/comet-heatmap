"""
A tool for finding all observable comets in a given semester at a given site.
"""

from astropy.time import Time
import astropy.config as _config

__version__ = '0.0.1.0'


class Conf(_config.ConfigNameSpace):


    path_to_library = _config.ConfigItem(
        os.path.join(os.path.split(__file__)[0], 'lib'), cfgtype='string', 
        description='Path to library.'
    )

    location = _config.ConfigItem(
        '327', cfgtype='string', 
        description=(
            'MPC observatory code.'
        )
    )

    timezone = _config.ConfigItem(
        8, cfgtype='integer', 
        description=(
            'Timezone of the observatory.'
        )
    )

    minimum_elevation = _config.ConfigItem(
        20., cfgtype='float', 
        description=(
            'Minimum elevation of targets.'
        )
    )

    limiting_magnitude = _config.ConfigItem(
        20., cfgtype='float', 
        description=(
            'Limiting magnitude.'
        )
    )

    database = _config.ConfigItem(
        'jpl', cfgtype='option(jpl, mpc)', 
        description=(
            'Database.'
        )
    )

    start = _config.ConfigItem(
        Time.now().to_value('iso', 'date'), cfgtype='string', 
        description=(
            'Starting date.'
        )
    )

    end =  _config.ConfigItem(
        Time(Time.now().mjd + 183, format='mjd').to_value('iso', 'date'), 
        cfgtype='string', 
        description=(
            'Ending date.'
        )
    )


conf = Conf()

del Time, _config