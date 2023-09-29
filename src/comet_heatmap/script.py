import time, argparse

# AstroPy
from astropy.time import Time

from . import conf
from .catalog import catalog
from .ephemerides import ephemerides
from .observability import observability
from .plot import plot


def run():
    """
    """

    parser = argparse.ArgumentParser()
    parser.add_argument(
        '-o', '--output_dir', default='', type=str, 
        help='Output (saving) directory.'
    )
    parser.add_argument(
        '-v', '--verbose', action='store_true', 
        help='Verbose or not.'
    )
    
    # Parse
    args = parser.parse_args()
    save_dir = os.path.abspath(args.output_dir)
    verbose = args.verbose

    # Timing
    t0 = time.time()

    # Update catalog
    if verbose:
        print('Updating comet catalog...')
    catalog()
    if verbose:
        print(f'{time.time() - t0} (total: {time.time() - t0})')

    # Timing
    t1 = time.time()

    # Calculate ephemerides
    if verbose:
        print('Calculating comet ephemerides...')
    ephemerides()
    if verbose:
        print(f'{time.time() - t1} (total: {time.time() - t0})')

    # Timing
    t2 = time.time()

    # Calculate observability
    if verbose:
        print('Calculating comet observability...')
    observability()
    if verbose:
        print(f'{time.time() - t2} (total: {time.time() - t0})')

    # Timing
    t3 = time.time()

    # Plot
    if verbose:
        print('Plotting...')
    plot(save_dir)
    if verbose:
        print(f'{time.time() - t3} (total: {time.time() - t0})')


if __name__ == "__main__":
    main()
