from astroplan import Observer, Constraint, AltitudeConstraint
from astroplan.plots import plot_altitude, plot_airmass

from astropy.coordinates import SkyCoord, get_body
import astropy.units as u
from astropy.visualization import time_support, quantity_support
import json

import matplotlib.pyplot as plt
import numpy as np

from .astroplan_like import my_observability_table

time_support()
quantity_support()


EL_MIN_DEFAULT = 0. * u.deg
EL_MAX_DEFAULT = 90. * u.deg
DAZ_MIN_DEFAULT = 0. * u.deg
DAZ_MAX_DEFAULT = 360. * u.deg
ALT_DEFAULT = 37000. * u.m
LOC_DEFAULT = (-77.861 * u.deg, 167.061 * u.deg)


def wrap360(ang):
    return ((ang + 360.) % 360.)


class SunRelativeAzConstraint(Constraint):
    '''
    Constrain an asymmetrical sun-relative azimuth pointing of the boresight to
    target.
    These angles are the angles that the boresight must point away from the sun
    in order to maintain the instrument and electronics in the shade of the
    sunshade. As always, azimuth is measured +clockwise from the reference
    point, in this case the sun azimuth. The range is [0, 360] deg.
    '''
    def __init__(self, min=None, max=None):
        '''
        min : `~astropy.units.Quantity` or `None` (optional)
            Minimum acceptable azimuth angle between sun and boresight. `None`
            indicates no limit.
        max : `~astropy.units.Quantity` or `None` (optional)
            Maximum acceptable azimuth angle between sun and boresight. `None`
            indicates no limit.
        '''
        self.min = min if min is not None else 0*u.deg
        self.max = max if max is not None else 360*u.deg


    def compute_constraint(self, times, observer, targets):
        sun_body = get_body('sun', times, location=observer.location)
        sun = SkyCoord(ra=sun_body.ra, dec=sun_body.dec, obstime=times,
                       location=observer.location)
        sun_altaz = sun.transform_to('altaz')

        # Calculate separation between boresight and sun
        # Targets are automatically converted to SkyCoord objects
        # by __call__ before compute_constraint is called.

        try:
            n_targets = len(targets)
        except TypeError as e:
            targets = SkyCoord([targets,])
            n_targets = len(targets)
        try:
            n_times = len(times)
        except TypeError as e:
            n_times = len(times)

        delta_az = np.atleast_2d(np.empty((n_targets, n_times)) * u.deg)
        for i, target in enumerate(targets):
            # TODO: not sure if squeeze is necessary
            target_with_loc = SkyCoord(
                ra=np.squeeze(target.ra),
                dec=np.squeeze(target.dec),
                obstime=np.squeeze(times),
                location=np.squeeze(observer.location)
            )
            target_altaz = target_with_loc.transform_to('altaz')
            delta_az[i,:] = wrap360(target_altaz.az.deg - sun_altaz.az.deg) * u.deg # time axis

        if self.min is None and self.max is not None:
            mask = self.max >= delta_az
        elif self.max is None and self.min is not None:
            mask = self.min <= delta_az
        elif self.min is not None and self.max is not None:
            mask = ((self.min <= delta_az) &
                    (delta_az <= self.max))
        else:
            raise ValueError("No max and/or min specified in "
                             "SunRelativeAzConstraint.")

        return mask


def get_observer(
        launch_lat,
        launch_lon,
        float_alt,
        times,
        stationary=False,
        name='observer'
    ):
    '''
    Generate an observer object with lat/lon that change over time like a
    long duration Antarctic balloon.

    This is not a high fidelity model, but is based loosely on parameters
    observed from ground track of a previous LDB flight.

    All parameters can be 1D arrays, and should match in shape.

    Parameters
    ----------
    launch_lat : astropy.Quantity
        initial latitude position of Observer
    launch_lon : astropy.Quantity
        initial longitude position of Observer
    float_alt : astropy.Quantity
        initial altitude (HAE) position of Observer
    times : astropy.Quantity
        times over which the observation campaign takes place
    stationary : bool (optional)
        If False, applies a simplified constant altitude model of an Antarctic 
        LDB orbit
    name : str
        Name for Observer. Typically this is a site or observing facility.

    Returns
    -------
    observer : astroplan.Observer
    '''
    t = (times - times[0]).to(u.hr)
    ldb = (launch_lat, launch_lon)
    if not stationary:
        # ground track:
        # Salter Test Flight Universal completed ~1 circuit in 11 days, 6 hr, 57 min.
        # Average the longitudinal velocity. Really, ballons have a spatial velocity,
        # so this constant lat assumption is an oversimplification.
        dlon_dt = -(360 * u.deg / (11 * u.day + 6 * u.hr + 57 * u.min)).to(u.deg / u.hr)
        lon = ldb[1] + dlon_dt * t
        lon = wrap360(lon.to(u.deg).value) * u.deg
        # add a little lat wobble, eyeballed from STFU flight track
        lat = 0 * u.deg + np.ones_like(lon.value) * ldb[0] + .2 * u.deg * np.sin(lon.to(u.rad) * 15)
        observer = Observer(
            longitude=lon,
            latitude=lat,
            elevation=float_alt,
            name=name
        )
    else:
        observer = Observer(
            longitude=launch_lon*np.ones_like(t.value),
            latitude=launch_lat*np.ones_like(t.value),
            elevation=float_alt,
            name=name
        )

    return observer


def observability(
        targets:list,
        observer:Observer,
        times,
        el_min=EL_MIN_DEFAULT,
        el_max=EL_MAX_DEFAULT,
        daz_min=DAZ_MIN_DEFAULT,
        daz_max=DAZ_MAX_DEFAULT,
        grid_times_targets=True,
        plot=True
    ):

    constraints = [
        AltitudeConstraint(el_min, el_max),
        SunRelativeAzConstraint(min=daz_min, max=daz_max)
    ]
    table = my_observability_table(constraints, observer, targets, times, grid_times_targets=grid_times_targets)

    if plot:
        # https://astroplan.readthedocs.io/en/latest/tutorials/constraints.html
        for j, target in enumerate(targets):
            observability_grid = np.zeros((len(constraints), len(times)))
            for i, constraint in enumerate(constraints):
                # Evaluate each constraint
                observability_grid[i, :] = constraint(
                    observer,
                    target,
                    times=times,
                    grid_times_targets=grid_times_targets
                )

            # Create plot showing observability of the target:
            extent = [-0.5, -0.5+len(times), -0.5, len(constraints) - 0.5]

            fig, ax = plt.subplots(figsize=(10,4))
            ax.imshow(observability_grid, extent=extent, cmap='bone_r', vmin=0, vmax=1, origin='lower')

            ax.set_yticks(range(0, len(constraints)))
            ax.set_yticklabels([c.__class__.__name__ for c in constraints])

            ax.set_xticks(range(len(times)))
            ax.set_xticklabels([t.datetime.strftime("%H:%M") for t in times])

            ax.set_xticks(np.arange(extent[0], extent[1]), minor=True)
            ax.set_yticks(np.arange(extent[2], extent[3]), minor=True)

            ax.grid(which='minor', color='w', linestyle='-', linewidth=1)
            ax.tick_params(axis='x', which='minor', bottom='off')
            plt.setp(ax.get_xticklabels(), rotation=30, ha='right')

            ax.tick_params(axis='y', which='minor', left='off')
            ax.set_xlabel('Time on {0} (UTC)'.format(times[0].datetime.date()))
            fig.subplots_adjust(left=0.25, right=0.9, top=0.9, bottom=0.1)
            ax.set_title(f'{target.name}: {table["fraction of time observable"][j] * (times[-1] - times[0]).to(u.hr):.2f}' + 
                         '\nBlack = Observable')
            fig.tight_layout()
        plt.show()

    return table


def time_vs_altitude(
        targets:list,
        observer:Observer,
        times,
        el_min=EL_MIN_DEFAULT,
        el_max=EL_MAX_DEFAULT,
    ):
    fig, ax = plt.subplots(figsize=(12,4))
    plot_altitude(targets, observer, times, ax=ax)
    ax.axhline(el_min, color='limegreen')
    ax.axhline(el_max, color='limegreen')
    ax.axhspan(el_min, el_max, color='limegreen', alpha=0.3)
    ax.set_ylabel('Altitude (deg)')
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    plt.show()

    return fig, ax


def time_vs_airmass(
        targets:list,
        observer:Observer,
        times
    ):
    fig, ax = plt.subplots(figsize=(12,4))
    plot_airmass(targets, observer, times, ax=ax, brightness_shading=False)
    ax.legend()
    ax.grid(True)
    fig.tight_layout()
    plt.show()

    return fig, ax


def time_vs_sun_relative_az(
        targets:list,
        observer:Observer,
        times,
        daz_min=DAZ_MIN_DEFAULT,
        daz_max=DAZ_MAX_DEFAULT,
    ):
    '''
    Plot the sun-relative azimuth of the observer when observing the target.
    This is a delta azimuth, target az - sun az. Positive angles denote the
    observer is pointed to a greater azimuth angle than the sun,
    clockwise viewed from above.
    '''
    fig, ax = plt.subplots(figsize=(12,4))
    for target in targets:
        ax.plot(
            times,
            wrap360(
                observer.altaz(times, target=target).az.deg -
                observer.sun_altaz(times).az.deg
            ),
            marker='.',
            label=target.name
        )
    ax.axhline(daz_min, color='limegreen')
    ax.axhline(daz_max, color='limegreen')
    ax.axhspan(daz_min, daz_max, color='limegreen', alpha=0.3)
    ax.set_ylabel('Sun-Relative Azimuth\n$AZ_{sun} - AZ_{tgt}$ (deg)\n+clockwise from local North')
    ax.grid(True)
    ax.legend()
    fig.tight_layout()
    plt.show()

    return fig, ax


def ground_track(
        observer,
        times
    ):
    '''Plot the observer lat/lon over time, polar plot. Axes limits are assumed
    for an Antarctic LDB flight.'''
    this_lon = observer.longitude
    this_lat = observer.latitude
    t = (times - times[0]).to(u.hr)
    fig, ax = plt.subplots(figsize=(4, 4), subplot_kw=dict(projection='polar'))
    mappable = ax.scatter(this_lon.to(u.rad), this_lat, c=range(len(this_lon)))
    ax.set_rmin(-90)
    ax.set_rmax(-70)
    ax.set_theta_offset(np.pi / 2)
    ax.set_theta_direction(-1)
    ax.set_thetamin(0)
    ax.set_thetamax(360)
    ax.grid(True)
    ax.set_title(f'Looking Down:\nGround Track: {(t[-1] - t[0]).to(u.day).value:.1f} days')
    fig.colorbar(mappable, ax=ax, label='Time Since Launch (hr)')
    fig.tight_layout()
    ax.set_facecolor('lightgrey')
    fig.set_facecolor('lightgrey')
    plt.show()

    return fig, ax


def load_config(filename):
    with open(filename, 'r') as f:
        param_dict = json.load(f)
    el_min = param_dict['el_min'] * u.deg
    el_max = param_dict['el_max'] * u.deg
    daz_min = param_dict['daz_min'] * u.deg
    daz_max = param_dict['daz_max'] * u.deg

    return el_min, el_max, daz_min, daz_max


def save_config(filename, param_dict):
    with open(filename, 'w') as f:
        json.dump(param_dict, f)