from astroplan import FixedTarget
from astropy.coordinates import SkyCoord, EarthLocation, get_body
from astropy.time import Time
import astropy.units as u
from astropy.visualization import time_support, quantity_support
import numpy as np

from libballoonplanner import (get_observer, observability, time_vs_altitude,
                               time_vs_sun_relative_az, save_config)

time_support()
quantity_support()


# Mission-specific properties
EL_MIN = 20 * u.deg
EL_MAX = 52 * u.deg
DAZ_MIN = 90 * u.deg
DAZ_MAX = 225 * u.deg

LDB = (-77.861 * u.deg, 167.061 * u.deg)
FLOAT_ALT = 37000 * u.m


if __name__ == '__main__':
    # Example code using library: TIM-like mission profile and target list

    # save off mission-specific constraints for GUI usage
    param_dict = dict(
        el_min=EL_MIN.to(u.deg).value,
        el_max=EL_MAX.to(u.deg).value,
        daz_min=DAZ_MIN.to(u.deg).value,
        daz_max=DAZ_MAX.to(u.deg).value,
    )
    save_config('tim_constraints.json', param_dict)

    launch_location = EarthLocation(lat=LDB[0], lon=LDB[1], height=FLOAT_ALT)
    launch_time = Time('2027-12-25 00:00:00', scale='utc', location=launch_location)

    step_hr = 1.0 * u.hr
    my_duration = 1 * 24 * u.hr
    timespan = np.arange(0, my_duration.value + step_hr.value, step_hr.value) * u.hr
    times = launch_time + timespan

    tim = get_observer(
        launch_location.lat,
        launch_location.lon,
        launch_location.height,
        times,
        stationary=False,
        name='TIM'
    )

    # Science targets
    target_names = ['RCW 38', 'RCW 36', 'RCW 19', 'Vy CMa']
    targets = [FixedTarget.from_name(target_name) for target_name in target_names]
    goods_s = FixedTarget(
        SkyCoord(
            ra='3h32m36.51s',
            dec='-27d47m33.74s',
            obstime=times,
            location=tim.location
        ),
        name='GOODS-S'
    )
    targets += [goods_s,]

    table = observability(
        targets,
        tim,
        times,
        el_min=EL_MIN,
        el_max=EL_MAX,
        daz_min=DAZ_MIN,
        daz_max=DAZ_MAX,
        plot=False
    )
    print(table)

    fig, ax = time_vs_altitude(targets, tim, times, el_min=EL_MIN, el_max=EL_MAX)
    ax.set_title('Elevation Axis Limits')

    fig, ax = time_vs_sun_relative_az(targets, tim, times, daz_min=DAZ_MIN, daz_max=DAZ_MAX)
    ax.set_title('Sun-Relative Azimuth Angle Limits')


    # Planets
    target_names = ['mars', 'jupiter', 'saturn', 'uranus', 'neptune']
    targets = []
    for name in target_names:
        foo = get_body(name, times, location=tim.location)
        # notimplemented: targets with more than one ra/dec vs. time (planets)
        coord = SkyCoord(ra=foo.ra, dec=foo.dec, obstime=times,
            location=tim.location)
        coord.name = name
        targets.append(coord)

    # Plotting objects with time-varying sky coordinates depends on
    # https://github.com/astropy/astroplan/pull/634, but tabular works ok
    table = observability(
        targets,
        tim,
        times,
        grid_times_targets=False,
        el_min=EL_MIN,
        el_max=EL_MAX,
        daz_min=DAZ_MIN,
        daz_max=DAZ_MAX,
        plot=False
    )
    print(table)

    fig, ax = time_vs_altitude(targets, tim, times, el_min=EL_MIN, el_max=EL_MAX)
    ax.set_title('Elevation Axis Limits')

    fig, ax = time_vs_sun_relative_az(targets, tim, times, daz_min=DAZ_MIN, daz_max=DAZ_MAX)
    ax.set_title('Sun-Relative Azimuth Angle Limits')
