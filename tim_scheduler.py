from astroplan import FixedTarget, AirmassConstraint, AltitudeConstraint, ObservingBlock, Transitioner
from astroplan.scheduling import SequentialScheduler, PriorityScheduler
from astroplan.scheduling import Schedule
from astropy.coordinates import EarthLocation, SkyCoord
from astropy.time import Time
import astropy.units as u
from astropy.visualization import time_support, quantity_support
import numpy as np

from balloonplanner.libballoonplanner import (SunRelativeAzConstraint, get_observer, load_config, observability)

time_support()
quantity_support()

LDB = (-77.861 * u.deg, 167.061 * u.deg)
FLOAT_ALT = 37000 * u.m


from astroplan.scheduling import Scheduler, Scorer
from astroplan.utils import time_grid_from_range
from astroplan.constraints import AltitudeConstraint
from astropy import units as u

import numpy as np


if __name__ == '__main__':
    el_min, el_max, daz_min, daz_max = load_config('tim_constraints.json')

    launch_location = EarthLocation(lat=LDB[0], lon=LDB[1], height=FLOAT_ALT)
    launch_time = Time('2027-12-25 00:00:00', scale='utc', location=launch_location)

    step_hr = 1.0 * u.hr
    my_duration = 3 * 24 * u.hr
    timespan = np.arange(0, my_duration.value + step_hr.value, step_hr.value) * u.hr
    times = launch_time + timespan

    tim = get_observer(
        launch_location.lat,
        launch_location.lon,
        launch_location.height,
        times,
        stationary=True,
        name='TIM'
    )

    # Science targets
    target_names = ['RCW 38', 'RCW 36', 'RCW 19', 'Vy CMa']
    targets = [FixedTarget.from_name(target_name) for target_name in target_names]
    targets_dict = {target_names[i]: targets[i] for i in range(len(target_names))}
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

    airmass = AirmassConstraint(max=10, boolean_constraint=False)

    table = observability(
        targets,
        tim,
        times,
        el_min=el_min,
        el_max=el_max,
        daz_min=daz_min,
        daz_max=daz_max,
        # additional_constraints=[airmass,],
        plot=False
    )
    print(table)

    global_constraints = [
        AltitudeConstraint(el_min, el_max),
        # SunRelativeAzConstraint(min=daz_min, max=daz_max),
        # airmass
    ]

    goods_exptime = 1 * u.hr
    science_field = ObservingBlock(
        goods_s,
        goods_exptime,
        2,
        configuration={'fridge_state': 'cycleA'},
        name='GOODS-S'
    )

    calibrator_exptime = 0.1 * u.hr
    calibrator = ObservingBlock(
        targets_dict['Vy CMa'],
        calibrator_exptime,
        2,
        configuration={'fridge_state': 'cycleA'},
        name='Vy CMa'
    )

    opportunity_exptime = .25 * u.hr
    opportunity = ObservingBlock(
        targets_dict['RCW 19'],
        opportunity_exptime,
        99,
        configuration={'fridge_state': 'cycleA'},
        name='RCW 19'
    )

    cycle_exptime = 7 * u.hr
    fridge_cycle = ObservingBlock(
        targets_dict['RCW 19'],
        cycle_exptime,
        1,
        name='Fridge cycle'
    )

    slew_rate = 0.1*u.deg/u.second
    transitioner = Transitioner(
        slew_rate,
        {
            # model fridge cycles as transitions between two vacuous states:
            # every time a fridge cycle is required, 8 hours are lost.
            # Otherwise, target transitions shall be dictated only by slew rates
            'fridge_state':{
                ('cycleA', 'cycleB'): 8 * u.hr,
                ('cycleB', 'cycleA'): 8 * u.hr,
                'default': 1*u.s
            }
        }
    )

    blocks = []
    # attempt 18 hours of science per day, with calibrators in between
    blocks += [science_field, calibrator, opportunity] * 18
    # after a day of science, cycle the fridge
    # blocks += [fridge_cycle,]

    # prio_scheduler = PriorityScheduler(
    #     constraints=global_constraints,
    #     observer=tim,
    #     transitioner=transitioner
    # )
    # priority_schedule = Schedule(times[0], times[-1])
    # prio_scheduler(blocks, priority_schedule)
    seq_scheduler = SequentialScheduler(
        constraints=global_constraints,
        observer=tim,
        transitioner=transitioner
    )
    seq_schedule = Schedule(times[0], times[-1])
    seq_scheduler(blocks, seq_schedule)

    print(seq_schedule.to_table().pprint_all())

    from astroplan.plots import plot_schedule_airmass
    import matplotlib.pyplot as plt
    plt.figure(figsize = (14,6))
    plot_schedule_airmass(seq_schedule)
    plt.legend(loc = "upper right")
    plt.show()