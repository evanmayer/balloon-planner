from astroplan import FixedTarget

from astropy.coordinates import Angle, SkyCoord, EarthLocation, get_body
from astropy.coordinates.name_resolve import NameResolveError
from astropy.time import Time
import astropy.units as u

import matplotlib.pyplot as plt
import numpy as np
import sys

from tkinter import (Tk, DoubleVar, BooleanVar, StringVar, N, E, S, W)
from tkinter import font
from tkinter import ttk

from balloonplanner.libballoonplanner import (get_observer, observability, time_vs_altitude,
                               time_vs_airmass, time_vs_sun_relative_az,
                               ground_track, load_config)
from balloonplanner.libballoonplanner import (EL_MIN_DEFAULT, EL_MAX_DEFAULT, DAZ_MIN_DEFAULT,
                               DAZ_MAX_DEFAULT, LOC_DEFAULT, ALT_DEFAULT)


def dispatch_analysis():
    '''Do input checking/setup, then run the selected function.'''
    # Gather input
    my_lat = lat.get() * u.deg
    my_lon = lon.get() * u.deg
    my_alt = alt.get() * u.m

    my_launch_date = launch_date.get()
    my_duration = duration.get() * u.hr
    my_step_hr = dt.get() * u.hr

    my_resolver = resolver.get()
    my_method = method_var.get()

    launch_location = EarthLocation(lat=my_lat, lon=my_lon, height=my_alt)
    launch_time = Time(my_launch_date, scale='utc', location=launch_location)
    timespan = np.arange(0, my_duration.value + my_step_hr.value, my_step_hr.value) * u.hr
    times = launch_time + timespan

    my_observer = get_observer(
        launch_location.lat,
        launch_location.lon,
        launch_location.height,
        times,
        stationary=stationaryFlag.get()
    )

    if 'name' in my_resolver:
        my_label = tgt_name.get()
        is_solar_system = False
        try:
            foo = SkyCoord.from_name(tgt_name.get())
            coord = SkyCoord(ra=foo.ra, dec=foo.dec, obstime=times, location=my_observer.location)
        except NameResolveError as e:
            print(tgt_name.get(), 'not found in CDS database. Trying solar system bodies...')
            foo = get_body(tgt_name.get(), times, location=my_observer.location)
            coord = SkyCoord(ra=foo.ra, dec=foo.dec, obstime=times,
                       location=my_observer.location)
            is_solar_system = True
    else:
        my_label = tgt_label.get()
        ra_str, dec_str = tgt_radec.get().split(', ')
        # If not sexagesimal, assume decimal degrees. Let Angle and SkyCoord 
        # handle input sanitizing
        if ':' not in ra_str:
            ra_str += 'd'
        if ':' not in dec_str:
            dec_str += 'd'
        ra = Angle(ra_str)
        dec = Angle(dec_str)
        coord = SkyCoord(ra, dec, obstime=times, location=my_observer.location)
    coord.name = my_label
    # print(f'{coord.name}: {coord}')

    if is_solar_system:
        targets = [coord,]
    else:
        targets = [FixedTarget(coord, my_label),]

    # Dispatch call
    if 'constraints' in my_method:
        if (is_solar_system):
            table = observability(
                targets,
                my_observer,
                times,
                el_min=el_min,
                el_max=el_max,
                daz_min=daz_min,
                daz_max=daz_max,
                grid_times_targets=(not is_solar_system),
                plot=False
            )
            print(table)
            err_msg = ('Plotting objects with time-varying sky coordinates ' +
                'depends on https://github.com/astropy/astroplan/pull/634')
            raise NotImplementedError(err_msg)
        table = observability(
            targets,
            my_observer,
            times,
            el_min=el_min,
            el_max=el_max,
            daz_min=daz_min,
            daz_max=daz_max,
            grid_times_targets=(not is_solar_system),
        )
        print(table)
    elif 'time_vs_altitude' in my_method:
        _, _ = time_vs_altitude(
            targets,
            my_observer,
            times,
            el_min=el_min,
            el_max=el_max
        )
    elif 'time_vs_airmass' in my_method:
        _, _ = time_vs_airmass(targets, my_observer, times)
    elif 'time_vs_sun_relative_az' in my_method:
        _, _ = time_vs_sun_relative_az(
            targets,
            my_observer,
            times,
            daz_min=daz_min,
            daz_max=daz_max
        )
    elif 'ground_track' in my_method:
        _, _ = ground_track(my_observer, times)
    else:
        pass


def cleanup():
    plt.close('all')
    root.destroy()


if __name__ == '__main__':
    # check for user-specified config json to tell us constraints info
    if len(sys.argv) > 1:
        config_file = sys.argv[1]
        el_min, el_max, daz_min, daz_max = load_config(config_file)
    else:
        el_min = EL_MIN_DEFAULT
        el_max = EL_MAX_DEFAULT
        daz_min = DAZ_MIN_DEFAULT
        daz_max = DAZ_MAX_DEFAULT

    # Build GUI
    root = Tk()

    style = ttk.Style(root)
    style.theme_use('classic')

    default_font = font.nametofont("TkDefaultFont")
    bold_font = font.nametofont("TkDefaultFont").copy()
    default_font.configure(size=16)
    bold_font.configure(weight='bold', size=16)

    root.title("evanmayer's Balloon Mission Planner")
    mainframe = ttk.Frame(root, padding=(3,3,3,3))
    mainframe.grid(column=0, row=0, sticky=(N, W, E, S))

    # Launch location
    ttk.Label(mainframe, text='Location', font=bold_font).grid(column=0, row=0, sticky=W)
    ttk.Label(mainframe, text='Latitude (deg)').grid(column=0, row=1, sticky=W)
    lat = DoubleVar()
    lat_entry = ttk.Entry(mainframe, width=6, textvariable=lat)
    lat_entry.grid(column=0, row=2, sticky=(W, E))

    ttk.Label(mainframe, text='Longitude (deg)').grid(column=0, row=3, sticky=W)
    lon = DoubleVar()
    lon_entry = ttk.Entry(mainframe, width=6, textvariable=lon)
    lon_entry.grid(column=0, row=4, sticky=(W, E))

    ttk.Label(mainframe, text='Altitude (m)').grid(column=0, row=5, sticky=W)
    alt = DoubleVar()
    alt_entry = ttk.Entry(mainframe, width=6, textvariable=alt)
    alt_entry.grid(column=0, row=6, sticky=(W, E))

    stationaryFlag = BooleanVar(value=False)
    is_stationary = ttk.Checkbutton(mainframe, text="Stationary lat/lon/alt",
        variable=stationaryFlag)
    is_stationary.grid(column=0, row=7, sticky=(W,E))

    # Mission duration
    ttk.Label(mainframe, text='Mission Time', font=bold_font).grid(column=0, row=8, sticky=W)
    ttk.Label(mainframe, text='Launch Date (UTC)\n[YYYY-MM-DD HH:MM:SS]').grid(column=0, row=9, sticky=W)
    launch_date = StringVar()
    launch_date_entry = ttk.Entry(mainframe, width=19, textvariable=launch_date)
    launch_date_entry.grid(column=0, row=10, sticky=(W, E))

    ttk.Label(mainframe, text='Duration (hr)').grid(column=0, row=11, sticky=W)
    duration = DoubleVar()
    duration_entry = ttk.Entry(mainframe, width=3, textvariable=duration)
    duration_entry.grid(column=0, row=12, sticky=(W, E))

    ttk.Label(mainframe, text='Time step (hr)').grid(column=0, row=13, sticky=W)
    dt = DoubleVar()
    dt_entry = ttk.Entry(mainframe, width=3, textvariable=dt)
    dt_entry.grid(column=0, row=14, sticky=(W, E))


    # Target
    ttk.Label(mainframe, text='Target', font=bold_font).grid(column=1, row=0, sticky=W)

    # Either name resolution or RA/Dec query are allowed. Enable entry for the selected mode.
    tgt_name = StringVar()
    tgt_name_entry = ttk.Entry(mainframe, width=20, textvariable=tgt_name)
    tgt_name_entry.grid(column=2, row=1, sticky=(W, E))

    tgt_radec = StringVar()
    tgt_radec_entry = ttk.Entry(mainframe, width=20, textvariable=tgt_radec)
    tgt_radec_entry.grid(column=2, row=2, sticky=(W, E))

    # Target name for plot labels
    ttk.Label(mainframe, text='Plot Label').grid(column=1, row=3, sticky=W)
    tgt_label = StringVar()
    tgt_label_entry = ttk.Entry(mainframe, width=20, textvariable=tgt_label)
    tgt_label_entry.grid(column=2, row=3, sticky=(W, E))


    resolver = StringVar()

    def enable_resolver_field():
        if 'name' in resolver.get():
            tgt_name_entry.state(['!disabled',])
            tgt_radec_entry.state(['disabled',])
            tgt_label_entry.state(['disabled',])
        else:
            tgt_name_entry.state(['disabled',])
            tgt_radec_entry.state(['!disabled',])
            tgt_label_entry.state(['!disabled',])

    name_button = ttk.Radiobutton(mainframe, text='CDS Object Name', command=enable_resolver_field, variable=resolver, value='name')
    name_button.grid(sticky=W, column=1, row=1)
    radec_button = ttk.Radiobutton(mainframe, text='Object RA/Dec [RA,Dec]', command=enable_resolver_field, variable=resolver, value='radec')
    radec_button.grid(sticky=W, column=1, row=2)

    # Analysis options
    ttk.Label(mainframe, text='Analysis', font=bold_font).grid(column=1, row=7, sticky=W)
    method_var = StringVar()
    method = ttk.Combobox(mainframe, textvariable=method_var)
    method.state(['readonly'])
    method.grid(stick=(W, E), column=1, row=8)
    method['values'] = ('constraints', 'time_vs_altitude', 'time_vs_airmass', 'time_vs_sun_relative_az', 'ground_track')
    root.option_add("*TCombobox*Listbox*Font", default_font)

    style = ttk.Style()
    style.theme_use('alt')
    style.configure('My.TButton', font=default_font, background='#77ddff', foreground='#1111ff')
    go = ttk.Button(mainframe, text='Go', style='My.TButton', command=dispatch_analysis)
    go.grid(sticky=(W,E), row=14, column=2)

    exit_button = ttk.Button(mainframe, text='Exit', command=cleanup)
    exit_button.grid(sticky=(W,E), row=14, column=1)

    # Prepare for startup
    root.columnconfigure(0, weight=1)
    root.rowconfigure(0, weight=1)
    mainframe.columnconfigure(1, weight=1)
    for child in mainframe.winfo_children():
        child.grid_configure(padx=5, pady=5)
    entries = [lat_entry, lon_entry, alt_entry, launch_date_entry,
               duration_entry, dt_entry, tgt_name_entry, tgt_radec_entry,
               tgt_label_entry, method]
    for entry in entries:
        entry.config(font=default_font)
    lat_entry.focus()

    def set_defaults():
        lat.set(LOC_DEFAULT[0].value)
        lon.set(LOC_DEFAULT[1].value)
        alt.set(ALT_DEFAULT.value)
        launch_date.set('2027-12-25 00:00:00')
        duration.set(24.0)
        dt.set(1.0)
        name_button.invoke()
        tgt_name.set('RCW 38')
        method_var.set('constraints')

    set_defaults()

    root.bind('<Return>', dispatch_analysis)
    root.mainloop()