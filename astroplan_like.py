# astroplan BSD3 copyright notice, because my_observability_table is modified
# astroplan source code.
# This is temporary until objects with time-varying coordinates are better
# supported: https://github.com/astropy/astroplan/pull/634

# Copyright (c) 2015-2017, Astroplan Developers All rights reserved.

# Redistribution and use in source and binary forms, with or without
# modification, are permitted provided that the following conditions are met:

#     * Redistributions of source code must retain the above copyright notice,
#       this list of conditions and the following disclaimer.
#     * Redistributions in binary form must reproduce the above copyright
#       notice, this list of conditions and the following disclaimer in the
#       documentation and/or other materials provided with the distribution.
#     * Neither the name of the Astropy Team nor the names of its contributors
#       may be used to endorse or promote products derived from this software
#     * without specific prior written permission.

# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" 
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
# DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR 
# CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF 
# SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
# INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN 
# CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) 
# ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE 
# POSSIBILITY OF SUCH DAMAGE.


from astropy import table
from astropy import units as u
import numpy as np

from astroplan import time_grid_from_range


def my_observability_table(constraints, observer, targets, times=None,
                        time_range=None, time_grid_resolution=0.5*u.hour,
                        grid_times_targets=True):
    """
    ECM: Redefined from astroplan to allow ``grid_times_targets`` control. Objects
    with changing ra/dec (e.g. planets) cannot be handled without
    ``grid_times_targets=False``.

    Creates a table with information about observability for all  the ``targets``
    over the requested ``time_range``, given the constraints in
    ``constraints_list`` for ``observer``.

    Parameters
    ----------
    constraints : list or `~astroplan.constraints.Constraint`
        Observational constraint(s)

    observer : `~astroplan.Observer`
        The observer who has constraints ``constraints``

    targets : {list, `~astropy.coordinates.SkyCoord`, `~astroplan.FixedTarget`}
        Target or list of targets

    times : `~astropy.time.Time` (optional)
        Array of times on which to test the constraint

    time_range : `~astropy.time.Time` (optional)
        Lower and upper bounds on time sequence, with spacing
        ``time_resolution``. This will be passed as the first argument into
        `~astroplan.time_grid_from_range`. If a single (scalar) time, the table
        will be for a 24 hour period centered on that time.

    time_grid_resolution : `~astropy.units.Quantity` (optional)
        If ``time_range`` is specified, determine whether constraints are met
        between test times in ``time_range`` by checking constraint at
        linearly-spaced times separated by ``time_resolution``. Default is 0.5
        hours.

    Returns
    -------
    observability_table : `~astropy.table.Table`
        A Table containing the observability information for each of the
        ``targets``. The table contains four columns with information about the
        target and it's observability: ``'target name'``, ``'ever observable'``,
        ``'always observable'``, and ``'fraction of time observable'``. The
        column ``'time observable'`` will also be present if the ``time_range``
        is given as a scalar. It also contains metadata entries ``'times'``
        (with an array of all the times), ``'observer'`` (the
        `~astroplan.Observer` object), and ``'constraints'`` (containing the
        supplied ``constraints``).
    """
    if not hasattr(constraints, '__len__'):
        constraints = [constraints]

    is_24hr_table = False
    if hasattr(time_range, 'isscalar') and time_range.isscalar:
        time_range = (time_range-12*u.hour, time_range+12*u.hour)
        is_24hr_table = True

    applied_constraints = [constraint(observer, targets, times=times,
                                      time_range=time_range,
                                      time_grid_resolution=time_grid_resolution,
                                      grid_times_targets=grid_times_targets)
                           for constraint in constraints]
    constraint_arr = np.logical_and.reduce(applied_constraints)

    colnames = ['target name', 'ever observable', 'always observable',
                'fraction of time observable']

    target_names = [target.name for target in targets]
    ever_obs = np.any(constraint_arr, axis=1)
    always_obs = np.all(constraint_arr, axis=1)
    frac_obs = np.sum(constraint_arr, axis=1) / constraint_arr.shape[1]

    tab = table.Table(names=colnames, data=[target_names, ever_obs, always_obs,
                                            frac_obs])

    if times is None and time_range is not None:
        times = time_grid_from_range(time_range,
                                     time_resolution=time_grid_resolution)

    if is_24hr_table:
        tab['time observable'] = tab['fraction of time observable'] * 24*u.hour

    tab.meta['times'] = times.datetime
    tab.meta['observer'] = observer
    tab.meta['constraints'] = constraints

    return tab

