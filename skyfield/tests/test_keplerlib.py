from numpy import pi, seterr, linspace
from skyfield.keplerlib import KeplerOrbit, propagate
from skyfield.elementslib import OsculatingElements
from skyfield.units import Angle, Distance, Velocity
from skyfield.tests.test_elementslib import compare, ele_to_vec
from skyfield.api import load
import os

seterr(all='raise')

# Test against HORIZONS.

def test_against_horizons_for_ceres():
    # See the following files in the Skyfield repository:
    #
    # horizons/ceres-orbital-elements
    # horizons/ceres-position

    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(2458886.500000000)

    a = 2.768873850275102E+00 # A
    e = 7.705857791518426E-02 # EC
    p_au = a * (1 - e ** 2)  # Wikipedia

    k = KeplerOrbit.from_mean_anomaly(
        p=Distance(au=p_au),  # see above
        e=e,
        i=Angle(degrees=2.718528770987308E+01),
        Om=Angle(degrees=2.336112629072238E+01),
        w=Angle(degrees=1.328964361683606E+02),
        M=Angle(degrees=1.382501360489816E+02),
        epoch=t,  #?
        mu_km_s=None,
        mu_au_d=2.9591220828559093E-04,
        center=None,
        target=None,
        center_name=None,
        target_name=None,
    )
    r, v = k._at(t)[:2]
    sun_au = [
        -0.004105894975783999, 0.006739680703224941, 0.002956344702049446,
    ]
    horizons_au = [
        1.334875927366032E+00, -2.239607658161781E+00, -1.328895183461897E+00,
    ]
    assert max(abs(r + sun_au - horizons_au)) < 2e-15

def test_propagating_from_vectors_to_vectors():
    # See the following file in the Skyfield repository:
    # horizons/hale-bopp-vector

    mu_m_s = 1.327124400189e20
    mu_km_s = mu_m_s * 1e-9
    mu = mu_km_s

    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(2454724.5)
    print(t.utc_jpl())

    from numpy import array
    r0 = array([1.777310651689592E+00, 1.638390146876578E+00, -2.712743223120575E+01])
    v0 = array([4.707733989610805E-04, -5.688697324947830E-04, -4.422633506777067E-03])
    from skyfield.units import Distance, Velocity
    orbit = KeplerOrbit(Distance(au=r0),
                        Velocity(au_per_d=v0),
                        t, mu_km_s=mu)

    # Propagating to the epoch itself should simply return the original
    # position and velocity vectors nearly unchanged.

    r, v = orbit._at(t)[:2]
    assert(max(abs(r0 - r)) < 1e-18)
    assert(max(abs(v0 - v)) < 1e-18)

    # Propagating to the date of perihelion should produce the vector
    # provided by HORIZONS.

    t = ts.tdb_jd(2450538.437848276)
    print(t.utc_jpl())
    r, v = orbit._at(t)[:2]
    # want_r
    # want_v
    print(r)
    print(v)

    asdf
    return

    t = ts.tdb_jd(2450538.4378482755)
    print(t.utc_jpl())

    e = 0.9949607008417696
    qr = 0.9174143409263262  # Perihelion distance
    semimajor_axis_au = qr / (1.0 - e)

    a = semimajor_axis_au
    p_au = a * (1 - e ** 2)  # Wikipedia

    k = KeplerOrbit.from_mean_anomaly(
        p=Distance(au=p_au),  # see above
        e=e,
        i=Angle(degrees=89.21708989130315),
        Om=Angle(degrees=282.9487539423989),
        w=Angle(degrees=130.662020526416),
        M=Angle(degrees=0.0),
        epoch=t,
        #mu_km_s=None,
        mu_km_s=mu,    #6.6743015e-20,
        #mu_au_d=2.9591220828559093E-04,
        center=None,
        target=None,
        center_name=None,
        target_name=None,
    )
    r, v = k._at(t)[:2]

    from skyfield.data.spice import inertial_frames
    #r1 = inertial_frames['ECLIPJ2000'].dot(r)
    r2 = inertial_frames['ECLIPJ2000'].T.dot(r)

    #print(v)

    planets = load('de421.bsp')
    sun = planets['sun']
    sunp = sun.at(t)

    from numpy import array
    horizons_au = array([
        -1.232674024434804E-01, 2.349174352473917E-01, 8.796973894528012E-01,
    ])
    print('have', r2)
    # print('have2', r2, sunp.position.au)
    # print('have2', r2 - sunp.position.au)
    # print('have2', r2 + sunp.position.au)
    print('want', horizons_au)
    print('delta:', r2 - horizons_au)
    print('sun:  ', sunp.position.au)

    # assert (r2.round(1) == [-0.1, 0.3, 0.9]).all()
    # assert (horizons_au.round(1) ==  [-0.1, 0.2, 0.9]).all()

    t = ts.tdb_jd(2454724.5)
    print(t.utc_jpl())
    r2008a, v = k._at(t)[:2]
    r2008b = inertial_frames['ECLIPJ2000'].T.dot(r2008a)

    HORI = 'X= 1.777310651+00  Y= 1.638390146+00  Z=-2.712743223+01'

    print('2008 have', r2008b)
    print('2008 want', HORI)
    sunp = sun.at(t)
    print('sun:  ', sunp.position.au)

    rr = k.position_at_epoch.au
    rr = inertial_frames['ECLIPJ2000'].T.dot(rr)

    print('Internal have', rr)
    print('Internal want', horizons_au)

    asdf

def OFF_test_against_horizons_for_hale_bopp():
    # See the following file in the Skyfield repository:
    # horizons/hale-bopp-vector

    ts = load.timescale(builtin=True)
    t = ts.tdb_jd(2450538.4378482755)
    print(t.utc_jpl())

    e = 0.9949607008417696
    qr = 0.9174143409263262  # Perihelion distance
    semimajor_axis_au = qr / (1.0 - e)

    a = semimajor_axis_au
    p_au = a * (1 - e ** 2)  # Wikipedia

    mu_m_s = 1.327124400189e20
    mu_km_s = mu_m_s * 1e-9
    mu = mu_km_s

    k = KeplerOrbit.from_mean_anomaly(
        p=Distance(au=p_au),  # see above
        e=e,
        i=Angle(degrees=89.21708989130315),
        Om=Angle(degrees=282.9487539423989),
        w=Angle(degrees=130.662020526416),
        M=Angle(degrees=0.0),
        epoch=t,
        #mu_km_s=None,
        mu_km_s=mu,    #6.6743015e-20,
        #mu_au_d=2.9591220828559093E-04,
        center=None,
        target=None,
        center_name=None,
        target_name=None,
    )
    r, v = k._at(t)[:2]

    from skyfield.data.spice import inertial_frames
    #r1 = inertial_frames['ECLIPJ2000'].dot(r)
    r2 = inertial_frames['ECLIPJ2000'].T.dot(r)

    #print(v)

    planets = load('de421.bsp')
    sun = planets['sun']
    sunp = sun.at(t)

    from numpy import array
    horizons_au = array([
        -1.232674024434804E-01, 2.349174352473917E-01, 8.796973894528012E-01,
    ])
    print('have', r2)
    # print('have2', r2, sunp.position.au)
    # print('have2', r2 - sunp.position.au)
    # print('have2', r2 + sunp.position.au)
    print('want', horizons_au)
    print('delta:', r2 - horizons_au)
    print('sun:  ', sunp.position.au)

    # assert (r2.round(1) == [-0.1, 0.3, 0.9]).all()
    # assert (horizons_au.round(1) ==  [-0.1, 0.2, 0.9]).all()

    t = ts.tdb_jd(2454724.5)
    print(t.utc_jpl())
    r2008a, v = k._at(t)[:2]
    r2008b = inertial_frames['ECLIPJ2000'].T.dot(r2008a)

    HORI = 'X= 1.777310651+00  Y= 1.638390146+00  Z=-2.712743223+01'

    print('2008 have', r2008b)
    print('2008 want', HORI)
    sunp = sun.at(t)
    print('sun:  ', sunp.position.au)

    rr = k.position_at_epoch.au
    rr = inertial_frames['ECLIPJ2000'].T.dot(rr)

    print('Internal have', rr)
    print('Internal want', horizons_au)

    asdf

# Test various round-trips through the kepler orbit object.

def _data_path(filename):
    return os.path.join(os.path.dirname(__file__), 'data', filename)

def check_orbit(p, e, i, Om, w, v,
                p_eps=None, e_eps=None, i_eps=None, Om_eps=None, w_eps=None, v_eps=None):
    pos0, vel0 = ele_to_vec(p, e, i, Om, w, v, mu)

    pos1, vel1 = propagate(pos0, vel0, 0, times, mu)
    ele = OsculatingElements(Distance(km=pos1), Velocity(km_per_s=vel1), dummy_time, mu)

    if p_eps: compare(p, ele.semi_latus_rectum.km, p_eps)
    if e_eps: compare(e, ele.eccentricity, e_eps)
    if i_eps: compare(i, ele.inclination.radians, i_eps, mod=True)
    if Om_eps: compare(Om, ele.longitude_of_ascending_node.radians, Om_eps, mod=True)
    if w_eps: compare(w, ele.argument_of_periapsis.radians, w_eps, mod=True)
    if v_eps: compare(v, ele.true_anomaly.radians, v_eps, mod=True)


times = linspace(-1e11, 1e11, 1001) # -3170 years to +3170 years, including 0
mu = 403503.2355022598
dummy_time = load.timescale().utc(2018)

def test_circular():
    check_orbit(p=300000, e=0, i=.5, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_circular_equatorial():
    check_orbit(p=300000, e=0, i=0, Om=0, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15)


def test_circular_polar():
    check_orbit(p=300000, e=0, i=pi/2, Om=1, w=0, v=1,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15)


def test_elliptical():
    check_orbit(p=300000, e=.3, i=1, Om=0, w=4, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)


def test_elliptical_equatorial():
    check_orbit(p=300000, e=.3, i=0, Om=0, w=1, v=5,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-7)


def test_elliptical_polar():
    check_orbit(p=300000, e=.2, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-2, e_eps=1e-8, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-8)


def test_parabolic():
    check_orbit(p=300000, e=1, i=1, Om=0, w=4, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-13, Om_eps=1e-13, w_eps=1e-13)


def test_parabolic_equatorial():
    check_orbit(p=300000, e=1, i=0, Om=0, w=1, v=2,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-13)


def test_parabolic_polar():
    check_orbit(p=300000, e=1, i=pi/2, Om=1, w=2, v=3,
                p_eps=1e-5, e_eps=1e-14, i_eps=1e-14, Om_eps=1e-13, w_eps=1e-13)


def test_hyperbolic():
    check_orbit(p=300000, e=1.3, i=1, Om=0, w=4, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)


def test_hyperbolic_equatorial():
    check_orbit(p=300000, e=1.3, i=0, Om=0, w=1, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-15, Om_eps=1e-15, w_eps=1e-6)


def test_hyperbolic_polar():
    check_orbit(p=300000, e=1.3, i=pi/2, Om=1, w=2, v=.5,
                p_eps=1e0, e_eps=1e-6, i_eps=1e-10, Om_eps=1e-10, w_eps=1e-6)
