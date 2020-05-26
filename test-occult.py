from skyfield.api import EarthSatellite, load, position_of_radec
from skyfield.nutationlib import iau2000b
import numpy as np
import tqdm






de421 = load('de421.bsp')
earth = de421['earth']
ts = load.timescale()

this_sat = EarthSatellite('1 00000U 21  0A   20146.43738583  .00000000  00000-0  00000-0 0  0000',
                          '2 00000   0.2000 000.0000 0000000   0.0000 000.0000 14.88336848949720', 'Test_600kmAlt', ts)

#100K times between 2022 and 2024

n_steps = 10000
jdrange = np.linspace(2459690.5,2460055.75,n_steps)


utrange = ts.ut1_jd(jdrange)
utrange._nutation_angles = iau2000b(utrange.tt)

#Crab Nebula
target_ra = 83.63308333
target_dec = 22.0145

#My RA is in degrees, hence the division by 15. 
target_pos = position_of_radec(target_ra / 15., target_dec)


#Currently not optimized for speed
bool_list = []
t_list = []
for this_ut in tqdm.tqdm(utrange):
    this_bool = this_sat.is_target_occulted(this_ut, target_pos, earth)
    dist_intersect = this_sat.intersect_line_with_sphere(this_ut, target_pos, earth)
    t_list.append(dist_intersect)
    bool_list.append(this_bool)

bool_mask = np.array(bool_list)
t_arr = np.array(t_list)

t_bool1 = np.isfinite(t_arr)
t_bool2 = t_arr > 0.

t_bool = np.logical_and(t_bool1,t_bool2)

print(np.array_equal(t_bool, bool_mask), np.array_equiv(t_bool,bool_mask))

t_finite = t_arr[np.isfinite(t_arr)]

t_intersect = t_finite[np.where(t_finite > 0.)]






n_occulted = np.shape(np.where(bool_mask))[1]
n_occulted2 = np.shape(t_intersect)[0]

frac_occulted = float(n_occulted)/float(n_steps)
frac_occulted2 = float(n_occulted2)/float(n_steps)
print("Steve's Method")
print("A total of %.2f percent of the time steps were occulted" %(frac_occulted * 100.))

print("Brandon's Method")
print("A total of %.2f percent of the time steps were occulted" %(frac_occulted2 * 100.))



