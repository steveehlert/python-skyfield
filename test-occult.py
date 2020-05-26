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
for this_ut in tqdm.tqdm(utrange):
    this_bool = this_sat.is_target_occulted(this_ut, target_pos, earth)
    bool_list.append(this_bool)

bool_mask = np.array(bool_list)

n_occulted = np.shape(np.where(bool_mask))[1]

frac_occulted = float(n_occulted)/float(n_steps)

print("A total of %.2f percent of the time steps were occulted" %(frac_occulted * 100.))

