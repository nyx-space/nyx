from nyx_space import cosmic
from nyx_space import time


cosm = cosmic.Cosm.de438()
print(cosm.frames_get_names())

eme2k = cosm.frame("EME2000")

e = time.Epoch.system_now()
print(e)

orbit = cosmic.Orbit.from_keplerian_altitude(
    400,
    ecc=1e-4,
    inc_deg=30.5,
    raan_deg=35.0,
    aop_deg=65.0,
    ta_deg=590,
    epoch=e,
    frame=eme2k,
)
print(orbit)
print(repr(orbit))

# Build a spacecraft
sc = cosmic.Spacecraft(
    orbit,
    dry_mass_kg=500.0,
    fuel_mass_kg=15.0,
    srp_area_m2=2.0,
    drag_area_m2=2.0,
    cr=1.8,
    cd=2.1,
)

print(sc)

print(sc.orbit)

# NOTE: This does not return a pointer to the object, but a new object!
orbit = sc.orbit

orbit.sma_km = sc.orbit.sma() + 100.0
sc.orbit = orbit

print(sc.orbit)