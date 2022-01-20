# import nyx_space
# from nyx_space import dynamics
# SpacecraftDynamics = dynamics.spacecraft.SpacecraftDynamics
# SolarPressure = dynamics.solarpressure.SolarPressure
# Drag = dynamics.drag.Drag

from nyx_space import cosmic
# # from nyx_space import propagators

# Bodies = cosmic.Bodies
# # from nyx_space.dynamics.spacecraft import SpacecraftDynamics
# # from nyx_space.dynamics.solarpressure import SolarPressure
# # from nyx_space.dynamics.drag import Drag
# # from nyx_space.cosmic import Bodies
# from nyx_space import time
# from nyx_space import io
# # from nyx_space import Spacecraft

# Initialize the cosm which stores the ephemeris
cosm = cosmic.Cosm.de438()
# print(cosm.frame("ecliptic"))
print(cosm.frames_get_names())

# # Grab the frames we'll use
# eme2k = cosm.frame("EME2000")
# iau_earth = cosm.frame("IAU Earth")

# # Define the epoch
# epoch = time.Epoch.from_gregorian_utc(2014, 7, 22, 11, 29, 10, 811_000)

# # Define the initial orbit
# orbit = cosmic.Orbit.cartesian(
#     -137380.1984338506,
#     75679.87867537055,
#     21487.63875187856,
#     -0.2324532014235503,
#     -0.4462753967758019,
#     0.08561205662877103,
#     epoch,
#     eme2k,
# )

# # Define the spacecraft
# # sat = Spacecraft(orbit, 1000.0, 0.0, 1.0, 15.0, 1.7, 2.2)


# # Set up the harmonics first because we need to pass them to the overarching orbital dynamics
# # Load the harmonics from the JGM3 file (GMAT uses the JGM2 in this case).
# # It's gunzipped (hence `true` as the last parameter)
# stor = io.gravity.HarmonicsMem.from_cof("data/JGM3.cof.gz", 20, 20, True)
# # Set up the orbital dynamics: we need to specify the models one by one here
# # because the usual functions wrap the dynamics so that they can be used in a Monte Carlo
# # setup.
# # orbital_dyn = OrbitalDynamics::new(vec![
# #     # Note that we are only accounting for Sun, Moon and Jupiter, in addition to the integration frame's GM
# #     PointMasses::new(
# #         eme2k,
# #         &[Bodies::Sun, Bodies::Luna, Bodies::JupiterBarycenter],
# #         cosm.clone(),
# #     ),
# #     # Specify that these harmonics are valid only in the IAU Earth frame. We're using the
# #     Harmonics::from_stor(iau_earth, stor, cosm.clone()),
# # ]);
# orbital_dyn = dynamics.orbital.OrbitalDynamics.point_masses(
#     # Note that we are only accounting for Sun, Moon and Jupiter, in addition to the integration frame's GM
#     [Bodies.Sun, Bodies.Luna, Bodies.JupiterBarycenter],
#     cosm,
# )


# # Set up SRP and Drag second, because we need to pass them to the overarching spacecraft dynamics
# srp = SolarPressure.default(eme2k, cosm).force_model()
# drag = Drag.std_atm1976(cosm).force_model()
# # Set up the spacecraft dynamics
# sc_dyn = SpacecraftDynamics.from_models(orbital_dyn, [srp, drag])

# # prop = propagators.Propagator.default(sc_dyn)

# # (out, traj) = prop \
# #     .with(sat) \
# #     .for_duration_sec(0.5 * 3_600 * 24) # Half a day

# from nyx_space import md
# from nyx_space import io

# import logging

# def init_logging():
#     FORMAT = '%(levelname)s %(name)s %(asctime)-15s %(filename)s:%(lineno)d %(message)s'
#     logging.basicConfig(format=FORMAT)
#     logging.getLogger().setLevel(logging.INFO)

# with open('data/simple-scenario.toml', 'r') as f:
#     scen_data = f.read()

# scenario = io.ScenarioSerde.from_toml_str(scen_data)

# md.MDProcess.execute_all_in_scenario(scenario, cosm)
