from nyx_space.orbit_determination import GroundStation, GaussMarkov, FltResid
from nyx_space.time import Unit

import pickle

gs_data = {
    "one": {
        "name": "Santiago, CL",
        "frame": "IAU Earth",
        "elevation_mask_deg": 5.0,
        "range_noise_km": {
            "tau": "24 h",
            "bias_sigma": 0.005,
            "steady_state_sigma": 0.0001,
        },
        "doppler_noise_km_s": {
            "tau": "24h",
            "bias_sigma": 5e-05,
            "steady_state_sigma": 1.5e-06,
        },
        "light_time_correction": False,
        "latitude_deg": -33.447487,
        "longitude_deg": -70.673676,
        "height_km": 0.5,
    },
    "two": {
        "name": "South Point, Hawaii, US",
        "frame": "IAU Earth",
        "latitude_deg": 18.911057,
        "longitude_deg": -155.681022,
        "height_km": 0.5,
        "elevation_mask_deg": 5.0,
        "range_noise_km": {
            "tau": "24 h",
            "bias_sigma": 0.005,
            "steady_state_sigma": 0.0001,
        },
        "doppler_noise_km_s": {
            "tau": "24 h",
            "bias_sigma": 5e-05,
            "steady_state_sigma": 1.5e-06,
        },
        "light_time_correction": False,
    },
    "three": {
        "name": "Dongara, AU",
        "frame": "IAU Earth",
        "latitude_deg": -29.251281,
        "longitude_deg": 114.934621,
        "height_km": 0.5,
        "elevation_mask_deg": 5.0,
        "range_noise_km": {
            "tau": "24 h",
            "bias_sigma": 0.005,
            "steady_state_sigma": 0.0001,
        },
        "doppler_noise_km_s": {
            "tau": "24 h",
            "bias_sigma": 5e-05,
            "steady_state_sigma": 1.5e-06,
        },
        "light_time_correction": False,
    },
    "four": {
        "name": "Maspalomas, ES",
        "frame": "IAU Earth",
        "latitude_deg": 27.760562,
        "longitude_deg": -15.586017,
        "height_km": 0.5,
        "elevation_mask_deg": 5.0,
        "range_noise_km": {
            "tau": "24 h",
            "bias_sigma": 0.005,
            "steady_state_sigma": 0.0001,
        },
        "doppler_noise_km_s": {
            "tau": "24 h",
            "bias_sigma": 5e-05,
            "steady_state_sigma": 1.5e-06,
        },
        "light_time_correction": False,
    },
}

gm_data = {
    "range_noise_km": {
        "tau": "24 h",
        "bias_sigma": 0.005,
        "steady_state_sigma": 0.0001,
    },
    "doppler_noise_km_s": {
        "tau": "24h",
        "bias_sigma": 5e-05,
        "steady_state_sigma": 1.5e-06,
    },
}


def test_ground_station():
    # Test that we can load a ground station from a list
    from_list = GroundStation.loads(list(gs_data.values()))
    assert len(from_list) == 4

    # Test that we can load a ground station from a dict
    from_dict = GroundStation.loads(gs_data)
    assert len(from_dict) == 4

    # Test that we can load a ground station from a single dict definition
    fourth = gs_data["four"]
    fourth.pop("range_noise_km")
    fourth.pop("doppler_noise_km_s")
    # Using the constructor
    unique = GroundStation(**fourth)
    assert unique.name == "Maspalomas, ES"

    # Check pickle
    pkld = pickle.dumps(unique)
    unpkld = pickle.loads(pkld)
    assert unpkld == unique


def test_gauss_markov():
    from_dict = GaussMarkov.loads(gm_data)
    assert len(from_dict) == 2

    from_list = GaussMarkov.loads(list(gm_data.values()))
    assert len(from_list) == 2

    unique_range = GaussMarkov.loads(gm_data["range_noise_km"])[0]
    assert unique_range.tau == Unit.Day * 1.0

    # NOTE: We can only pickle via the `dumps` and `loads` functions.
    # This looses the stateful info.

    pkld = pickle.dumps(unique_range)
    unpkl = pickle.loads(pkld)
    assert unpkl.tau == Unit.Day * 1.0


def test_flt_resid():
    # Test pickle
    flt = FltResid(min_accepted=10, num_sigmas=3.0)
    unplkd = pickle.loads(pickle.dumps(flt))
    assert unplkd == flt


if __name__ == "__main__":
    test_ground_station()
    test_gauss_markov()
    test_flt_resid()
