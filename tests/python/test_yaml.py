from nyx_space.orbit_determination import GroundStation, GaussMarkov
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
    # Using loads
    unique = GroundStation(**fourth)
    assert unique.name == "Maspalomas, ES"

    # Or the constructor
    unique = GroundStation(**fourth)
    assert unique.name == "Maspalomas, ES"

    # Check pickle
    pkld = pickle.dumps(unique)
    unpkld = pickle.loads(pkld)
    assert unpkld.name == "Maspalomas, ES"


def test_gauss_markov():
    from_list = GaussMarkov.loads(list(gm_data.values()))
    assert len(from_list) == 2

    unique_range = GaussMarkov.loads(gm_data["range_noise_km"])[0]
    assert unique_range.tau == Unit.Day * 1.0

    # NOTE: We cannot pickle GaussMarkov because it includes hifitime Duration which is of type `builtin.Duration` and I can't change that.


if __name__ == "__main__":
    test_ground_station()
    test_gauss_markov()
