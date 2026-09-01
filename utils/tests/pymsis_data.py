from datetime import datetime, timezone
import json
from pymsis import msis


def generate_vectors():
    # 7-element Ap array: [daily_ap, current_3hr, -3hr, -6hr, -9hr, -12to33hr_avg, -36to57hr_avg]
    storm_ap = [15.0, 130.0, 150.0, 50.0, 20.0, 15.0, 15.0]
    flat_ap = [15.0, 15.0, 15.0, 15.0, 15.0, 15.0, 15.0]

    alts_km = [100.0, 120.0, 483.0]
    lats_deg = [0.0, 15.0, 45.0, 80.0]
    f107 = 150.0
    f107a = 150.0
    lon = 0.0

    test_cases = []

    for hour in [6, 12, 18]:
        date = datetime(2025, 1, 1, hour, 0, tzinfo=timezone.utc)
        for alt in alts_km:
            for lat in lats_deg:
                for is_storm in [True, False]:
                    ap_array = storm_ap if is_storm else flat_ap
                    geom_act = -1 if is_storm else 1

                    # pymsis run expects arrays.
                    # Output densities are in kg/m^3 (total mass) and m^-3 (species).
                    output = msis.run(
                        [date],
                        [lon],
                        [lat],
                        [alt],
                        [f107],
                        [f107a],
                        [ap_array],
                        geomagnetic_activity=geom_act,
                        version=0,
                    )
                    # Extract densities. pymsis output shape: (ndates, nlons, nlats, nalts, 11)
                    # Indices: 0: Total mass, 1: N2, 2: O2, 3: O, 4: He, 5: H, 6: Ar, 7: N, 8: Anomalous O, 9: T_exo, 10: T_alt
                    result = output[0, :]

                    test_cases.append(
                        {
                            "altitude_km": alt,
                            "latitude_deg": lat,
                            "longitude_deg": lon,
                            "day_of_year": date.timetuple().tm_yday,
                            "ut_seconds": date.hour * 3600.0 + date.minute * 60.0,
                            "f107_daily": f107,
                            "f107_avg": f107a,
                            "ap_array": ap_array,
                            "is_storm": is_storm,
                            "expected_total_density_kg_m3": float(result[0]),
                            "expected_temperature_k": float(result[10]),
                        }
                    )

    with open("data/03_tests/nrlmsise00_validation.json", "w") as f:
        json.dump(test_cases, f, indent=4)


if __name__ == "__main__":
    generate_vectors()
