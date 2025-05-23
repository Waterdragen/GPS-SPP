import csv
from datetime import timedelta, datetime
from math import degrees, inf
import numpy as np
from os import path
from time import perf_counter

import math_fn
from core import ObsData, NavData
from helper_fn import parse_sat_prn, parse_time
from util import group_by_duplicates, list_iter, find_first, pretty_format_degree

OBS_FILE = "site0900.01o"
NAV_FILE = "site0900.01n"
OBS_CSV_FILE = "site0900-01o.csv"
NAV_CSV_FILE = "site0900-01n.csv"

if not path.exists(OBS_FILE):
    print(f"Please make sure '{OBS_FILE}' is placed in the same directory as main.py, aborting...")
    exit(1)
if not path.exists(NAV_FILE):
    print(f"Please make sure '{NAV_FILE}' is placed in the same directory as main.py, aborting...")
    exit(1)


def jump_to_suffix_line(file_obj, suffix: str) -> str:
    suffix += "\n"
    line: str
    for line in file_obj:
        if line.endswith(suffix):
            return line.rstrip()
    raise RuntimeError("Cannot find suffix " + suffix)


def jump_to_containing_line(file_obj, s: str) -> str:
    for line in file_obj:
        if s in line:
            return line.rstrip()
    raise RuntimeError("Cannot find text " + s)


def read_obs() -> int:
    with open(OBS_FILE) as f, open(OBS_CSV_FILE, 'w') as table:
        table.write(ObsData.header_row())

        text_types_of_satellites = jump_to_suffix_line(f, "# / TYPES OF OBSERV")
        type_obs_count = int(text_types_of_satellites[1:6])
        type_obs_names = [t for t in text_types_of_satellites[6:60].split()]
        c1_index = type_obs_names.index("C1")

        text_time_of_first_obs = jump_to_suffix_line(f, "TIME OF FIRST OBS")
        *first_obs_dt, second = parse_time(text_time_of_first_obs, four_digit_year=True)
        first_obs_time = datetime(*first_obs_dt, int(second))

        print(f"First observation time: {first_obs_time}")

        _ = jump_to_suffix_line(f, "END OF HEADER")

        current_line = f.readline()
        number_of_obs = 0  # for allocating exact size later

        while current_line != "":
            # Check the epoch flag is 0 (OK)
            epoch_flag = current_line[27:30]
            if epoch_flag != " 0 ":
                current_line = f.readline()
                continue

            year, month, day, hour, minute, second = parse_time(current_line[:26])
            dt = datetime(year, month, day, hour, minute, int(second))
            date_num = math_fn.date_num(dt)
            time_num = math_fn.time_num(dt)

            sat_prn = parse_sat_prn(current_line)

            for prn in sat_prn:
                lines = ""
                new_line = f.readline()
                # Add an arbitrary character to make it 80, multiple of 16
                lines += new_line + "\n"
                while len(lines) == 80:
                    new_line = f.readline()
                    lines += new_line + "\n"

                obs = [float(lines[i * 16: i * 16 + 14]) for i in range(type_obs_count)]

                obs_data = ObsData(prn, year, month, day, hour, minute, second,
                                   int(epoch_flag),
                                   obs[c1_index],
                                   date_num, time_num)
                table.write(obs_data.values_row())
                number_of_obs += 1

            current_line = f.readline()

    print(f"Successfully wrote observation data to {OBS_CSV_FILE}")
    return number_of_obs


def read_nav() -> list[NavData]:
    nav_data_list = []

    with open(NAV_FILE) as f, open(NAV_CSV_FILE, 'w') as table:
        table.write(NavData.header_row())

        _ = jump_to_containing_line(f, "END OF HEADER")
        try:
            while True:
                nav_data = NavData.parse(f)
                nav_data_list.append(nav_data)
                table.write(nav_data.values_row())
        except StopIteration:
            pass

    print(f"Successfully wrote navigation data to {NAV_CSV_FILE}")
    return nav_data_list


def process_gps_position(nav_data_list: list[NavData], obs_file, number_of_obs: int):
    # The navigation file is always sorted by time, record the ranges of timeframes with slices
    grouped_nav_data_list, nav_slice_list = group_by_duplicates(iter(nav_data_list),
                                                                lambda nav1, nav2: nav1.time() == nav2.time())
    nav_time_gen = (nav_data.time() for nav_data in grouped_nav_data_list)
    # Maps time to slice of original list
    time_map = [(t, s) for t, s in zip(nav_time_gen, nav_slice_list)]

    gps_pos_list = np.empty((number_of_obs, 3), dtype=float)
    nav_list = [None] * number_of_obs
    obs_list = [None] * number_of_obs
    c1_list = np.empty(number_of_obs, dtype=float)
    sat_clock_err_list = np.empty(number_of_obs, dtype=float)

    for index, obs_row in enumerate(obs_file):
        obs_data = ObsData.from_csv_row(obs_row)
        obs_time = obs_data.time()
        # Sort time-slice map by closest time first
        time_map.sort(key=lambda item: abs(item[0] - obs_time))
        match_nav_data = None
        for _, s in time_map:
            # Linear search within the same time frame
            match_nav_data = find_first(list_iter(nav_data_list, s), lambda nav_data: nav_data.prn == obs_data.prn)
            if match_nav_data is not None:
                break

        # No satellites with matching PRN found
        if match_nav_data is None:
            continue
        # Found a matching satellite with matching PRN, but expired
        if abs(match_nav_data.time() - obs_time) >= timedelta(hours=4):
            continue

        x, y, z, sat_clock_err = calculate_gps_position(match_nav_data, obs_data)
        gps_pos_list[index, 0] = x
        gps_pos_list[index, 1] = y
        gps_pos_list[index, 2] = z
        nav_list[index] = match_nav_data
        obs_list[index] = obs_data
        c1_list[index] = obs_data.c1
        sat_clock_err_list[index] = sat_clock_err

    ecef_pos, wgs84_pos = least_square_adjustment(gps_pos_list, nav_list, obs_list, c1_list, sat_clock_err_list)
    pretty_print_results(ecef_pos, wgs84_pos)


def calculate_gps_position(nav_data: NavData,
                           obs_data: ObsData,
                           rx_clock_error=0.0) -> tuple[float, float, float, float]:
    # Time from ephemeris reference epoch (tk)
    transmission_sec = math_fn.transmission_sec(obs_data.c1)
    t = obs_data.time_num - transmission_sec - rx_clock_error
    tk = math_fn.time_diff_k(t, nav_data.t0)
    satellite_clock_error = math_fn.satellite_clock_error(nav_data.sv_clock_bias,
                                                          nav_data.sv_clock_drift,
                                                          nav_data.sv_clock_drift_rate, tk)

    # Mean motion (n)
    n = math_fn.mean_motion(nav_data.sqrt_a, nav_data.delta_n)

    # Mean anomaly (M)
    M = math_fn.mean_anomaly(nav_data.m0, n, tk)

    # Eccentric anomaly (E)
    E = math_fn.ecc_anomaly(M, nav_data.ecc)

    # True anomaly (θ)
    theta = math_fn.true_anomaly(E, nav_data.ecc)

    # Argument of latitude (φ)
    phi = math_fn.argument_of_latitude(theta, nav_data.omega)

    # Orbit radius of GPS satellite position (r)
    r = math_fn.orbit_radius(nav_data.sqrt_a, nav_data.ecc, E,
                             nav_data.crs, nav_data.crc, phi)

    # Corrected GPS orbit inclination (i)
    i = math_fn.inclination(nav_data.i0, nav_data.d_i, tk,
                            nav_data.cis, nav_data.cic, phi)

    # Corrected argument of latitude
    phi_c = math_fn.argument_of_latitude_corrected(nav_data.cus, nav_data.cuc, phi)

    # GPS Positions in the orbital plane (φ)
    x0, y0 = math_fn.gps_position_orbital_plane(r, phi_c)

    # Longitude of ascending node (Ω_k)
    omega_k = math_fn.longitude_of_ascending_node(nav_data.omega_0, nav_data.d_omega, tk, nav_data.t0)

    # GPS Positions in ECEF
    pos = math_fn.gps_position_ecef(x0, y0, omega_k, i)

    # Correct GPS Satellite positions with earth's rotation
    x, y, z = math_fn.correct_earth_rotation(pos, transmission_sec)

    return x, y, z, satellite_clock_error


def least_square_adjustment(gps_pos_list: np.ndarray,
                            nav_list,
                            obs_list,
                            c1_list: np.ndarray,
                            sat_clock_err_list: np.ndarray) -> tuple[tuple[float, float, float],
                                                                     tuple[float, float, float]]:
    # We will progressively use more satellites as the accuracy increases
    num_sats = 4  # Minimum satellites required
    total_sats, _ = gps_pos_list.shape
    total_reached = 0

    # Initial guess
    approx_pos = np.array([1.0, 1.0, 1.0], dtype=float)

    rx_clock_dist_prev = 0.0
    rx_clock_dist_adjusted = inf

    iono_adjustments = None  # be aware this must be initialized before use
    adjustments = np.array([inf, inf, inf, 0.0], dtype=float)

    # dx, dy, dz, dt all converge to < 0.1mm, and have used all the least square parameters at least once
    while abs(rx_clock_dist_adjusted) > 1e-4 or total_reached < 2 or np.sum(np.abs(adjustments[:3])) > 1e-4:
        # Recalculate with new clock error
        rx_clock_error = math_fn.transmission_sec(adjustments[3])
        for i in range(num_sats):
            x, y, z, sat_clock_err = calculate_gps_position(nav_list[i], obs_list[i], rx_clock_error)
            gps_pos_list[i, 0] = x
            gps_pos_list[i, 1] = y
            gps_pos_list[i, 2] = z
            sat_clock_err_list[i] = sat_clock_err

        # Actual geometric distance to the satellites
        rho = math_fn.geometric_distance(approx_pos, gps_pos_list[:num_sats])

        # Design matrix
        b = math_fn.design_matrix(approx_pos, gps_pos_list[:num_sats], rho)

        # Error function
        f = c1_list[:num_sats] - rho + sat_clock_err_list[:num_sats]
        if total_reached:
            f += iono_adjustments

        # Least squares 4 unknowns: dx, dy, dz, dt
        # (B^T * B) ^ -1 * B^T * f
        adjustments = np.linalg.inv(b.T @ b) @ b.T @ f

        approx_pos += adjustments[:3]
        rx_clock_dist_adjusted = adjustments[3] - rx_clock_dist_prev
        rx_clock_dist_prev = adjustments[3]

        num_sats *= 4
        # Limit the number of satellites, or jump to total satellites if the error is below 1 meter
        if num_sats > total_sats or abs(rx_clock_dist_adjusted) < 1:
            num_sats = total_sats
            base_pos = approx_pos[0], approx_pos[1], approx_pos[2]
            base_pos_wgs = math_fn.xyz_to_phi_lambda_h(base_pos)
            iono_adjustments = math_fn.klobuchar_adjustment(base_pos,
                                                            base_pos_wgs,
                                                            gps_pos_list)
            total_reached += 1

    x, y, z = approx_pos[0], approx_pos[1], approx_pos[2]
    phi, lam, h = math_fn.xyz_to_phi_lambda_h((x, y, z))

    return (x, y, z), (phi, lam, h)


def pretty_print_results(ecef_pos: tuple[float, float, float],
                         wgs84_pos: tuple[float, float, float]):
    print("Successfully calculated the GPS receiver position!")
    print()
    x, y, z = ecef_pos
    print(f"ECEF Coordinates\n"
          f"({x:.3f}, {y:.3f}, {z:.3f})")
    print()
    phi, lam, h = wgs84_pos
    phi = degrees(phi)
    lam = degrees(lam)
    print(f"WGS84 Datum\n"
          f"Latitude: {phi:.9f}\n"
          f"Longitude: {lam:.9f}\n"
          f"Height: {h:.3f}\n"
          f"\n"
          f"Paste these coordinates in Google map:\n"
          f"({pretty_format_degree(phi)}, {pretty_format_degree(lam)})")


def main():
    t_start = perf_counter()

    nav_data_list = read_nav()
    number_of_obs = read_obs()
    with open(OBS_CSV_FILE, newline='') as obs_file:
        obs_file = csv.reader(obs_file)
        next(obs_file)  # exclude header from csv
        process_gps_position(nav_data_list, obs_file, number_of_obs)

    t_end = perf_counter()
    print(f"Time elapsed: {t_end - t_start:.3f}s")


if __name__ == '__main__':
    main()
