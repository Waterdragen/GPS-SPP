from datetime import timedelta, datetime
from math import atan, atan2, asin, cos, pi, sqrt, sin, tan

import numpy as np

from consts import Consts, WGS84


def date_num(dt: datetime) -> int:
    return dt.toordinal() + 366


def time_num(dt: datetime) -> int:
    return int(timedelta(days=dt.isoweekday() % 7,
                         hours=dt.hour,
                         minutes=dt.minute,
                         seconds=dt.second).total_seconds())


def transmission_sec(distance: float) -> float:
    return distance / Consts.c


def mean_motion(sqrt_a: float, delta_n: float) -> float:
    # n = sqrt(mu / a^3)
    #   = sqrt(mu) / sqrt(a^3)
    #   = sqrt(mu) / sqrt(a) ^ 3
    # n = n0 + Δn
    return sqrt(Consts.mu) / sqrt_a ** 3 + delta_n


def time_diff_k(t: float, t0: float) -> float:
    tk = t - t0
    if tk > Consts.half_week:
        tk -= 2 * Consts.half_week
    elif tk < - Consts.half_week:
        tk += 2 * Consts.half_week
    return tk


def mean_anomaly(m0: float, n: float, tk: float) -> float:
    # Mk = m0 + tk * (n0 + Δn)
    return m0 + tk * n


def ecc_anomaly(m: float, ecc: float) -> float:
    # E = M + e * sin(E)
    # firstly, assume e * sin(M) is very small ~= 0
    # E = M
    # then we can work bottom-up by substituting E to get the new E
    iterations = 30
    e = m
    for _ in range(iterations):
        e = m + ecc * sin(e)

    return e


def true_anomaly(E: float, ecc: float) -> float:
    # θ = 2 * atan(sqrt(1 + e) / sqrt(1 - e) * tan(E / 2))
    return 2 * atan(sqrt((1 + ecc) / (1 - ecc)) * tan(E / 2))


def argument_of_latitude(theta: float, omega: float) -> float:
    return theta + omega


def orbit_radius(sqrt_a, ecc, E, crs, crc, phi) -> float:
    r = sqrt_a * sqrt_a * (1 - ecc * cos(E))
    delta_r = crs * sin(2 * phi) + crc * cos(2 * phi)
    return r + delta_r


def inclination(i0, d_i, tk, cis, cic, phi) -> float:
    i = i0 + d_i * tk
    delta_i = cis * sin(2 * phi) + cic * cos(2 * phi)
    return i + delta_i


def argument_of_latitude_corrected(cus, cuc, phi) -> float:
    delta_phi = cus * sin(2 * phi) + cuc * cos(2 * phi)
    return phi + delta_phi


def gps_position_orbital_plane(r, phi_c) -> tuple[float, float]:
    x0 = r * cos(phi_c)
    y0 = r * sin(phi_c)
    return x0, y0


def longitude_of_ascending_node(omega_0, d_omega, tk, t0) -> float:
    return omega_0 + (d_omega - Consts.omega_e) * tk - Consts.omega_e * t0


def gps_position_ecef(x0, y0, omega_k, i) -> tuple[float, float, float]:
    x = x0 * cos(omega_k) - y0 * cos(i) * sin(omega_k)
    y = x0 * sin(omega_k) + y0 * cos(i) * cos(omega_k)
    z = y0 * sin(i)
    return x, y, z

def satellite_clock_error(a0, a1, a2, tk) -> float:
    return Consts.c * (a0 + a1 * tk + a2 * tk * tk)


def correct_earth_rotation(pos: tuple[float, float, float], time_transmission) -> tuple[float, float, float]:
    # Rotated angle ωτ (rad) = rotation speed Ω_e (rad/s) * time (s)
    rotated_angle = Consts.omega_e * time_transmission
    x, y, z = pos
    new_x = x * cos(rotated_angle) + y * sin(rotated_angle)
    new_y = x * -sin(rotated_angle) + y * cos(rotated_angle)
    return new_x, new_y, z


def geometric_distance(base_pos: np.array, gps_pos_list: np.ndarray) -> np.ndarray:
    return np.sqrt((base_pos[0] - gps_pos_list[:, 0]) ** 2 +
                   (base_pos[1] - gps_pos_list[:, 1]) ** 2 +
                   (base_pos[2] - gps_pos_list[:, 2]) ** 2)


def design_matrix(base_pos: np.array, gps_pos_list: np.ndarray, rho: np.ndarray) -> np.ndarray:
    # Or the "unit vector", â = a / |a| where a is dx, dy, dz, 1, and |a| is rho
    rows, _ = gps_pos_list.shape
    b = np.empty((rows, 4), dtype=float)
    b[:, 0] = (base_pos[0] - gps_pos_list[:, 0]) / rho
    b[:, 1] = (base_pos[1] - gps_pos_list[:, 1]) / rho
    b[:, 2] = (base_pos[2] - gps_pos_list[:, 2]) / rho
    b[:, 3] = np.ones_like(rho)
    return b


def xyz_to_phi_lambda_h(pos: tuple[float, float, float]) -> tuple[float, float, float]:
    a, f, ecc_2 = WGS84.a, WGS84.f, WGS84.ecc_2
    x, y, z = pos
    p = sqrt(x * x + y * y)
    lam = atan2(y, x)

    r = sqrt(p * p + z * z)
    sin_phi = z / r
    phi = asin(sin_phi)
    h = r - a * (1 - sin_phi * sin_phi * f)
    tolerance = 1e-10

    while True:
        sin_phi = sin(phi)
        cos_phi = cos(phi)
        n = a / sqrt(1 - ecc_2 * sin_phi * sin_phi)  # prime vertical
        dp = p - (n + h) * cos_phi
        dz = z - (n * (1 - ecc_2) + h) * sin_phi
        h += sin_phi * dz + cos_phi * dp
        phi += (cos_phi * dz - sin_phi * dp) / (n + h)
        if (dp * dp + dz * dz) < tolerance:
            break

    return phi, lam, h  # phi and lambda are in radians


def gps_satellites_in_enu(base_pos_xyz: np.array,
                          base_pos_phi_lambda_h: tuple[float, float, float],
                          gps_pos_list: np.ndarray) -> np.ndarray:
    # The transformation from global geodetic to local geodetic coordinate system
    # can be calculated as below:
    #
    # [ -sin(λ)        cos(λ)        0      ]   [X_sat - X0]
    # [ -sin(φ)cos(λ) -sin(φ)sin(λ)  cos(φ) ] * [Y_sat - Y0]
    # [  cos(φ)cos(λ)  cos(φ)sin(λ)  sin(φ) ]   [Z_sat - Z0]
    #
    # where X0, Y0, Z0, φ, λ are the base position parameters

    x, y, z = base_pos_xyz[0], base_pos_xyz[1], base_pos_xyz[2]
    phi, lam, _ = base_pos_phi_lambda_h
    sin_phi, cos_phi, sin_lam, cos_lam = sin(phi), cos(phi), sin(lam), cos(lam)
    rows, _ = gps_pos_list.shape
    gps_sat_enu = np.empty((rows, 3), dtype=float)
    gps_sat_enu[:, 0] = ((gps_pos_list[:, 0] - x) * -sin_lam +
                         (gps_pos_list[:, 1] - y) * cos_lam)
    gps_sat_enu[:, 1] = ((gps_pos_list[:, 0] - x) * -sin_phi * cos_lam +
                         (gps_pos_list[:, 1] - y) * -sin_phi * sin_lam +
                         (gps_pos_list[:, 2] - z) * cos_phi)
    gps_sat_enu[:, 2] = ((gps_pos_list[:, 0] - x) * cos_phi * cos_lam +
                         (gps_pos_list[:, 1] - y) * cos_phi * sin_lam +
                         (gps_pos_list[:, 2] - z) * sin_phi)
    return gps_sat_enu


def elevation_angle_from_enu(gps_sat_enu: np.ndarray) -> np.ndarray:
    # the elevation angle (radians) can be calculated as follows:
    # d = sqrt(E ^ 2 + N ^ 2)
    # theta = atan(U / d)
    return np.arctan2(gps_sat_enu[:, 2],
                      np.sqrt(gps_sat_enu[:, 0] ** 2 + gps_sat_enu[:, 1] ** 2))


def klobuchar_adjustment(base_pos_xyz: np.array,
                         base_pos_phi_lambda_h: tuple[float, float, float],
                         gps_pos_list: np.ndarray) -> np.ndarray:
    # This is a coarse klobuchar adjustments, which does not account for the time of the day
    # And we use the night model
    # therefore, klobuchar coefficients are not used here
    gps_sat_enu: np.ndarray = gps_satellites_in_enu(base_pos_xyz,
                                                    base_pos_phi_lambda_h,
                                                    gps_pos_list)
    elev: np.ndarray = elevation_angle_from_enu(gps_sat_enu)  # in radians

    # slant factor
    # (radians / pi) is the number of half-turns
    f = 1.0 + 16.0 * (0.53 - elev / pi) ** 3

    # 5 light-nanoseconds slant distance (m)
    return Consts.c * f * 5e-9
