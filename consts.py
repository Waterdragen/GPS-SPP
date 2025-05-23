class Consts:
    c = 299792458        # Speed of light (m/s)
    g = 9.80665          # Acceleration of earth's gravity (m/s^2)
    G = 6.67259e-11      # Universal gravitational constant (m^3/kg s^2)
    M = 5.972e24         # Mass of the Earth (kg)
    mu = 3.986004418e14  # Standard gravitational parameter (μ = G * M) (m^3/s^2)
    omega_e = 7.2921151467e-5  # Earth rotation rate (rad/s)
    
    half_week = 3.5 * 60 * 60 * 24

class WGS84:
    a = 6378137
    f = 1 / 298.257223563
    # e = sqrt(2f - f^2)
    ecc_2 = (2 - f) * f
