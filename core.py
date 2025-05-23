from dataclasses import dataclass
from datetime import datetime
from helper_fn import parse_sci_float, parse_time

FIELD1 = slice(3, 22)
FIELD2 = slice(22, 41)
FIELD3 = slice(41, 60)
FIELD4 = slice(60, 79)
    
class KeyValueCSV:
    @classmethod
    def header_row(cls) -> str:
        return ",".join(cls.__annotations__.keys()) + '\n'

    def values_row(self):
        return ",".join(str(i) for i in vars(self).values()) + '\n'


@dataclass
class ObsData(KeyValueCSV):
    prn: int
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float
    epoch_flag: int
    c1: float
    date_num: int
    time_num: int

    @classmethod
    def from_csv_row(cls, csv_row: list[str]):
        return cls(
            int(csv_row[0]),
            int(csv_row[1]),
            int(csv_row[2]),
            int(csv_row[3]),
            int(csv_row[4]),
            int(csv_row[5]),
            float(csv_row[6]),
            int(csv_row[7]),
            float(csv_row[8]),
            int(csv_row[9]),
            int(csv_row[10]),
        )

    def time(self) -> datetime:
        int_second = int(self.second)
        micro_second = round((self.second - int_second) * 1_000_000)
        return datetime(self.year, self.month, self.day, self.hour, self.minute, int_second, micro_second)


@dataclass
class NavData(KeyValueCSV):
    prn: int
    year: int
    month: int
    day: int
    hour: int
    minute: int
    second: float
    sv_clock_bias: float
    sv_clock_drift: float
    sv_clock_drift_rate: float

    iode: float     # Issue of Data, Ephemeris
    crs: float      # correction terms to orbital radius
    delta_n: float  # mean motion difference
    m0: float       # mean anomaly

    cuc: float      # latitude argument correction terms
    ecc: float      # eccentricity
    cus: float      # latitude argument correction terms
    sqrt_a: float   # sqrt of semi-major axis

    t0: float       # or toe, Time of Ephemeris
    cic: float      # inclination correction terms
    omega_0: float  # right ascension
    cis: float      # inclination correction terms

    i0: float       # inclination
    crc: float      # orbital radius correction terms
    omega: float    # perigee argument
    d_omega: float  # or "dΩ/dt", "omega dot", rate of right ascension

    d_i: float      # or "di/dt", rate of inclination
    code_l2: float  # codes for L2 channels
    week_no: float
    l2_p_data_flag: float

    sv_accuracy: float
    sv_health: float
    tgd: float
    iodc: float     # Issue of Data, Clock

    transmission_time_of_msg: float
    fit_interval: float
    spare_1: float
    spare_2: float

    @classmethod
    def parse(cls, f):
        line = next(f)
        prn = int(line[:2])
        year, month, day, hour, minute, second = parse_time(line[2:22])

        sv_clock_bias = parse_sci_float(line[FIELD2])
        sv_clock_drift = parse_sci_float(line[FIELD3])
        sv_clock_drift_rate = parse_sci_float(line[FIELD4])

        misc = []
        for _ in range(7):
            line = next(f)
            misc.append([
                parse_sci_float(line[FIELD1]),
                parse_sci_float(line[FIELD2]),
                parse_sci_float(line[FIELD3]),
                parse_sci_float(line[FIELD4]),
            ])
        return cls(
            prn, year, month, day, hour, minute, second, sv_clock_bias, sv_clock_drift, sv_clock_drift_rate,
            misc[0][0], misc[0][1], misc[0][2], misc[0][3],  # broadcast orbit 1
            misc[1][0], misc[1][1], misc[1][2], misc[1][3],  # broadcast orbit 2
            misc[2][0], misc[2][1], misc[2][2], misc[2][3],  # broadcast orbit 3
            misc[3][0], misc[3][1], misc[3][2], misc[3][3],  # broadcast orbit 4
            misc[4][0], misc[4][1], misc[4][2], misc[4][3],  # broadcast orbit 5
            misc[5][0], misc[5][1], misc[5][2], misc[5][3],  # broadcast orbit 6
            misc[6][0], misc[6][1], misc[6][2], misc[6][3],  # broadcast orbit 7
        )

    def time(self) -> datetime:
        int_second = int(self.second)
        micro_second = round((self.second - int_second) * 1_000_000)
        return datetime(self.year, self.month, self.day, self.hour, self.minute, int_second, micro_second)
