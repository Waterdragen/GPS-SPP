def parse_sci_float(s: str) -> float:
    return float(s[:15]) * 10 ** int(s[16:])

def parse_sat_prn(s: str) -> list[int]:
    # collect the numbers in the string into a list
    # for example: "10G10G22G17G15G18G 6G26G28G23G 3" becomes
    # [10, 22, 17, 15, 18, 6, 26, 28, 23, 3]
    # and there are 10 satellites
    num_satellites = int(s[30:32])
    return [int(s[33 + i * 3: 35 + i * 3]) for i in range(num_satellites)]

def parse_time(s: str, four_digit_year=False) -> tuple[int, int, int, int, int, float]:
    # helper function to parse a space delimited time
    tme = s.split()
    year = int(tme[0])
    if not four_digit_year:
        year += 1900 if year >= 80 else 2000
    return (year, int(tme[1]), int(tme[2]),
            int(tme[3]), int(tme[4]), float(tme[5]))