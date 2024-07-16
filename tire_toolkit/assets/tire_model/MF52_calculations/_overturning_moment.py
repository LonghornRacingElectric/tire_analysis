import numpy as np

def get_M_x(overturning_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma, Fy) -> float:

    # Overturning moment coeffs
    QSX1 = overturning_coeffs["QSX1"]
    QSX2 = overturning_coeffs["QSX2"]
    QSX3 = overturning_coeffs["QSX3"]

    # Scaling coeffs
    LFZO = scaling_coeffs["LFZO"]
    LCX = scaling_coeffs["LCX"]
    LMUX = scaling_coeffs["LMUX"]
    LEX = scaling_coeffs["LEX"]
    LKX = scaling_coeffs["LKX"]
    LHX = scaling_coeffs["LHX"]
    LVX = scaling_coeffs["LVX"]
    LXAL = scaling_coeffs["LXAL"]
    LGAX = scaling_coeffs["LGAX"]
    LCY = scaling_coeffs["LCY"]
    LMUY = scaling_coeffs["LMUY"]
    LEY = scaling_coeffs["LEY"]
    LKY = scaling_coeffs["LKY"]
    LHY = scaling_coeffs["LHY"]
    LVY = scaling_coeffs["LVY"]
    LGAY = scaling_coeffs["LGAY"]
    LKYG = scaling_coeffs["LKYG"]
    LTR = scaling_coeffs["LTR"]
    LRES = scaling_coeffs["LRES"]
    LCZ = scaling_coeffs["LCZ"]
    LGAZ = scaling_coeffs["LGAZ"]
    LYKA = scaling_coeffs["LYKA"]
    LVYKA = scaling_coeffs["LVYKA"]
    LS = scaling_coeffs["LS"]
    LSGKP = scaling_coeffs["LSGKP"]
    LSGAL = scaling_coeffs["LSGAL"]
    LGYR = scaling_coeffs["LGYR"]
    LMX = scaling_coeffs["LMX"]
    LVMX = scaling_coeffs["LVMX"]
    LMY = scaling_coeffs["LMY"]
    LIP = scaling_coeffs["LIP"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]

    # Dimension coeffs
    R0 = dimensions["UNLOADED_RADIUS"]

    M_x = R0 * FZ * (QSX1 * LVMX + (-1 * QSX2 * gamma + QSX3 * Fy / FNOMIN) * LMX)

    return M_x