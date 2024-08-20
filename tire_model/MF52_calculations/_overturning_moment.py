import numpy as np

def get_M_x(overturning_coeffs, scaling_coeffs, vertical_coeffs, dimensions, FZ, gamma, Fy) -> float:

    # Overturning moment coeffs
    QSX1 = overturning_coeffs["QSX1"]
    QSX2 = overturning_coeffs["QSX2"]
    QSX3 = overturning_coeffs["QSX3"]

    # Scaling coeffs
    LMX = scaling_coeffs["LMX"]
    LVMX = scaling_coeffs["LVMX"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]

    # Dimension coeffs
    R0 = dimensions["UNLOADED_RADIUS"]

    M_x = R0 * FZ * (QSX1 * LVMX + (-1 * QSX2 * gamma + QSX3 * Fy / FNOMIN) * LMX)

    return M_x