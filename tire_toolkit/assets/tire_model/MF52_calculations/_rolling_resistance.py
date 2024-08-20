import numpy as np

def get_M_y(rolling_coeffs, long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, FZ):

    # Rolling resistance coeffs
    QSY1 = rolling_coeffs["QSY1"]
    QSY2 = rolling_coeffs["QSY2"]
    QSY3 = rolling_coeffs["QSY3"]
    QSY4 = rolling_coeffs["QSY4"]

    # Pure long coeffs
    PKX1 = long_coeffs["PKX1"]
    PKX2 = long_coeffs["PKX2"]
    PKX3 = long_coeffs["PKX3"]
    PHX1 = long_coeffs["PHX1"]
    PHX2 = long_coeffs["PHX2"]
    PVX1 = long_coeffs["PVX1"]
    PVX2 = long_coeffs["PVX2"]

    # Scaling coeffs
    LFZO = scaling_coeffs["LFZO"]
    LMUX = scaling_coeffs["LMUX"]
    LKX = scaling_coeffs["LKX"]
    LHX = scaling_coeffs["LHX"]
    LVX = scaling_coeffs["LVX"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]

    # Dimension coeffs
    R0 = dimensions["UNLOADED_RADIUS"]

    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)

    K_x = FZ * (PKX1 + PKX2 * df_z) * np.exp(PKX3 * df_z) * LKX

    S_Hx = (PHX1 + PHX2 * df_z) * LHX
    S_Vx = FZ * (PVX1 + PVX2 * df_z) * LVX * LMUX

    M_y = R0 * (S_Vx + K_x * S_Hx)

    return M_y