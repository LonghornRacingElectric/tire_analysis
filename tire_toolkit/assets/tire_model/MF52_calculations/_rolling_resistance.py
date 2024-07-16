import numpy as np

def get_M_y(rolling_coeffs, long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma):

    # Rolling resistance coeffs
    QSY1 = rolling_coeffs["QSY1"]
    QSY2 = rolling_coeffs["QSY2"]
    QSY3 = rolling_coeffs["QSY3"]
    QSY4 = rolling_coeffs["QSY4"]

    # Pure long coeffs
    CFX1 = long_coeffs["PCX1"]
    CFX2 = long_coeffs["PDX1"]
    CFX3 = long_coeffs["PDX2"]
    CFX4 = long_coeffs["PDX3"]
    CFX5 = long_coeffs["PEX1"]
    CFX6 = long_coeffs["PEX2"]
    CFX7 = long_coeffs["PEX3"]
    CFX8 = long_coeffs["PEX4"]
    CFX9 = long_coeffs["PKX1"]
    CFX10 = long_coeffs["PKX2"]
    CFX11 = long_coeffs["PKX3"]
    CFX12 = long_coeffs["PHX1"]
    CFX13 = long_coeffs["PHX2"]
    CFX14 = long_coeffs["PVX1"]
    CFX15 = long_coeffs["PVX2"]

    # Combined long coeffs
    CCFX1 = long_coeffs["RBX1"]
    CCFX2 = long_coeffs["RBX2"]
    CCFX3 = long_coeffs["RCX1"]
    CCFX4 = long_coeffs["REX1"]
    CCFX5 = long_coeffs["REX2"]
    CCFX6 = long_coeffs["RHX1"]

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

    IA_x = gamma * LGAX
    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    mu_x = (CFX2 + CFX3 * df_z) * (1 - CFX4 * IA_x**2) * LMUX

    C_x = CFX1 * LCX
    D_x = mu_x * FZ
    K_x = FZ * (CFX9 + CFX10 * df_z) * np.exp(CFX11 * df_z) * LKX
    B_x = K_x / (C_x * D_x)

    S_Hx = (CFX12 + CFX13 * df_z) * LHX
    S_Vx = FZ * (CFX14 + CFX15 * df_z) * LVX * LMUX

    M_y = R0 * (S_Vx + K_x * S_Hx)

    return M_y