import numpy as np

from LHR_tire_toolkit.MF52_calculations._longitudinal_force import get_F_x
from LHR_tire_toolkit.MF52_calculations._lateral_force import get_F_y

def get_M_z(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma):
    
    combined_aligning = _combined_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma)

    return combined_aligning

def _combined_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma):
    
    # Pure aligning coeffs
    CMZ1 = aligning_coeffs["QBZ1"]
    CMZ2 = aligning_coeffs["QBZ2"]
    CMZ3 = aligning_coeffs["QBZ3"]
    CMZ4 = aligning_coeffs["QBZ4"]
    CMZ5 = aligning_coeffs["QBZ5"]
    CMZ6 = aligning_coeffs["QBZ9"]
    CMZ7 = aligning_coeffs["QBZ10"]
    CMZ8 = aligning_coeffs["QCZ1"]
    CMZ9 = aligning_coeffs["QDZ1"]
    CMZ10 = aligning_coeffs["QDZ2"]
    CMZ11 = aligning_coeffs["QDZ3"]
    CMZ12 = aligning_coeffs["QDZ4"]
    CMZ13 = aligning_coeffs["QDZ6"]
    CMZ14 = aligning_coeffs["QDZ7"]
    CMZ15 = aligning_coeffs["QDZ8"]
    CMZ16 = aligning_coeffs["QDZ9"]
    CMZ17 = aligning_coeffs["QEZ1"]
    CMZ18 = aligning_coeffs["QEZ2"]
    CMZ19 = aligning_coeffs["QEZ3"]
    CMZ20 = aligning_coeffs["QEZ4"]
    CMZ21 = aligning_coeffs["QEZ5"]
    CMZ22 = aligning_coeffs["QHZ1"]
    CMZ23 = aligning_coeffs["QHZ2"]
    CMZ24 = aligning_coeffs["QHZ3"]
    CMZ25 = aligning_coeffs["QHZ4"]

    # Combined aligning coeffs
    CCMZ1 = aligning_coeffs["SSZ1"]
    CCMZ2 = aligning_coeffs["SSZ2"]
    CCMZ3 = aligning_coeffs["SSZ3"]
    CCMZ4 = aligning_coeffs["SSZ4"]

    # Pure lat coeffs
    CFY1 = lat_coeffs["PCY1"]
    CFY2 = lat_coeffs["PDY1"]
    CFY3 = lat_coeffs["PDY2"]
    CFY4 = lat_coeffs["PDY3"]
    CFY5 = lat_coeffs["PEY1"]
    CFY6 = lat_coeffs["PEY2"]
    CFY7 = lat_coeffs["PEY3"]
    CFY8 = lat_coeffs["PEY4"]
    CFY9 = lat_coeffs["PKY1"]
    CFY10 = lat_coeffs["PKY2"]
    CFY11 = lat_coeffs["PKY3"]
    CFY12 = lat_coeffs["PHY1"]
    CFY13 = lat_coeffs["PHY2"]
    CFY14 = lat_coeffs["PHY3"]
    CFY15 = lat_coeffs["PVY1"]
    CFY16 = lat_coeffs["PVY2"]
    CFY17 = lat_coeffs["PVY3"]
    CFY18 = lat_coeffs["PVY4"]

    # Combined lat coeffs
    CCFY1 = lat_coeffs["RBY1"]
    CCFY2 = lat_coeffs["RBY2"]
    CCFY3 = lat_coeffs["RBY3"]
    CCFY4 = lat_coeffs["RCY1"]
    CCFY5 = lat_coeffs["REY1"]
    CCFY6 = lat_coeffs["REY2"]
    CCFY7 = lat_coeffs["RHY1"]
    CCFY8 = lat_coeffs["RHY2"]
    CCFY9 = lat_coeffs["RVY1"]
    CCFY10 = lat_coeffs["RVY2"]
    CCFY11 = lat_coeffs["RVY3"]
    CCFY12 = lat_coeffs["RVY4"]
    CCFY13 = lat_coeffs["RVY5"]
    CCFY14 = lat_coeffs["RVY6"]

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

    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    IA_z = gamma * LGAZ

    ### Pure slip dependencies

    # Lateral
    IA_y = gamma * LGAY
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * LMUY
    S_Hy = (CFY12 + CFY13 * df_z) * LHY + CFY14 * IA_y
    S_Vy = FZ * ((CFY15 + CFY16 * df_z) * LVY + (CFY17 + CFY18 * df_z) * IA_y) * LMUY

    K_y = CFY9 * FNOMIN * np.sin(2 * np.arctan(FZ / (CFY10 * FNOMIN * LFZO))) * \
        (1 - CFY11 * abs(IA_y)) * LFZO * LKY
    C_y = CFY1 * LCY
    D_y = mu_y * FZ
    B_y = K_y / (C_y * D_y)

    # Longitudinal
    K_x = FZ * (CFX9 + CFX10 * df_z) * np.exp(CFX11 * df_z) * LKX
    
    # Aligning
    D_t = FZ * (CMZ9 + CMZ10 * df_z) * (1 + CMZ11 * IA_z + CMZ12 * IA_z**2) * (R0 / FNOMIN) * LTR
    C_t = CMZ8
    B_t = (CMZ1 + CMZ2 * df_z + CMZ3 * df_z**2) * (1 + CMZ4 * IA_z + CMZ5 * abs(IA_z)) * LKY / LMUY
    
    D_VySR = mu_y * FZ * (CCFY9 + CCFY10 * df_z + CCFY11 * gamma) * np.cos(np.arctan(CCFY12 * alpha))
    S_VySR = D_VySR * np.sin(CCFY13 * np.arctan(CCFY14 * kappa)) * LVYKA

    S_Ht = CMZ22 + CMZ23 * df_z + (CMZ24 + CMZ25 * df_z) * IA_z
    SA_t = alpha + S_Ht

    E_t = (CMZ17 + CMZ18 * df_z + CMZ19 * df_z**2) * (1 + (CMZ20 + CMZ21 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))

    S_Hf = S_Hy + S_Vy / K_y
    SA_r = alpha + S_Hf

    # Adjusted SA values
    SA_t_eq = np.arctan(np.sqrt((np.tan(SA_t))**2 + (K_x / K_y)**2 * kappa**2)) * np.sign(SA_t)
    SA_r_eq = np.arctan(np.sqrt((np.tan(SA_r))**2 + (K_x / K_y)**2 * kappa**2)) * np.sign(SA_r)

    # Pneumatic trail
    t_adj = D_t * np.cos(C_t * np.arctan(B_t * SA_t_eq - E_t * (B_t * SA_t_eq - np.arctan(B_t * SA_t_eq)))) * np.cos(alpha)

    # FX = self._combined_long([FZ, SA, SR, IA])
    FX = get_F_x(
        long_coeffs,
        scaling_coeffs, 
        vertical_coeffs,
        dimensions,
        operating_conditions,
        FZ,
        alpha,
        kappa,
        gamma)[1]
    
    FY = get_F_y(
        lat_coeffs = lat_coeffs,
        scaling_coeffs = scaling_coeffs,
        vertical_coeffs = vertical_coeffs,
        dimensions = dimensions,
        operating_conditions = operating_conditions,
        FZ = FZ,
        alpha = alpha,
        kappa = kappa,
        gamma = gamma)[1]

    F_y_IA_adj = FY - S_VySR

    D_r = FZ * ((CMZ13 + CMZ14 * df_z) * LRES + (CMZ15 + CMZ16 * df_z) * IA_z) * R0 * LMUY
    B_r = CMZ6 * LKY / LMUY + CMZ7 * B_y * C_y

    M_zr = D_r * np.cos(np.arctan(B_r * SA_r_eq)) * np.cos(alpha)

    s = (CCMZ1 + CCMZ2 * (FY / FNOMIN) + (CCMZ3 + CCMZ4 * df_z) * gamma) * R0 * LS

    MZ_adj = -t_adj * F_y_IA_adj + M_zr + s * FX

    return [t_adj, MZ_adj]

def _pure_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma) -> float:
    
    # Pure aligning coeffs
    CMZ1 = aligning_coeffs["QBZ1"]
    CMZ2 = aligning_coeffs["QBZ2"]
    CMZ3 = aligning_coeffs["QBZ3"]
    CMZ4 = aligning_coeffs["QBZ4"]
    CMZ5 = aligning_coeffs["QBZ5"]
    CMZ6 = aligning_coeffs["QBZ9"]
    CMZ7 = aligning_coeffs["QBZ10"]
    CMZ8 = aligning_coeffs["QCZ1"]
    CMZ9 = aligning_coeffs["QDZ1"]
    CMZ10 = aligning_coeffs["QDZ2"]
    CMZ11 = aligning_coeffs["QDZ3"]
    CMZ12 = aligning_coeffs["QDZ4"]
    CMZ13 = aligning_coeffs["QDZ6"]
    CMZ14 = aligning_coeffs["QDZ7"]
    CMZ15 = aligning_coeffs["QDZ8"]
    CMZ16 = aligning_coeffs["QDZ9"]
    CMZ17 = aligning_coeffs["QEZ1"]
    CMZ18 = aligning_coeffs["QEZ2"]
    CMZ19 = aligning_coeffs["QEZ3"]
    CMZ20 = aligning_coeffs["QEZ4"]
    CMZ21 = aligning_coeffs["QEZ5"]
    CMZ22 = aligning_coeffs["QHZ1"]
    CMZ23 = aligning_coeffs["QHZ2"]
    CMZ24 = aligning_coeffs["QHZ3"]
    CMZ25 = aligning_coeffs["QHZ4"]

    # Pure lat coeffs
    CFY1 = lat_coeffs["PCY1"]
    CFY2 = lat_coeffs["PDY1"]
    CFY3 = lat_coeffs["PDY2"]
    CFY4 = lat_coeffs["PDY3"]
    CFY5 = lat_coeffs["PEY1"]
    CFY6 = lat_coeffs["PEY2"]
    CFY7 = lat_coeffs["PEY3"]
    CFY8 = lat_coeffs["PEY4"]
    CFY9 = lat_coeffs["PKY1"]
    CFY10 = lat_coeffs["PKY2"]
    CFY11 = lat_coeffs["PKY3"]
    CFY12 = lat_coeffs["PHY1"]
    CFY13 = lat_coeffs["PHY2"]
    CFY14 = lat_coeffs["PHY3"]
    CFY15 = lat_coeffs["PVY1"]
    CFY16 = lat_coeffs["PVY2"]
    CFY17 = lat_coeffs["PVY3"]
    CFY18 = lat_coeffs["PVY4"]

    # Combined lat coeffs
    CCFY1 = lat_coeffs["RBY1"]
    CCFY2 = lat_coeffs["RBY2"]
    CCFY3 = lat_coeffs["RBY3"]
    CCFY4 = lat_coeffs["RCY1"]
    CCFY5 = lat_coeffs["REY1"]
    CCFY6 = lat_coeffs["REY2"]
    CCFY7 = lat_coeffs["RHY1"]
    CCFY8 = lat_coeffs["RHY2"]
    CCFY9 = lat_coeffs["RVY1"]
    CCFY10 = lat_coeffs["RVY2"]
    CCFY11 = lat_coeffs["RVY3"]
    CCFY12 = lat_coeffs["RVY4"]
    CCFY13 = lat_coeffs["RVY5"]
    CCFY14 = lat_coeffs["RVY6"]

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

    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    IA_z = gamma * LGAZ

    # Lateral Dependencies
    IA_y = gamma * LGAY
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * LMUY
    S_Hy = (CFY12 + CFY13 * df_z) * LHY + CFY14 * IA_y
    S_Vy = FZ * ((CFY15 + CFY16 * df_z) * LVY + (CFY17 + CFY18 * df_z) * IA_y) * LMUY

    K_y = CFY9 * FNOMIN * np.sin(2 * np.arctan(FZ / (CFY10 * FNOMIN * LFZO))) * \
        (1 - CFY11 * abs(IA_y)) * LFZO * LKY
    C_y = CFY1 * LCY
    D_y = mu_y * FZ
    B_y = K_y / (C_y * D_y)

    # Pure Aligning Moment
    S_Ht = CMZ22 + CMZ23 * df_z + (CMZ24 + CMZ25 * df_z) * IA_z
    SA_t = alpha + S_Ht

    D_r = FZ * ((CMZ13 + CMZ14 * df_z) * LRES + (CMZ15 + CMZ16 * df_z) * IA_z) * R0 * LMUY
    B_r = CMZ6 * LKY / LMUY + CMZ7 * B_y * C_y

    D_t = FZ * (CMZ9 + CMZ10 * df_z) * (1 + CMZ11 * IA_z + CMZ12 * IA_z**2) * (R0 / FNOMIN) * LTR
    C_t = CMZ8
    B_t = (CMZ1 + CMZ2 * df_z + CMZ3 * df_z**2) * (1 + CMZ4 * IA_z + CMZ5 * abs(IA_z)) * LKY / LMUY

    E_t = (CMZ17 + CMZ18 * df_z + CMZ19 * df_z**2) * (1 + (CMZ20 + CMZ21 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))

    S_Hf = S_Hy + S_Vy / K_y
    
    # Residual Torque

    # Found typo in the delft paper. Hr should be Hf
    SA_r = alpha + S_Hf
    M_zr = D_r * np.cos(np.arctan(B_r * SA_r)) * np.cos(alpha)

    # Pneumatic trail
    t = D_t * np.cos(C_t * np.arctan(B_t * SA_t - E_t * (B_t * SA_t - np.arctan(B_t * SA_t)))) * np.cos(alpha)

    FY = get_F_y(
        lat_coeffs = lat_coeffs,
        scaling_coeffs = scaling_coeffs,
        vertical_coeffs = vertical_coeffs,
        dimensions = dimensions,
        operating_conditions = operating_conditions,
        FZ = FZ,
        alpha = alpha,
        kappa = kappa,
        gamma = gamma)[1]
    
    M_Z0 = -t * FY + M_zr
    M_Z = M_Z0

    return M_Z