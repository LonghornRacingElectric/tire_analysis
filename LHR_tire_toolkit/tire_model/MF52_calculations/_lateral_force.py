import numpy as np

def get_F_y(lat_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma):
    
    combined_lat = _combined_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma)

    return combined_lat

def _combined_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma) -> float:

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
    IA_y = gamma * LGAY
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * LMUY

    C_ySR = CCFY4
    B_ySR = CCFY1 * np.cos(np.arctan(CCFY2 * (alpha - CCFY3))) * LYKA
    E_ySR = CCFY5 + CCFY6 * df_z
    S_HySR = CCFY7 + CCFY8 * df_z

    D_VySR = mu_y * FZ * (CCFY9 + CCFY10 * df_z + CCFY11 * gamma) * np.cos(np.arctan(CCFY12 * alpha))

    S_VySR = D_VySR * np.sin(CCFY13 * np.arctan(CCFY14 * kappa)) * LVYKA

    SR_s = kappa + S_HySR
        
    # The Delft paper defines FY_adj this way, but we can decompose this into force output from the pure fit * some scaling coefficient
    # D_ySR = self.pure_lat_coeffs([self.FZ_nom, FZ, SA, IA]) / (np.cos(C_ySR * np.arctan(B_ySR * S_HySR - E_ySR * (B_ySR * S_HySR - np.arctan(B_ySR * S_HySR)))))
    # FY_adj = D_ySR * np.cos(C_ySR * np.arctan(B_ySR * SR_s - E_ySR * (B_ySR * SR_s - np.arctan(B_ySR * SR_s)))) + S_VySR

    # This is the same calculation, but the variables are shown a more intuitive way
    G_ySR = (np.cos(C_ySR * np.arctan(B_ySR * SR_s - E_ySR * (B_ySR * SR_s - np.arctan(B_ySR * SR_s))))) / (np.cos(C_ySR * np.arctan(B_ySR * S_HySR - E_ySR * (B_ySR * S_HySR - np.arctan(B_ySR * S_HySR)))))
    FY_0 = _pure_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma)

    FY_adj = FY_0 * G_ySR + S_VySR

    return [abs(mu_y), FY_adj]

def _pure_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma) -> float:
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
    
    IA_y = gamma * LGAY
    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * LMUY

    C_y = CFY1 * LCY
    D_y = mu_y * FZ
    K_y = CFY9 * FNOMIN * np.sin(2 * np.arctan(FZ / (CFY10 * FNOMIN * LFZO))) * \
        (1 - CFY11 * abs(IA_y)) * LFZO * LKY
    B_y = K_y / (C_y * D_y)

    S_Hy = (CFY12 + CFY13 * df_z) * LHY + CFY14 * IA_y
    S_Vy = FZ * ((CFY15 + CFY16 * df_z) * LVY + (CFY17 + CFY18 * df_z) * IA_y) * LMUY
    SA_y = alpha + S_Hy

    E_y = (CFY5 + CFY6 * df_z) * (1 - (CFY7 + CFY8 * IA_y) * np.sign(SA_y)) * LEY

    F_Y0 = D_y * np.sin(C_y * np.arctan(B_y * SA_y - E_y * (B_y * SA_y - np.arctan(B_y * SA_y)))) + S_Vy
    F_Y = F_Y0

    return F_Y
