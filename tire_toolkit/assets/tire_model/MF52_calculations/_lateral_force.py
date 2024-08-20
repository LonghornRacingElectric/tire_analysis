import numpy as np

def get_F_y(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma):
    
    combined_lat = _combined_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma)

    return combined_lat

def _combined_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma) -> float:

    if kappa != 0:
        # Combined lat coeffs
        RBY1 = lat_coeffs["RBY1"]
        RBY2 = lat_coeffs["RBY2"]
        RBY3 = lat_coeffs["RBY3"]
        RCY1 = lat_coeffs["RCY1"]
        REY1 = lat_coeffs["REY1"]
        REY2 = lat_coeffs["REY2"]
        RHY1 = lat_coeffs["RHY1"]
        RHY2 = lat_coeffs["RHY2"]
        RVY1 = lat_coeffs["RVY1"]
        RVY2 = lat_coeffs["RVY2"]
        RVY3 = lat_coeffs["RVY3"]
        RVY4 = lat_coeffs["RVY4"]
        RVY5 = lat_coeffs["RVY5"]
        RVY6 = lat_coeffs["RVY6"]

        # Scaling coeffs
        LFZO = scaling_coeffs["LFZO"]
        LYKA = scaling_coeffs["LYKA"]
        LVYKA = scaling_coeffs["LVYKA"]

        # Vertical coeffs
        FNOMIN = vertical_coeffs["FNOMIN"]

        # Pure slip conditions
        mu_y, FY_0 = _pure_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, kappa, gamma)
        
        df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)

        C_ySR = RCY1
        B_ySR = RBY1 * np.cos(np.arctan(RBY2 * (alpha - RBY3))) * LYKA
        E_ySR = REY1 + REY2 * df_z
        S_HySR = RHY1 + RHY2 * df_z

        D_VySR = mu_y * FZ * (RVY1 + RVY2 * df_z + RVY3 * gamma) * np.cos(np.arctan(RVY4 * alpha))
        S_VySR = D_VySR * np.sin(RVY5 * np.arctan(RVY6 * kappa)) * LVYKA

        SR_s = kappa + S_HySR

        # The Delft paper defines FY_adj this way, but we can decompose this into force output from the pure fit * some scaling coefficient
        # D_ySR = self.pure_lat_coeffs([self.FZ_nom, FZ, SA, IA]) / (np.cos(C_ySR * np.arctan(B_ySR * S_HySR - E_ySR * (B_ySR * S_HySR - np.arctan(B_ySR * S_HySR)))))
        # FY_adj = D_ySR * np.cos(C_ySR * np.arctan(B_ySR * SR_s - E_ySR * (B_ySR * SR_s - np.arctan(B_ySR * SR_s)))) + S_VySR

        # This is the same calculation, but the variables are shown a more intuitive way
        G_ySR = (np.cos(C_ySR * np.arctan(B_ySR * SR_s - E_ySR * (B_ySR * SR_s - np.arctan(B_ySR * SR_s))))) / (np.cos(C_ySR * np.arctan(B_ySR * S_HySR - E_ySR * (B_ySR * S_HySR - np.arctan(B_ySR * S_HySR)))))

        FY_adj = FY_0 * G_ySR + S_VySR
        mu_adj = mu_y * G_ySR

        return [abs(mu_adj), FY_adj]

    else:
        mu_y, FY_0 = _pure_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, gamma)
        return [abs(mu_y), FY_0]

def _pure_lat(lat_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, gamma) -> float:
    # Pure lat coeffs
    PCY1 = lat_coeffs["PCY1"]
    PDY1 = lat_coeffs["PDY1"]
    PDY2 = lat_coeffs["PDY2"]
    PDY3 = lat_coeffs["PDY3"]
    PEY1 = lat_coeffs["PEY1"]
    PEY2 = lat_coeffs["PEY2"]
    PEY3 = lat_coeffs["PEY3"]
    PEY4 = lat_coeffs["PEY4"]
    PKY1 = lat_coeffs["PKY1"]
    PKY2 = lat_coeffs["PKY2"]
    PKY3 = lat_coeffs["PKY3"]
    PHY1 = lat_coeffs["PHY1"]
    PHY2 = lat_coeffs["PHY2"]
    PHY3 = lat_coeffs["PHY3"]
    PVY1 = lat_coeffs["PVY1"]
    PVY2 = lat_coeffs["PVY2"]
    PVY3 = lat_coeffs["PVY3"]
    PVY4 = lat_coeffs["PVY4"]

    # Scaling coeffs
    LFZO = scaling_coeffs["LFZO"]
    LCY = scaling_coeffs["LCY"]
    LMUY = scaling_coeffs["LMUY"]
    LEY = scaling_coeffs["LEY"]
    LKY = scaling_coeffs["LKY"]
    LHY = scaling_coeffs["LHY"]
    LVY = scaling_coeffs["LVY"]
    LGAY = scaling_coeffs["LGAY"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]
    
    IA_y = gamma * LGAY
    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    mu_y = (PDY1 + PDY2 * df_z) * (1 - PDY3 * IA_y**2) * LMUY

    C_y = PCY1 * LCY
    D_y = mu_y * FZ
    K_y = PKY1 * FNOMIN * np.sin(2 * np.arctan(FZ / (PKY2 * FNOMIN * LFZO))) * \
        (1 - PKY3 * abs(IA_y)) * LFZO * LKY
    B_y = K_y / (C_y * D_y)

    S_Hy = (PHY1 + PHY2 * df_z) * LHY + PHY3 * IA_y
    S_Vy = FZ * ((PVY1 + PVY2 * df_z) * LVY + (PVY3 + PVY4 * df_z) * IA_y) * LMUY
    SA_y = alpha + S_Hy

    E_y = (PEY1 + PEY2 * df_z) * (1 - (PEY3 + PEY4 * IA_y) * np.sign(SA_y)) * LEY

    assert E_y <= 1, f"E_y doesn't meet requirements of <= 1 | Curent value: {E_y}"

    F_Y0 = D_y * np.sin(C_y * np.arctan(B_y * SA_y - E_y * (B_y * SA_y - np.arctan(B_y * SA_y)))) + S_Vy
    F_Y = F_Y0

    return mu_y, F_Y
