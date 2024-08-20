import numpy as np

from tire_toolkit.assets.tire_model.MF52_calculations._longitudinal_force import get_F_x
from tire_toolkit.assets.tire_model.MF52_calculations._lateral_force import get_F_y

def get_M_z(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, FZ, alpha, kappa, gamma):
    
    combined_aligning = _combined_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, FZ, alpha, kappa, gamma)

    return combined_aligning

def _combined_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, long_coeffs, vertical_coeffs, dimensions, FZ, alpha, kappa, gamma):
    
    if kappa != 0:
        # Pure aligning coeffs
        QBZ1 = aligning_coeffs["QBZ1"]
        QBZ2 = aligning_coeffs["QBZ2"]
        QBZ3 = aligning_coeffs["QBZ3"]
        QBZ4 = aligning_coeffs["QBZ4"]
        QBZ5 = aligning_coeffs["QBZ5"]
        QBZ9 = aligning_coeffs["QBZ9"]
        QBZ10 = aligning_coeffs["QBZ10"]
        QCZ1 = aligning_coeffs["QCZ1"]
        QDZ1 = aligning_coeffs["QDZ1"]
        QDZ2 = aligning_coeffs["QDZ2"]
        QDZ3 = aligning_coeffs["QDZ3"]
        QDZ4 = aligning_coeffs["QDZ4"]
        QDZ6 = aligning_coeffs["QDZ6"]
        QDZ7 = aligning_coeffs["QDZ7"]
        QDZ8 = aligning_coeffs["QDZ8"]
        QDZ9 = aligning_coeffs["QDZ9"]
        QEZ1 = aligning_coeffs["QEZ1"]
        QEZ2 = aligning_coeffs["QEZ2"]
        QEZ3 = aligning_coeffs["QEZ3"]
        QEZ4 = aligning_coeffs["QEZ4"]
        QEZ5 = aligning_coeffs["QEZ5"]
        QHZ1 = aligning_coeffs["QHZ1"]
        QHZ2 = aligning_coeffs["QHZ2"]
        QHZ3 = aligning_coeffs["QHZ3"]
        QHZ4 = aligning_coeffs["QHZ4"]

        # Combined aligning coeffs
        SSZ1 = aligning_coeffs["SSZ1"]
        SSZ2 = aligning_coeffs["SSZ2"]
        SSZ3 = aligning_coeffs["SSZ3"]
        SSZ4 = aligning_coeffs["SSZ4"]

        # Pure lat coeffs
        PCY1 = lat_coeffs["PCY1"]
        PDY1 = lat_coeffs["PDY1"]
        PDY2 = lat_coeffs["PDY2"]
        PDY3 = lat_coeffs["PDY3"]
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

        # Combined lat coeffs
        RVY1 = lat_coeffs["RVY1"]
        RVY2 = lat_coeffs["RVY2"]
        RVY3 = lat_coeffs["RVY3"]
        RVY4 = lat_coeffs["RVY4"]
        RVY5 = lat_coeffs["RVY5"]
        RVY6 = lat_coeffs["RVY6"]

        # Pure long coeffs
        PKX1 = long_coeffs["PKX1"]
        PKX2 = long_coeffs["PKX2"]
        PKX3 = long_coeffs["PKX3"]

        # Scaling coeffs
        LFZO = scaling_coeffs["LFZO"]
        LKX = scaling_coeffs["LKX"]
        LCY = scaling_coeffs["LCY"]
        LMUY = scaling_coeffs["LMUY"]
        LKY = scaling_coeffs["LKY"]
        LHY = scaling_coeffs["LHY"]
        LVY = scaling_coeffs["LVY"]
        LGAY = scaling_coeffs["LGAY"]
        LTR = scaling_coeffs["LTR"]
        LRES = scaling_coeffs["LRES"]
        LGAZ = scaling_coeffs["LGAZ"]
        LVYKA = scaling_coeffs["LVYKA"]
        LS = scaling_coeffs["LS"]

        # Vertical coeffs
        FNOMIN = vertical_coeffs["FNOMIN"]

        # Dimension coeffs
        R0 = dimensions["UNLOADED_RADIUS"]

        df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
        IA_z = gamma * LGAZ

        # Lateral Dependencies
        IA_y = gamma * LGAY
        mu_y = (PDY1 + PDY2 * df_z) * (1 - PDY3 * IA_y**2) * LMUY
        S_Hy = (PHY1 + PHY2 * df_z) * LHY + PHY3 * IA_y
        S_Vy = FZ * ((PVY1 + PVY2 * df_z) * LVY + (PVY3 + PVY4 * df_z) * IA_y) * LMUY

        K_y = PKY1 * FNOMIN * np.sin(2 * np.arctan(FZ / (PKY2 * FNOMIN * LFZO))) * \
            (1 - PKY3 * abs(IA_y)) * LFZO * LKY
        C_y = PCY1 * LCY
        D_y = mu_y * FZ
        B_y = K_y / (C_y * D_y)

        # Longitudinal
        K_x = FZ * (PKX1 + PKX2 * df_z) * np.exp(PKX3 * df_z) * LKX
        
        # Aligning
        D_t = FZ * (QDZ1 + QDZ2 * df_z) * (1 + QDZ3 * IA_z + QDZ4 * IA_z**2) * (R0 / FNOMIN) * LTR
        C_t = QCZ1
        B_t = (QBZ1 + QBZ2 * df_z + QBZ3 * df_z**2) * (1 + QBZ4 * IA_z + QBZ5 * abs(IA_z)) * LKY / LMUY
        
        D_VySR = mu_y * FZ * (RVY1 + RVY2 * df_z + RVY3 * gamma) * np.cos(np.arctan(RVY4 * alpha))
        S_VySR = D_VySR * np.sin(RVY5 * np.arctan(RVY6 * kappa)) * LVYKA

        S_Ht = QHZ1 + QHZ2 * df_z + (QHZ3 + QHZ4 * df_z) * IA_z
        SA_t = alpha + S_Ht

        E_t = (QEZ1 + QEZ2 * df_z + QEZ3 * df_z**2) * (1 + (QEZ4 + QEZ5 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))

        S_Hf = S_Hy + S_Vy / K_y
        SA_r = alpha + S_Hf

        # Adjusted SA values
        SA_t_eq = np.arctan(np.sqrt((np.tan(SA_t))**2 + (K_x / K_y)**2 * kappa**2)) * np.sign(SA_t)
        SA_r_eq = np.arctan(np.sqrt((np.tan(SA_r))**2 + (K_x / K_y)**2 * kappa**2)) * np.sign(SA_r)

        # Pneumatic trail
        t_adj = D_t * np.cos(C_t * np.arctan(B_t * SA_t_eq - E_t * (B_t * SA_t_eq - np.arctan(B_t * SA_t_eq)))) * np.cos(alpha)

        # Pure lat and long forces
        FX = get_F_x(
            long_coeffs,
            scaling_coeffs, 
            vertical_coeffs,
            FZ,
            alpha,
            kappa,
            gamma)[1]
        
        FY = get_F_y(
            lat_coeffs = lat_coeffs,
            scaling_coeffs = scaling_coeffs,
            vertical_coeffs = vertical_coeffs,
            FZ = FZ,
            alpha = alpha,
            kappa = kappa,
            gamma = 0)[1]

        F_y_IA_adj = FY - S_VySR

        D_r = FZ * ((QDZ6 + QDZ7 * df_z) * LRES + (QDZ8 + QDZ9 * df_z) * IA_z) * R0 * LMUY
        B_r = QBZ9 * LKY / LMUY + QBZ10 * B_y * C_y

        M_zr = D_r * np.cos(np.arctan(B_r * SA_r_eq)) * np.cos(alpha)

        s = (SSZ1 + SSZ2 * (FY / FNOMIN) + (SSZ3 + SSZ4 * df_z) * gamma) * R0 * LS

        MZ_adj = -t_adj * F_y_IA_adj + M_zr + s * FX

        return [s, t_adj, MZ_adj]

    else:
        t, Mz = _pure_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, vertical_coeffs, dimensions, FZ, alpha, kappa, gamma)
        return [0, t, Mz]

def _pure_aligning(aligning_coeffs, scaling_coeffs, lat_coeffs, vertical_coeffs, dimensions, FZ, alpha, kappa, gamma) -> float:
    
    # Pure aligning coeffs
    QBZ1 = aligning_coeffs["QBZ1"]
    QBZ2 = aligning_coeffs["QBZ2"]
    QBZ3 = aligning_coeffs["QBZ3"]
    QBZ4 = aligning_coeffs["QBZ4"]
    QBZ5 = aligning_coeffs["QBZ5"]
    QBZ9 = aligning_coeffs["QBZ9"]
    QBZ10 = aligning_coeffs["QBZ10"]
    QCZ1 = aligning_coeffs["QCZ1"]
    QDZ1 = aligning_coeffs["QDZ1"]
    QDZ2 = aligning_coeffs["QDZ2"]
    QDZ3 = aligning_coeffs["QDZ3"]
    QDZ4 = aligning_coeffs["QDZ4"]
    QDZ6 = aligning_coeffs["QDZ6"]
    QDZ7 = aligning_coeffs["QDZ7"]
    QDZ8 = aligning_coeffs["QDZ8"]
    QDZ9 = aligning_coeffs["QDZ9"]
    QEZ1 = aligning_coeffs["QEZ1"]
    QEZ2 = aligning_coeffs["QEZ2"]
    QEZ3 = aligning_coeffs["QEZ3"]
    QEZ4 = aligning_coeffs["QEZ4"]
    QEZ5 = aligning_coeffs["QEZ5"]
    QHZ1 = aligning_coeffs["QHZ1"]
    QHZ2 = aligning_coeffs["QHZ2"]
    QHZ3 = aligning_coeffs["QHZ3"]
    QHZ4 = aligning_coeffs["QHZ4"]

    # Pure lat coeffs
    PCY1 = lat_coeffs["PCY1"]
    PDY1 = lat_coeffs["PDY1"]
    PDY2 = lat_coeffs["PDY2"]
    PDY3 = lat_coeffs["PDY3"]
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
    LKY = scaling_coeffs["LKY"]
    LHY = scaling_coeffs["LHY"]
    LVY = scaling_coeffs["LVY"]
    LGAY = scaling_coeffs["LGAY"]
    LTR = scaling_coeffs["LTR"]
    LRES = scaling_coeffs["LRES"]
    LGAZ = scaling_coeffs["LGAZ"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]

    # Dimension coeffs
    R0 = dimensions["UNLOADED_RADIUS"]

    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    IA_z = gamma * LGAZ

    # Lateral Dependencies
    IA_y = gamma * LGAY
    mu_y = (PDY1 + PDY2 * df_z) * (1 - PDY3 * IA_y**2) * LMUY
    S_Hy = (PHY1 + PHY2 * df_z) * LHY + PHY3 * IA_y
    S_Vy = FZ * ((PVY1 + PVY2 * df_z) * LVY + (PVY3 + PVY4 * df_z) * IA_y) * LMUY

    K_y = PKY1 * FNOMIN * np.sin(2 * np.arctan(FZ / (PKY2 * FNOMIN * LFZO))) * \
        (1 - PKY3 * abs(IA_y)) * LFZO * LKY
    C_y = PCY1 * LCY
    D_y = mu_y * FZ
    B_y = K_y / (C_y * D_y)

    # Pure Aligning Moment
    S_Ht = QHZ1 + QHZ2 * df_z + (QHZ3 + QHZ4 * df_z) * IA_z
    SA_t = alpha + S_Ht

    D_r = FZ * ((QDZ6 + QDZ7 * df_z) * LRES + (QDZ8 + QDZ9 * df_z) * IA_z) * R0 * LMUY
    B_r = QBZ9 * LKY / LMUY + QBZ10 * B_y * C_y

    D_t = FZ * (QDZ1 + QDZ2 * df_z) * (1 + QDZ3 * IA_z + QDZ4 * IA_z**2) * (R0 / FNOMIN) * LTR
    C_t = QCZ1
    B_t = (QBZ1 + QBZ2 * df_z + QBZ3 * df_z**2) * (1 + QBZ4 * IA_z + QBZ5 * abs(IA_z)) * LKY / LMUY

    E_t = (QEZ1 + QEZ2 * df_z + QEZ3 * df_z**2) * (1 + (QEZ4 + QEZ5 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))

    assert E_t <= 1, f"E_t doesn't meet requirements of <= 1 | Curent value: {E_t}"
    
    # Residual Torque

    # Found typo in the delft paper. Hf should be Hr
    S_Hr = S_Hy + S_Vy / K_y

    SA_r = alpha + S_Hr
    M_zr = D_r * np.cos(np.arctan(B_r * SA_r)) * np.cos(alpha)

    # Pneumatic trail
    t = D_t * np.cos(C_t * np.arctan(B_t * SA_t - E_t * (B_t * SA_t - np.arctan(B_t * SA_t)))) * np.cos(alpha)

    FY_0 = get_F_y(
        lat_coeffs = lat_coeffs,
        scaling_coeffs = scaling_coeffs,
        vertical_coeffs = vertical_coeffs,
        FZ = FZ,
        alpha = alpha,
        kappa = kappa,
        gamma = gamma)[1]
    
    M_Z0 = -1 * t * FY_0 + M_zr
    M_Z = M_Z0

    return t, M_Z