import numpy as np
import warnings
warnings.filterwarnings("ignore")

def get_F_x(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma):
    combined_long = _combined_long(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma)

    return combined_long

def _combined_long(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, alpha, kappa, gamma) -> float:

    if alpha != 0:
        # Combined long coeffs
        RBX1 = long_coeffs["RBX1"]
        RBX2 = long_coeffs["RBX2"]
        RCX1 = long_coeffs["RCX1"]
        REX1 = long_coeffs["REX1"]
        REX2 = long_coeffs["REX2"]
        RHX1 = long_coeffs["RHX1"]

        # Scaling coeffs
        LFZO = scaling_coeffs["LFZO"]
        LXAL = scaling_coeffs["LXAL"]

        # Vertical coeffs
        FNOMIN = vertical_coeffs["FNOMIN"]
        
        df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
        
        C_xSA = RCX1
        B_xSA = RBX1 * np.cos(np.arctan(RBX2 * kappa)) * LXAL
        E_xSA = REX1 + REX2 * df_z
        S_HxSA = RHX1

        SA_s = alpha + S_HxSA

        # The Delft paper defines FX_adj this way, but we can decompose this into force output from the pure fit * some scaling coefficient
        # D_xSA = _pure_long([FZ, SR, IA]) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
        # FX_adj = D_xSA * np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))

        # This is the same calculation, but the variables are shown a more intuitive way
        G_xSA = (np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
        mu_x, FX_0 = _pure_long(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, kappa, gamma)

        FX_adj = FX_0 * G_xSA
        mu_x_adj = mu_x * G_xSA
    
        return mu_x_adj, FX_adj

    else:
        return _pure_long(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, kappa, gamma)

def _pure_long(long_coeffs, scaling_coeffs, vertical_coeffs, FZ, kappa, gamma) -> float:

    # Pure long coeffs
    PCX1 = long_coeffs["PCX1"]
    PDX1 = long_coeffs["PDX1"]
    PDX2 = long_coeffs["PDX2"]
    PDX3 = long_coeffs["PDX3"]
    PEX1 = long_coeffs["PEX1"]
    PEX2 = long_coeffs["PEX2"]
    PEX3 = long_coeffs["PEX3"]
    PEX4 = long_coeffs["PEX4"]
    PKX1 = long_coeffs["PKX1"]
    PKX2 = long_coeffs["PKX2"]
    PKX3 = long_coeffs["PKX3"]
    PHX1 = long_coeffs["PHX1"]
    PHX2 = long_coeffs["PHX2"]
    PVX1 = long_coeffs["PVX1"]
    PVX2 = long_coeffs["PVX2"]

    # Scaling coeffs
    LFZO = scaling_coeffs["LFZO"]
    LCX = scaling_coeffs["LCX"]
    LMUX = scaling_coeffs["LMUX"]
    LEX = scaling_coeffs["LEX"]
    LKX = scaling_coeffs["LKX"]
    LHX = scaling_coeffs["LHX"]
    LVX = scaling_coeffs["LVX"]
    LGAX = scaling_coeffs["LGAX"]

    # Vertical coeffs
    FNOMIN = vertical_coeffs["FNOMIN"]

    IA_x = gamma * LGAX
    df_z = (FZ - FNOMIN * LFZO) / (FNOMIN * LFZO)
    mu_x = (PDX1 + PDX2 * df_z) * (1 - PDX3 * IA_x**2) * LMUX

    C_x = PCX1 * LCX
    D_x = mu_x * FZ
    K_x = FZ * (PKX1 + PKX2 * df_z) * np.exp(PKX3 * df_z) * LKX
    B_x = K_x / (C_x * D_x)

    S_Hx = (PHX1 + PHX2 * df_z) * LHX
    S_Vx = FZ * (PVX1 + PVX2 * df_z) * LVX * LMUX
    SR_x = kappa + S_Hx

    E_x = (PEX1 + PEX2 * df_z + PEX3 * df_z**2) * (1 - PEX4 * np.sign(SR_x)) * LEX

    assert E_x <= 1, f"E_x doesn't meet requirements of <= 1 | Curent value: {E_x}"
    
    F_X0 = D_x * np.sin(C_x * np.arctan(B_x * SR_x - E_x * (B_x * SR_x - np.arctan(B_x * SR_x)))) + S_Vx
    F_X = F_X0
    
    return mu_x, F_X