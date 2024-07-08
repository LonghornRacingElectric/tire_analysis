import numpy as np
import warnings
warnings.filterwarnings("ignore")

def get_F_x(long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma):
    combined_long = _combined_long(long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma)

    return combined_long

def _combined_long(long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma) -> float:

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
    
    C_xSA = CCFX3
    B_xSA = CCFX1 * np.cos(np.arctan(CCFX2 * kappa)) * LXAL
    E_xSA = CCFX4 + CCFX5 * df_z
    S_HxSA = CCFX6

    SA_s = alpha + S_HxSA

    # The Delft paper defines FX_adj this way, but we can decompose this into force output from the pure fit * some scaling coefficient
    # D_xSA = _pure_long([FZ, SR, IA]) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
    # FX_adj = D_xSA * np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))

    # This is the same calculation, but the variables are shown a more intuitive way
    G_xSA = (np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
    mu_x, FX_0 = _pure_long(long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma)

    FX_adj = FX_0 * G_xSA
    
    return mu_x, FX_adj

def _pure_long(long_coeffs, scaling_coeffs, vertical_coeffs, dimensions, operating_conditions, FZ, alpha, kappa, gamma) -> float:

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
    SR_x = kappa + S_Hx

    E_x = (CFX5 + CFX6 * df_z + CFX7 * df_z**2) * (1 - CFX8 * np.sign(kappa)) * LEX
    
    F_X0 = D_x * np.sin(C_x * np.arctan(B_x * SR_x - E_x * (B_x * SR_x - np.arctan(B_x * SR_x)))) + S_Vx
    F_X = F_X0
    
    return mu_x, F_X