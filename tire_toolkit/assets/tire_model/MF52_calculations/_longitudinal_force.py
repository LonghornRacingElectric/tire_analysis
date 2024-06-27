def _combined_long(self, data: list[float]) -> float:
    FZ, SA, SR, IA = data

    [CCFX1, CCFX2, CCFX3, CCFX4, CCFX5, CCFX6] = self.combined_long_coeffs

    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])
    
    C_xSA = CCFX3
    B_xSA = CCFX1 * np.cos(np.arctan(CCFX2 * SR)) * self.scaling_coeffs[21]
    E_xSA = CCFX4 + CCFX5 * df_z
    S_HxSA = CCFX6

    SA_s = SA + S_HxSA

    # The Delft paper defines FX_adj this way, but we can decompose this into force output from the pure fit * some scaling coefficient
    # D_xSA = self._pure_long([FZ, SR, IA]) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
    # FX_adj = D_xSA * np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))

    # This is the same calculation, but the variables are shown a more intuitive way
    G_xSA = (np.cos(C_xSA * np.arctan(B_xSA * SA_s - E_xSA * (B_xSA * SA_s - np.arctan(B_xSA * SA_s))))) / (np.cos(C_xSA * np.arctan(B_xSA * S_HxSA - E_xSA * (B_xSA * S_HxSA - np.arctan(B_xSA * S_HxSA)))))
    FX_0 = self._pure_long([FZ, SR, IA])

    FX_adj = FX_0 * G_xSA
    
    return FX_adj

def _pure_long(self, data: list[float], sign_condition: bool = False) -> float:
    FZ, SR, IA = data

    [CFX1, CFX2, CFX3, CFX4, CFX5, CFX6, CFX7, CFX8, CFX9, CFX10, CFX11, CFX12, CFX13, CFX14, CFX15] = self.pure_long_coeffs

    IA_x = IA * self.scaling_coeffs[7]
    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])
    mu_x = (CFX2 + CFX3 * df_z) * (1 - CFX4 * IA_x**2) * self.scaling_coeffs[2]

    self.C_x = CFX1 * self.scaling_coeffs[1]
    self.D_x = mu_x * FZ
    K_x = FZ * (CFX9 + CFX10 * df_z) * np.exp(CFX11 * df_z) * self.scaling_coeffs[4]
    self.B_x = K_x / (self.C_x * self.D_x)

    S_Hx = (CFX12 + CFX13 * df_z) * self.scaling_coeffs[5]
    S_Vx = FZ * (CFX14 + CFX15 * df_z) * self.scaling_coeffs[6] * self.scaling_coeffs[2]
    SR_x = SR + S_Hx

    self.E_x = (CFX5 + CFX6 * df_z + CFX7 * df_z**2) * (1 - CFX8 * np.sign(SR_x)) * self.scaling_coeffs[3]
    
    self.F_X0 = self.D_x * np.sin(self.C_x * np.arctan(self.B_x * SR_x - self.E_x * (self.B_x * SR_x - np.arctan(self.B_x * SR_x)))) + S_Vx
    F_X = self.F_X0
    
    return self.E_x > 1 if sign_condition else F_X