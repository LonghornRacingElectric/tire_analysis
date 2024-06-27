def _combined_aligning(self, data: list[float]) -> float:
    FZ, SA, SR, IA = data

    CCMZ1, CCMZ2, CCMZ3, CCMZ4 = self.combined_aligning_coeffs
    CMZ1, CMZ2, CMZ3, CMZ4, CMZ5, CMZ6, CMZ7, CMZ8, CMZ9, CMZ10, CMZ11, CMZ12, CMZ13, CMZ14, CMZ15, CMZ16, CMZ17, CMZ18, CMZ19, CMZ20, CMZ21, CMZ22, CMZ23, CMZ24, CMZ25 = self.pure_aligning_coeffs
    CCFY1, CCFY2, CCFY3, CCFY4, CCFY5, CCFY6, CCFY7, CCFY8, CCFY9, CCFY10, CCFY11, CCFY12, CCFY13, CCFY14 = self.combined_lat_coeffs
    CFY1, CFY2, CFY3, CFY4, CFY5, CFY6, CFY7, CFY8, CFY9, CFY10, CFY11, CFY12, CFY13, CFY14, CFY15, CFY16, CFY17, CFY18 = self.pure_lat_coeffs
    CFX1, CFX2, CFX3, CFX4, CFX5, CFX6, CFX7, CFX8, CFX9, CFX10, CFX11, CFX12, CFX13, CFX14, CFX15 = self.pure_long_coeffs

    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])
    IA_z = IA * self.scaling_coeffs[17]

    # Pure slip dependencies
    # Lateral
    IA_y = IA * self.scaling_coeffs[14]
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * self.scaling_coeffs[9]
    S_Hy = (CFY12 + CFY13 * df_z) * self.scaling_coeffs[12] + CFY14 * IA_y
    S_Vy = FZ * ((CFY15 + CFY16 * df_z) * self.scaling_coeffs[13] + (CFY17 + CFY18 * df_z) * IA_y) * self.scaling_coeffs[9]

    K_y = CFY9 * self.FZ_nom * np.sin(2 * np.arctan(FZ / (CFY10 * self.FZ_nom * self.scaling_coeffs[0]))) * \
        (1 - CFY11 * abs(IA_y)) * self.scaling_coeffs[0] * self.scaling_coeffs[11]
    C_y = CFY1 * self.scaling_coeffs[8]
    D_y = mu_y * FZ
    B_y = K_y / (C_y * D_y)

    # Longitudinal
    K_x = FZ * (CFX9 + CFX10 * df_z) * np.exp(CFX11 * df_z) * self.scaling_coeffs[4]
    
    # Aligning
    D_t = FZ * (CMZ9 + CMZ10 * df_z) * (1 + CMZ11 * IA_z + CMZ12 * IA_z**2) * (self.R_nom / self.FZ_nom) * self.scaling_coeffs[15]
    C_t = CMZ8
    B_t = (CMZ1 + CMZ2 * df_z + CMZ3 * df_z**2) * (1 + CMZ4 * IA_z + CMZ5 * abs(IA_z)) * self.scaling_coeffs[11] / self.scaling_coeffs[9]

    E_t = (CMZ17 + CMZ18 * df_z + CMZ19 * df_z**2) * (1 + (CMZ20 + CMZ21 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))
    
    D_VySR = mu_y * FZ * (CCFY9 + CCFY10 * df_z + CCFY11 * IA) * np.cos(np.arctan(CCFY12 * SA))
    S_VySR = D_VySR * np.sin(CCFY13 * np.arctan(CCFY14 * SR)) * self.scaling_coeffs[23]

    S_Ht = CMZ22 + CMZ23 * df_z + (CMZ24 + CMZ25 * df_z) * IA_z
    SA_t = SA + S_Ht

    S_Hf = S_Hy + S_Vy / K_y
    SA_r = SA + S_Hf

    # Adjusted SA values
    SA_t_eq = np.arctan(np.sqrt((np.tan(SA_t))**2 + (K_x / K_y)**2 * SR**2)) * np.sign(SA_t)
    SA_r_eq = np.arctan(np.sqrt((np.tan(SA_r))**2 + (K_x / K_y)**2 * SR**2)) * np.sign(SA_r)

    # Pneumatic trail
    t_adj = D_t * np.cos(C_t * np.arctan(B_t * SA_t_eq - E_t * (B_t * SA_t_eq - np.arctan(B_t * SA_t_eq)))) * np.cos(SA)

    FX = self._combined_long([FZ, SA, SR, IA])
    FY = self._combined_lat([FZ, SA, SR, IA])

    F_y_IA_adj = FY - S_VySR

    D_r = FZ * ((CMZ13 + CMZ14 * df_z) * self.scaling_coeffs[16] + (CMZ15 + CMZ16 * df_z) * IA_z) * self.R_nom * self.scaling_coeffs[9]
    B_r = CMZ6 * self.scaling_coeffs[11] / self.scaling_coeffs[9] + CMZ7 * B_y * C_y

    M_zr = D_r * np.cos(np.arctan(B_r * SA_r_eq)) * np.cos(SA)

    s = (CCMZ1 + CCMZ2 * (FY / self.FZ_nom) + (CCMZ3 + CCMZ4 * df_z) * IA) * self.R_nom * self.scaling_coeffs[24]

    MZ_adj = -t_adj * F_y_IA_adj + M_zr + s * FX

    return MZ_adj

def _pure_aligning(self, data: list[float]) -> float:
    FZ, SA, IA = data

    # Lateral Dependencies
    [CFY1, CFY2, CFY3, CFY4, CFY5, CFY6, CFY7, CFY8, CFY9, CFY10, \
        CFY11, CFY12, CFY13, CFY14, CFY15, CFY16, CFY17, CFY18] = self.pure_lat_coeffs
    
    IA_y = IA * self.scaling_coeffs[14]
    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])
    mu_y = (CFY2 + CFY3 * df_z) * (1 - CFY4 * IA_y**2) * self.scaling_coeffs[9]

    C_y = CFY1 * self.scaling_coeffs[8]
    D_y = mu_y * FZ
    K_y = CFY9 * self.FZ_nom * np.sin(2 * np.arctan(FZ / (CFY10 * self.FZ_nom * self.scaling_coeffs[0]))) * \
        (1 - CFY11 * abs(IA_y)) * self.scaling_coeffs[0] * self.scaling_coeffs[11]
    B_y = K_y / (C_y * D_y)
    F_y0 = self._pure_lat([FZ, SA, IA])

    B_y = K_y / (C_y * D_y)
    S_Hy = (CFY12 + CFY13 * df_z) * self.scaling_coeffs[12] + CFY14 * IA_y
    S_Vy = FZ * ((CFY15 + CFY16 * df_z) * self.scaling_coeffs[13] + (CFY17 + CFY18 * df_z) * IA_y) * self.scaling_coeffs[9]

    # Pure Aligning Moment
    [CMZ1, CMZ2, CMZ3, CMZ4, CMZ5, CMZ6, CMZ7, CMZ8, CMZ9, CMZ10, CMZ11, CMZ12, \
        CMZ13, CMZ14, CMZ15, CMZ16, CMZ17, CMZ18, CMZ19, CMZ20, CMZ21, CMZ22, CMZ23, CMZ24, CMZ25] = self.pure_aligning_coeffs
    
    IA_z = IA * self.scaling_coeffs[17]
    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])

    S_Ht = CMZ22 + CMZ23 * df_z + (CMZ24 + CMZ25 * df_z) * IA_z
    SA_t = SA + S_Ht

    D_r = FZ * ((CMZ13 + CMZ14 * df_z) * self.scaling_coeffs[16] + (CMZ15 + CMZ16 * df_z) * IA_z) * self.R_nom * self.scaling_coeffs[9]
    B_r = CMZ6 * self.scaling_coeffs[11] / self.scaling_coeffs[9] + CMZ7 * B_y * C_y

    D_t = FZ * (CMZ9 + CMZ10 * df_z) * (1 + CMZ11 * IA_z + CMZ12 * IA_z**2) * (self.R_nom / self.FZ_nom) * self.scaling_coeffs[15]
    C_t = CMZ8
    B_t = (CMZ1 + CMZ2 * df_z + CMZ3 * df_z**2) * (1 + CMZ4 * IA_z + CMZ5 * abs(IA_z)) * self.scaling_coeffs[11] / self.scaling_coeffs[9]

    E_t = (CMZ17 + CMZ18 * df_z + CMZ19 * df_z**2) * (1 + (CMZ20 + CMZ21 * IA_z) * (2 / np.pi) * np.arctan(B_t * C_t * SA_t))

    S_Hf = S_Hy + S_Vy / K_y
    
    # Residual Torque

    # Found typo in the delft paper. Hr should be Hf
    SA_r = SA + S_Hf
    M_zr = D_r * np.cos(np.arctan(B_r * SA_r)) * np.cos(SA)

    # Pneumatic trail
    t = D_t * np.cos(C_t * np.arctan(B_t * SA_t - E_t * (B_t * SA_t - np.arctan(B_t * SA_t)))) * np.cos(SA)

    M_Z0 = -t * F_y0 + M_zr
    M_Z = M_Z0

    return M_Z