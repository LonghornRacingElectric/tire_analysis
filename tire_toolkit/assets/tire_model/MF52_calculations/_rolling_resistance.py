def _combined_rolling(self, data: list[float]) -> float:
    FZ, SA, SR, IA = data

    # TTC doesn't give FY data, so we'll make an estimate rather than using the calculations below
    # [CCMY1, CCMY2, CCMY3, CCMY4, VREF] = self.combined_rolling_coeffs
    # FX = self._combined_long([FZ, SA, SR, IA])
    # MY = self.R_nom * FZ * (CCMY1 + CCMY2 * FX / self.FZ_nom + CCMY3 * np.abs(V_x / VREF) + CCMY4 * (V_x / VREF)**4)

    [CFX1, CFX2, CFX3, CFX4, CFX5, CFX6, CFX7, CFX8, CFX9, CFX10, CFX11, CFX12, CFX13, CFX14, CFX15] = self.pure_long_coeffs

    df_z = (FZ - self.FZ_nom * self.scaling_coeffs[0]) / (self.FZ_nom * self.scaling_coeffs[0])

    K_x = FZ * (CFX9 + CFX10 * df_z) * np.exp(CFX11 * df_z) * self.scaling_coeffs[4]

    S_Hx = (CFX12 + CFX13 * df_z) * self.scaling_coeffs[5]
    S_Vx = FZ * (CFX14 + CFX15 * df_z) * self.scaling_coeffs[6] * self.scaling_coeffs[2]

    MY = self.R_nom * (S_Vx + K_x * S_Hx)

    return MY