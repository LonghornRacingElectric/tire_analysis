def _combined_overturning(self, data: list[float]) -> float:
    FZ, SA, SR, IA = data

    [CCMX1, CCMX2, CCMX3] = self.overturning_coeffs

    FY = self._combined_lat([FZ, SA, SR, IA])

    MX = self.R_nom * FZ * (CCMX1 * self.scaling_coeffs[19] + (-CCMX2 * IA + CCMX3 * FY / self.FZ_nom) * self.scaling_coeffs[18])
    
    return MX