"""Independent analytic oracle for the documented helium DIS fit-B checks."""
import math


def r1990_fit_b(x: float, q2: float) -> float:
    """Whitlow R1990 fit B, with Q2 expressed in GeV^2."""
    if x <= 0.0 or q2 <= 0.04:
        raise ValueError("R1990 requires x > 0 and Q2 > 0.04 GeV2")
    scale = 0.125**2
    theta = 1.0 + 12.0 * q2 / (q2 + 1.0) * scale / (scale + x * x)
    return 0.0635 / math.log(q2 / 0.04) * theta + 0.5747 / q2 - 0.3534 / (q2 * q2 + 0.09)
