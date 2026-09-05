"""Nominal fixed-beam support for the published COMPASS multiplicity cells."""
import math


def unsupported_cells(snapshot, dataset, species, target_weights):
    if "nu_window" not in snapshot.get("selection", {}):
        return []
    energy = float(snapshot.get("beam_energy_gev", snapshot.get("beam", {}).get("energy_gev")))
    mass = .493677 if species.startswith("k") else .13957039
    selection = snapshot["selection"]
    pmin, pmax = selection["hadron_momentum_gev"]
    q2min = float(selection["q2_min_gev2"])
    w2min = float(selection["w_min_gev"])**2
    unsupported = []
    for index, point in enumerate(dataset["points"]):
        supported = False
        for component, weight in target_weights.items():
            if weight == 0.:
                continue
            target_mass = {"P": .9382720813, "N": .9395654133}[component]
            scale = 2.*target_mass*energy
            # Existence of x in this rectangle with Q2>q2min and W2>w2min.
            lower = max(
                point["y_low"], selection["y"][0],
                q2min/(scale*point["x_high"]),
                (w2min-target_mass**2)/(scale*(1.-point["x_low"])),
                (q2min+w2min-target_mass**2)/scale,
                math.hypot(pmin, mass)/(energy*point["z_low"]),
            )
            upper = min(point["y_high"], selection["y"][1],
                        math.hypot(pmax, mass)/(energy*point["z_high"]))
            supported |= lower < upper
        if not supported:
            unsupported.append(index)
    return unsupported
