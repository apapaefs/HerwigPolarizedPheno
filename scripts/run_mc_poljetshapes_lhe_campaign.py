#!/usr/bin/env python3
"""Pinned entry point for the combined hard-spin-transfer campaign."""

from __future__ import annotations

import sys
from typing import Sequence

import run_phenomenology_campaign as campaign


MEASUREMENT = "MC_POLJETSHAPES_LHE"


def pinned_arguments(argv: Sequence[str]) -> list[str]:
    arguments = list(argv)
    if not arguments:
        return arguments
    if "--measurement" in arguments:
        raise campaign.CampaignError(
            "The MC_POLJETSHAPES_LHE runner pins its measurement; do not pass "
            "--measurement"
        )
    if arguments[0] == "list":
        raise campaign.CampaignError(
            "The dedicated runner accepts fetch-data, prepare, campaign, "
            "postprocess, plot, or full"
        )
    return [arguments[0], "--measurement", MEASUREMENT, *arguments[1:]]


def main(argv: Sequence[str] | None = None) -> int:
    try:
        return campaign.main(pinned_arguments(sys.argv[1:] if argv is None else argv))
    except campaign.CampaignError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 1


if __name__ == "__main__":
    raise SystemExit(main())
