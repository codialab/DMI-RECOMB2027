"""Fail-closed boundary-source checks for Stage-11 GEMs.

COBRA ``model.exchanges`` is not the complete set of boundary reactions.
Demand/sink/custom one-metabolite boundary reactions can also provide mass
when their stoichiometric orientation and bounds permit it. Medium QC must
therefore reason over ``model.boundary`` and reaction stoichiometry.
"""
from __future__ import annotations

from collections.abc import Iterable
from typing import Any, Mapping


def _boundary_coefficient(reaction) -> float:
    """Return the sole metabolite coefficient of a COBRA boundary reaction."""
    if not getattr(reaction, "boundary", False):
        raise ValueError(f"Reaction is not a COBRA boundary reaction: {reaction.id}")
    coefficients = list(reaction.metabolites.values())
    if len(coefficients) != 1:
        raise RuntimeError(
            f"Boundary reaction {reaction.id} has {len(coefficients)} metabolites; "
            "cannot determine source direction safely"
        )
    coefficient = float(coefficients[0])
    if coefficient == 0.0:
        raise RuntimeError(
            f"Boundary reaction {reaction.id} has zero stoichiometric coefficient"
        )
    return coefficient


def boundary_source_open(reaction, tolerance: float = 0.0) -> bool:
    """Return whether a boundary reaction can create its boundary metabolite."""
    tolerance = float(tolerance)
    if tolerance < 0.0:
        raise ValueError("tolerance must be non-negative")
    coefficient = _boundary_coefficient(reaction)
    if coefficient < 0.0:
        return float(reaction.lower_bound) < -tolerance
    return float(reaction.upper_bound) > tolerance


def boundary_source_direction(reaction) -> str:
    """Return the flux sign which supplies the boundary metabolite."""
    return "negative" if _boundary_coefficient(reaction) < 0.0 else "positive"


def boundary_source_capacity(reaction) -> float:
    """Return currently available capacity in the source-producing sign."""
    if boundary_source_direction(reaction) == "negative":
        return max(0.0, -float(reaction.lower_bound))
    return max(0.0, float(reaction.upper_bound))


def boundary_source_routes(model, tolerance: float = 0.0) -> dict[str, dict[str, Any]]:
    """Inventory open source routes with their flux sign and capacity."""
    routes: dict[str, dict[str, Any]] = {}
    for reaction in model.boundary:
        if boundary_source_open(reaction, tolerance=tolerance):
            routes[reaction.id] = {
                "source_direction": boundary_source_direction(reaction),
                "source_capacity": boundary_source_capacity(reaction),
            }
    return dict(sorted(routes.items()))


def active_boundary_source_fluxes(model, fluxes: Mapping[str, float], tolerance: float = 1.0e-9) -> dict[str, float]:
    """Report active nutrient/source fluxes using the same stoichiometric sign audit."""
    active: dict[str, float] = {}
    for reaction in model.boundary:
        value = float(fluxes.get(reaction.id, 0.0))
        if (boundary_source_direction(reaction) == "negative" and value < -tolerance) or (
            boundary_source_direction(reaction) == "positive" and value > tolerance
        ):
            active[reaction.id] = value
    return dict(sorted(active.items()))


def close_boundary_source_direction(reaction) -> None:
    """Close only source flux, preserving the reaction's legitimate drain sign."""
    if boundary_source_direction(reaction) == "negative":
        reaction.lower_bound = max(0.0, float(reaction.lower_bound))
    else:
        reaction.upper_bound = min(0.0, float(reaction.upper_bound))


def close_boundary_source_direction_for_negative_control(reaction) -> None:
    """Close source flux on a disposable negative-control model.

    RIPTiDe ``set_bounds=True`` can tighten a boundary reaction to a strictly
    source-only interval, for example ``(-0.05, -0.002)`` for a conventional
    exchange.  Closing only the source sign must then make zero flux available;
    otherwise setting the lower bound first would create an invalid interval.

    Reversible/drain-capable intervals retain their drain sign.  Strictly
    source-only intervals collapse to ``(0, 0)``.  This helper is intentionally
    separate from :func:`close_boundary_source_direction`, whose semantics are
    part of frozen-medium construction.
    """
    lower = float(reaction.lower_bound)
    upper = float(reaction.upper_bound)
    if boundary_source_direction(reaction) == "negative":
        reaction.bounds = (max(0.0, lower), max(0.0, upper))
    else:
        reaction.bounds = (min(0.0, lower), min(0.0, upper))


def assert_resolved_source_routes(
    model,
    expected_routes: Mapping[str, Mapping[str, Any]],
    *,
    context: str,
    tolerance: float = 1.0e-9,
) -> None:
    """Require exact resolved-medium route IDs, signs, and capacities."""
    expected = {
        str(rid): {
            "source_direction": str(route["source_direction"]),
            "source_capacity": float(route["source_capacity"]),
        }
        for rid, route in expected_routes.items()
    }
    actual = boundary_source_routes(model, tolerance=tolerance)
    if set(actual) != set(expected):
        raise RuntimeError(
            f"{context} source-route IDs differ: expected={sorted(expected)}, actual={sorted(actual)}"
        )
    for rid, route in expected.items():
        observed = actual[rid]
        if observed["source_direction"] != route["source_direction"]:
            raise RuntimeError(
                f"{context} source direction differs for {rid}: "
                f"expected={route['source_direction']}, actual={observed['source_direction']}"
            )
        if abs(observed["source_capacity"] - route["source_capacity"]) > tolerance:
            raise RuntimeError(
                f"{context} source capacity differs for {rid}: "
                f"expected={route['source_capacity']}, actual={observed['source_capacity']}"
            )


def assert_source_routes_within_resolved_medium(
    model,
    expected_routes: Mapping[str, Mapping[str, Any]],
    *,
    context: str,
    tolerance: float = 1.0e-9,
) -> None:
    """Allow contextualization to prune/tighten medium routes, never widen them.

    This is intentionally distinct from ``assert_resolved_source_routes``:
    the parent must exactly reproduce the frozen medium, whereas an extracted
    expert can omit an unused allowed exchange.
    """
    expected = {
        str(rid): {
            "source_direction": str(route["source_direction"]),
            "source_capacity": float(route["source_capacity"]),
        }
        for rid, route in expected_routes.items()
    }
    actual = boundary_source_routes(model, tolerance=tolerance)
    unauthorized = sorted(set(actual) - set(expected))
    if unauthorized:
        raise RuntimeError(f"{context} contains unauthorized source routes: {unauthorized}")
    for rid, observed in actual.items():
        route = expected[rid]
        if observed["source_direction"] != route["source_direction"]:
            raise RuntimeError(
                f"{context} source direction differs for {rid}: "
                f"expected={route['source_direction']}, actual={observed['source_direction']}"
            )
        if observed["source_capacity"] > route["source_capacity"] + tolerance:
            raise RuntimeError(
                f"{context} source capacity widened for {rid}: "
                f"configured={route['source_capacity']}, actual={observed['source_capacity']}"
            )


def open_boundary_source_ids(model, tolerance: float = 0.0) -> list[str]:
    """Return sorted IDs of all currently source-capable boundary reactions."""
    return sorted(
        reaction.id
        for reaction in model.boundary
        if boundary_source_open(reaction, tolerance=tolerance)
    )


def unexpected_boundary_sources(
    model,
    allowed_source_ids: Iterable[str],
    tolerance: float = 0.0,
) -> list[str]:
    """Return source-capable boundary IDs outside an explicit allow-list."""
    allowed = {str(value) for value in allowed_source_ids}
    return sorted(
        reaction_id
        for reaction_id in open_boundary_source_ids(model, tolerance=tolerance)
        if reaction_id not in allowed
    )


def assert_no_unexpected_boundary_sources(
    model,
    allowed_source_ids: Iterable[str],
    *,
    context: str,
    tolerance: float = 0.0,
) -> None:
    """Fail closed when a model has an unauthorized boundary source route."""
    unexpected = unexpected_boundary_sources(
        model, allowed_source_ids, tolerance=tolerance
    )
    if unexpected:
        raise RuntimeError(
            f"{context} contains unauthorized boundary source routes: {unexpected}"
        )
