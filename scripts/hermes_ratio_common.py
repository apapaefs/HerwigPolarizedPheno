"""Display-only MC/data ratios, retaining signed and poorly determined data.

Experimental uncertainties are shown separately about unity. They are not
Gaussian confidence intervals for a quotient with an uncertain denominator.
"""
from __future__ import annotations
import math

COLOR = '#d55e00'
FOCUS = (-1., 3.)


def ratio_record(data, stat, total, mc, mcstat, supported=True):
    """Compute a fixed-central-denominator display, without physics corrections."""
    if not isinstance(supported, bool):
        raise ValueError('Ratio support must be boolean')
    if not supported and all(v is None for v in (data, stat, total, mc, mcstat)):
        return {'ratio': None, 'ratio_mc_stat': None, 'data_relative_stat': None,
                'data_relative_total': None, 'denominator_zero_compatible': False,
                'ratio_missing_reason': 'No physical projection support'}
    for name, value in (('data', data), ('stat', stat), ('total', total)):
        if isinstance(value, bool) or value is None or not math.isfinite(float(value)):
            raise ValueError(f'Invalid ratio {name}')
    data, stat, total = map(float, (data, stat, total))
    if stat < 0 or total < 0 or total < stat - 1e-12:
        raise ValueError('Ratio errors require total >= statistical >= 0')
    if (mc is None) != (mcstat is None):
        raise ValueError('MC value/error must be present or missing together')
    if mc is not None:
        if any(isinstance(v, bool) or not math.isfinite(float(v)) for v in (mc, mcstat)) or float(mcstat) < 0:
            raise ValueError('Invalid MC ratio input')
        if not supported:
            raise ValueError('Unsupported cells cannot have an MC prediction')
    reason = ('No MC support' if not supported else 'Missing MC prediction' if mc is None
              else 'Data central value is zero' if data == 0 else None)
    denominator = abs(data)
    return {'ratio': float(mc)/data if reason is None else None,
            'ratio_mc_stat': float(mcstat)/denominator if reason is None else None,
            'data_relative_stat': stat/denominator if data else None,
            'data_relative_total': total/denominator if data else None,
            'denominator_zero_compatible': abs(data) <= total,
            'ratio_missing_reason': reason}


def ratio_limits(rows, focus=True):
    if focus:
        return list(FOCUS)
    extents = [0., 1., 2.]
    for row in rows:
        if row.get('ratio') is None:
            continue
        extents.extend([row['ratio']-row['ratio_mc_stat'], row['ratio']+row['ratio_mc_stat'],
                        1-row['data_relative_total'], 1+row['data_relative_total']])
    lo, hi = min(extents), max(extents)
    pad = .09*max(hi-lo, 1.)
    return [lo-pad, hi+pad]


def draw_ratio(axis, rows, ylim=FOCUS):
    """Draw separate MC-stat and data-relative displays; mark every truncation.

    The caller owns x axes/titles. Orange triangles mark clipped MC centers
    (labelled numerically) or MC error endpoints. Gray triangles mark clipped
    experimental reference bands. Hollow MC markers flag zero-compatible data.
    """
    from matplotlib.patches import Rectangle
    lo, hi = map(float, ylim)
    if not lo < hi:
        raise ValueError('Invalid ratio axis limits')
    axis.set_ylim(lo, hi)
    axis.axhline(1, color='0.25', linewidth=.65, linestyle='--', zorder=3)
    axis.axhline(0, color='0.75', linewidth=.45, zorder=0)
    axis.grid(axis='y', color='0.92', linewidth=.4, zorder=0)
    span = hi-lo
    for row in rows:
        left, right, marker = row['low'], row['high'], row['marker']
        if not (0 < left < right and left <= marker <= right):
            raise ValueError('Invalid ratio bin geometry')
        value = row.get('ratio')
        if value is None:
            if row.get('supported', True):
                label = 'data = 0' if row['data'] == 0 else 'MC missing'
                axis.text(marker, lo+.10*span, label, ha='center', va='bottom',
                          fontsize=6, color='0.4', rotation=90)
            continue
        stat, total = row['data_relative_stat'], row['data_relative_total']
        for error, color in ((total, '0.88'), (stat, '0.73')):
            bottom, top = max(lo, 1-error), min(hi, 1+error)
            if top > bottom:
                axis.add_patch(Rectangle((left, bottom), right-left, top-bottom,
                                         facecolor=color, edgecolor='none', zorder=1))
        for edge, direction, glyph in ((1-total, -1, 'v'), (1+total, 1, '^')):
            if edge < lo or edge > hi:
                y = lo+.012*span if direction < 0 else hi-.012*span
                axis.plot(marker, y, marker=glyph, color='0.48', markersize=3,
                          linestyle='none', zorder=5)
        error = row['ratio_mc_stat']
        bottom, top = max(lo, value-error), min(hi, value+error)
        if top > bottom:
            axis.add_patch(Rectangle((left, bottom), right-left, top-bottom,
                                     facecolor=COLOR, alpha=.20, edgecolor='none', zorder=2))
        if lo <= value <= hi:
            axis.hlines(value, left, right, color=COLOR, linewidth=1.0, zorder=4)
            axis.plot(marker, value, marker='o', markersize=3.1,
                      markerfacecolor='white' if row['denominator_zero_compatible'] else COLOR,
                      markeredgecolor=COLOR, markeredgewidth=.8, linestyle='none', zorder=6)
            for edge, direction, glyph in ((value-error,-1,'v'),(value+error,1,'^')):
                if edge < lo or edge > hi:
                    edge_marker = math.sqrt(left*marker) if direction < 0 else math.sqrt(marker*right)
                    axis.plot(edge_marker, lo+.025*span if direction < 0 else hi-.025*span,
                              marker=glyph, color=COLOR, markersize=3.2, linestyle='none', zorder=6)
        else:
            # A wide interval may also overflow the boundary opposite its center.
            if value > hi and value-error < lo:
                axis.plot(marker,lo+.025*span,marker='v',color=COLOR,markersize=3.2,linestyle='none',zorder=6)
            if value < lo and value+error > hi:
                axis.plot(marker,hi-.025*span,marker='^',color=COLOR,markersize=3.2,linestyle='none',zorder=6)
            upper = value > hi
            y = hi-.055*span if upper else lo+.055*span
            axis.plot(marker, y, marker='^' if upper else 'v', color=COLOR,
                      markerfacecolor='white' if row['denominator_zero_compatible'] else COLOR,
                      markeredgewidth=.8, markersize=4, linestyle='none', zorder=6)
            axis.annotate(f'{value:.3g}', (marker,y), xytext=(0,-6 if upper else 6),
                          textcoords='offset points', ha='center',
                          va='top' if upper else 'bottom', fontsize=6.5, color=COLOR, zorder=7)
    return axis
