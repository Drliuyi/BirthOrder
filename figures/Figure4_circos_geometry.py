"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import hashlib
import json
import math
import os
import platform
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path
import matplotlib
import matplotlib.pyplot as plt
from matplotlib import font_manager
from matplotlib.font_manager import FontProperties
from matplotlib.textpath import TextToPath
from matplotlib.patches import Circle, Rectangle, Wedge, FancyBboxPatch, Polygon
from matplotlib.path import Path as MPath
from matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb
import numpy as np
import pandas as pd
DISEASES = ['250', '250.2', '411', '411.3', '411.4', '351']
SHORT = ['DM', 'T2D', 'IHD', 'Angina', 'Cor. ather.', '351']
FULL = ['DM', 'T2D-related PheCode', 'IHD', 'Angina', 'Coronary atherosclerosis-related PheCode', 'Other peripheral nervous-system disorders']
STATE_COLOURS = {'SUPPORTED_NEGATIVE': '#2C6E96', 'SUPPORTED_POSITIVE': '#CE6C63', 'ESTIMABLE_NOT_SUPPORTED': '#DEE7EB', 'LIMITED_ESTIMATE': '#EAE5DC', 'OTHER_REVIEW': '#D6CFD9'}
DOMAINS = []
NAVY = '#243E4C'
INK = '#34434A'
EFFECT_ANCHORS = ['#1F5D88', '#9CBFD6', '#FFFFFF', '#E7B1A8', '#AF4540']
NO_ESTIMATE_COLOUR = '#E3E5E7'

def apply_effect_palette(d):
    """Colour derives only from finite original product; never from CI/status/P."""
    d['finite_ab'] = np.isfinite(d.ab)
    assert int(d.finite_ab.sum()) == 532, 'Finite count differs from frozen 532: stop, do not fabricate.'
    assert not (d.supported & ~d.finite_ab).any(), 'Original supported flag has nonfinite product.'
    c = float(np.quantile(np.abs(d.loc[d.finite_ab, 'ab']), 0.98, method='linear'))
    assert math.isfinite(c) and c > 0
    d['color_ab'] = np.where(d.finite_ab, np.clip(d.ab, -c, c), np.nan)
    d['clipped_for_color'] = d.finite_ab & d.ab.abs().gt(c)
    d['fill'] = [to_hex(EFFECT_CMAP((value + c) / (2 * c))).upper() if finite else NO_ESTIMATE_COLOUR for value, finite in zip(d.color_ab, d.finite_ab)]
    rgb = np.array([to_rgb(colour) for colour in d.fill])
    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    luminance = linear @ np.array([0.2126, 0.7152, 0.0722])
    d['relative_luminance'] = luminance
    d['contrast_black'] = (luminance + 0.05) / 0.05
    d['contrast_white'] = 1.05 / (luminance + 0.05)
    d['star_color'] = ['#FFFFFF' if support and white > black else '#000000' if support else '' for white, black, support in zip(d.contrast_white, d.contrast_black, d.supported)]
    d.loc[~d.supported, 'star_color'] = ''
    d['point_estimate_warning'] = ''
    review = d.finite_ab & d['估计状态'].isin(['BOOT_UNSTABLE', 'FAILURE_REVIEW', 'MODEL_LIMITATION'])
    d.loc[review, 'point_estimate_warning'] = 'Finite exploratory point estimate; diagnostic/bootstrap limitation retained; colour is not reliable inference.'
    d.loc[~d.finite_ab, 'point_estimate_warning'] = 'No finite product; not filled with zero.'
    scale = dict(quantile=0.98, quantile_method='NumPy linear; R type 7', finite_source_n=int(d.finite_ab.sum()), center=0, c=c, vmin=-c, vmax=c, palette_anchors=EFFECT_ANCHORS, palette_steps=1025, missing_colour=NO_ESTIMATE_COLOUR, raw_min=float(d.loc[d.finite_ab, 'ab'].min()), raw_max=float(d.loc[d.finite_ab, 'ab'].max()), clipped_negative=int((d.finite_ab & d.ab.lt(-c)).sum()), clipped_positive=int((d.finite_ab & d.ab.gt(c)).sum()), original_ab_values_changed=False, applies_to='Both main80 and supplement140', label='a × b (Cox log-hazard coefficient product)', star_rule='Original-batch evidence flag only; supported + nonfinite is a stop condition', contrast_rule='Choose black/white with greater contrast against actual sRGB hex fill; linearized relative luminance')
    return (d, scale)

def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()

def save_csv(df, path):
    if path.exists():
        raise FileExistsError(path)
    df.to_csv(path, index=False, encoding='utf-8-sig')

def load_data(source, mapping, results, display_mapping):
    raw = pd.read_csv(source, dtype=str, keep_default_na=False)
    assert raw.shape == (1080, 28), raw.shape
    assert raw['候选ID'].is_unique
    assert (raw['类型'] == '分类关联_非连续乘积').sum() == 240
    d = raw.loc[raw['类型'].eq('连续路径')].copy()
    d['source_row'] = d.index + 1
    assert len(d) == 840 and d['原指标ID'].nunique() == 140
    assert not d.duplicated(['原指标ID', '疾病编码']).any()
    assert set(d['疾病编码']) == set(DISEASES)
    assert d.groupby('原指标ID')['疾病编码'].apply(lambda x: set(x) == set(DISEASES)).all()
    for source_col, col in [('a乘b', 'ab'), ('乘积CI下限', 'ci_low'), ('乘积CI上限', 'ci_high'), ('原批次联合BH_q', 'frozen_q')]:
        d[col] = pd.to_numeric(d[source_col], errors='coerce')
    valid = np.isfinite(d[['ab', 'ci_low', 'ci_high', 'frozen_q']]).all(axis=1)
    same_direction = d.ab.gt(0) & d.ci_low.gt(0) | d.ab.lt(0) & d.ci_high.lt(0)
    flag = d['原批次证据'].eq('EXPLORATORY_PATH_SUPPORTED')
    supported = flag
    assert not (supported & ~np.isfinite(d.ab)).any(), 'Frozen supported position has no finite product.'
    assert (flag == flag & d['估计状态'].eq('CI_AVAILABLE') & d.frozen_q.lt(0.05) & valid & same_direction).all(), 'Frozen source audit conflict.'
    d['supported'] = supported
    d['state'] = 'OTHER_REVIEW'
    d.loc[d['估计状态'].isin(['MODEL_LIMITATION', 'BOOT_UNSTABLE', 'FAILURE_REVIEW']), 'state'] = 'LIMITED_ESTIMATE'
    d.loc[d['估计状态'].eq('CI_AVAILABLE') & ~supported, 'state'] = 'ESTIMABLE_NOT_SUPPORTED'
    d.loc[supported & d.ab.lt(0), 'state'] = 'SUPPORTED_NEGATIVE'
    d.loc[supported & d.ab.gt(0), 'state'] = 'SUPPORTED_POSITIVE'
    assert (d.state == 'OTHER_REVIEW').sum() == 0, 'Unmapped state: explicit review required.'
    assert supported.sum() == 291
    mids = set(d.loc[supported, '原指标ID'])
    assert len(mids) == 80
    d['in_main80'] = d['原指标ID'].isin(mids)
    assert d.in_main80.sum() == 480
    m = pd.read_csv(mapping, dtype=str, keep_default_na=False)
    assert m.marker_id.is_unique and len(m) == 140
    assert set(m.marker_id) == set(d['原指标ID'])
    source_domains = d.groupby('原指标ID')['领域'].agg(lambda x: list(set(x)))
    assert source_domains.apply(len).eq(1).all()
    assert all((source_domains[row.marker_id][0] == row.source_domain for row in m.itertuples()))
    assert set(m.source_domain) == {x[0] for x in DOMAINS}
    dm = pd.read_csv(display_mapping, dtype=str, keep_default_na=False)
    assert len(dm) == 140 and dm.marker_id.is_unique
    assert set(dm.marker_id) == set(m.marker_id)
    assert set(dm.display_domain_id) == {row[0] for row in DOMAINS}
    assert m.set_index('marker_id').source_domain.eq(dm.set_index('marker_id').source_domain).all()
    changed = dm.loc[dm.source_domain.ne(dm.display_domain_id), 'marker_id'].tolist()
    assert changed == ['commute_miles'], changed
    m['source_marker_order'] = pd.to_numeric(m.marker_order)
    m['legacy_domain_label'] = m.domain
    m['source_domain_order'] = pd.to_numeric(m.domain_order)
    add = dm[['marker_id', 'source_domain', 'display_domain_id', 'display_label', 'colour', 'display_order', 'change_reason']].rename(columns={'display_label': 'display_domain', 'colour': 'display_domain_colour', 'display_order': 'display_domain_order'})
    m = m.merge(add, on=['marker_id', 'source_domain'], how='left', validate='one_to_one')
    assert m.display_domain_id.notna().all()
    m['display_domain_order'] = pd.to_numeric(m.display_domain_order)
    m = m.sort_values(['display_domain_order', 'source_marker_order']).reset_index(drop=True)
    m['marker_order'] = np.arange(1, len(m) + 1)
    m['domain'] = m.display_domain
    m['domain_order'] = m.display_domain_order
    m['within_domain_order'] = m.groupby('display_domain_id', sort=False).cumcount() + 1
    d['marker_id'] = d['原指标ID']
    d['source_domain'] = d['领域']
    d['display_domain_id'] = d.marker_id.map(m.set_index('marker_id').display_domain_id)
    d['display_domain'] = d.marker_id.map(m.set_index('marker_id').display_domain)
    d['display_domain_colour'] = d.marker_id.map(m.set_index('marker_id').display_domain_colour)
    d['display_change_reason'] = d.marker_id.map(m.set_index('marker_id').change_reason)
    d['domain_id'] = d.display_domain_id
    order = dict(zip(m.marker_id, m.marker_order))
    d['marker_order'] = d.marker_id.map(order)
    d['ring_order'] = d['疾病编码'].map({v: k + 1 for k, v in enumerate(DISEASES)})
    d, scale = apply_effect_palette(d)
    d = d.sort_values(['marker_order', 'ring_order'])
    counts = []
    for code, name, exp in zip(DISEASES, FULL, EXPECTED):
        z = d[d['疾病编码'].eq(code)]
        neg = int(z.state.eq('SUPPORTED_NEGATIVE').sum())
        pos = int(z.state.eq('SUPPORTED_POSITIVE').sum())
        assert (neg, pos) == exp, (code, neg, pos)
        counts.append(dict(disease_code=code, disease=name, negative=neg, positive=pos, total=neg + pos))
    counts = pd.DataFrame(counts)
    assert counts.negative.sum() == 258 and counts.positive.sum() == 33
    assert d.loc[d['疾病编码'].eq('351'), 'state'].eq('LIMITED_ESTIMATE').all()
    domain_rows = []
    main_expected = dict(NMR_METABOLOMICS=53, HEALTH_BEHAVIOUR=8, ADULT_BODY_COMPOSITION=7, CLINICAL_BIOCHEMISTRY=5, BLOOD_PRESSURE=3, EDUCATION_SOCIOECONOMIC=2, OCCUPATION=1, PSYCHOSOCIAL_SOCIAL=1)
    for sid, label, arc, colour in DOMAINS:
        n = int(m.display_domain_id.eq(sid).sum())
        nm = int((m.display_domain_id.eq(sid) & m.marker_id.isin(mids)).sum())
        assert nm == main_expected.get(sid, 0)
        domain_rows.append(dict(domain_id=sid, display_label=label, arc_label=arc, colour=colour, main80_markers=nm, supp140_markers=n, source_domain_markers=int(m.source_domain.eq(sid).sum())))
    assert m.display_domain_id.eq('HEALTH_BEHAVIOUR').sum() == 13
    assert m.display_domain_id.eq('OCCUPATION').sum() == 2
    original = raw.set_index('候选ID').loc[d['候选ID'], raw.columns.drop('候选ID')]
    current = d.set_index('候选ID')[raw.columns.drop('候选ID')]
    pd.testing.assert_frame_equal(current, original, check_dtype=True)
    save_csv(d, results / 'CONTINUOUS_840_REAL_DATA_PLOT.csv')
    save_csv(raw[raw['类型'].ne('连续路径')], results / 'CATEGORICAL_240_EXCLUDED.csv')
    save_csv(m[m.marker_id.isin(mids)], results / 'MAIN_80_MARKER_ORDER.csv')
    save_csv(m, results / 'SUPP_140_MARKER_ORDER.csv')
    save_csv(counts, results / 'DISEASE_SUPPORTED_COUNTS.csv')
    save_csv(pd.DataFrame(domain_rows), results / 'DOMAIN_MAPPING_AND_COUNTS.csv')
    save_csv(pd.DataFrame([dict(state=k, count=int(d.state.eq(k).sum()), role='Frozen support/status audit only; not used to choose cell fill') for k in STATE_COLOURS]), results / 'SUPPORT_STATE_AUDIT_ONLY.csv')
    sample = pd.concat([d.sample(n=10, random_state=20261022), d.groupby('state', sort=False).sample(n=2, random_state=20261022)])
    sample = sample.drop_duplicates('候选ID').copy()
    sample['source_values_exact'] = [all((str(row[c]) == str(raw.loc[int(row['source_row']) - 1, c]) for c in raw.columns)) for row in sample.to_dict('records')]
    sample['star_expected_from_original_flag'] = sample['原批次证据'].eq('EXPLORATORY_PATH_SUPPORTED')
    assert sample.source_values_exact.all()
    assert sample.star_expected_from_original_flag.eq(sample.supported).all()
    save_csv(sample, results / 'RANDOM_PATH_SPOTCHECK.csv')
    return (d, m, counts, scale)

def polar(r, theta):
    a = math.radians(theta)
    return (r * math.cos(a), r * math.sin(a))

def geometry(d, markers):
    """Budget every degree: 270 = cells + inter-domain gaps + boundary margins."""
    present = [sid for sid, *_ in DOMAINS if markers.display_domain_id.eq(sid).any()]
    gap, margin = (1.05, 0.65)
    width = (270 - gap * (len(present) - 1) - 2 * margin) / len(markers)
    theta = 180 + margin
    slots, domains = ({}, [])
    for sid in present:
        start = theta
        block = markers[markers.display_domain_id.eq(sid)]
        for row in block.itertuples():
            slots[row.marker_id] = (theta, theta + width)
            theta += width
        domains.append(dict(domain_id=sid, theta_start=start, theta_end=theta, markers=len(block)))
        theta += gap
    assert abs(theta - gap - (450 - margin)) < 1e-08
    rows = []
    for row in d.to_dict('records'):
        lo, hi = slots[row['marker_id']]
        ring = int(row['ring_order'])
        outer = 0.922 - (ring - 1) * 0.105
        inner = outer - 0.095
        radius = (outer + inner) / 2
        x, y = polar(radius, (lo + hi) / 2)
        rows.append({**row, 'theta_start': lo, 'theta_end': hi, 'theta_mid': (lo + hi) / 2, 'r_inner': inner, 'r_outer': outer, 'x': x, 'y': y})
    return (pd.DataFrame(rows), domains, width)

def draw_info(ax, counts, present, mode, info_artists, scale):
    """Latest user reference overrides circular clipping: a square upper-left card.

    The square extends outside the circle only in the unused quadrant. Its lower
    and right bounds remain >0 and <0 respectively, so no 270-degree data cell
    or supported star can be covered. Legends are below the bars, not beside them.
    """
    frame = FancyBboxPatch((-1.035, 0.015), 1.02, 1.02, boxstyle='round,pad=0,rounding_size=.007', facecolor='#FBFCFE', edgecolor='#B9CEDB', linewidth=1.1, zorder=1)
    ax.add_patch(frame)

    def txt(x, y, s, size=8, **kwargs):
        artist = ax.text(x, y, s, fontsize=size, color=INK, va='center', zorder=6, **kwargs)
        artist.set_gid(f'information-{mode}-text-{len(info_artists):03d}')
        info_artists.append(artist)
        return artist

    def rect(x, y, w, h, colour, gid):
        p = Rectangle((x, y), w, h, facecolor=colour, edgecolor='none', zorder=4)
        p.set_gid(gid)
        ax.add_patch(p)
        info_artists.append(p)
        return p
    txt(-0.985, 0.992, 'Supported pathways by disease', 15, weight='bold')
    for x, state, label in [(-0.985, 'SUPPORTED_NEGATIVE', 'Supported · negative'), (-0.525, 'SUPPORTED_POSITIVE', 'Supported · positive')]:
        rect(x, 0.919, 0.035, 0.035, STATE_COLOURS[state], 'bar-legend-' + state)
        txt(x + 0.057, 0.936, label, 9.0)
    txt(-0.065, 0.936, 'Total', 8.7, ha='center', weight='bold')
    x0, span, height = (-0.675, 0.535, 0.052)
    ys = np.array([0.922 - k * 0.105 - 0.095 / 2 for k in range(6)])
    for tick in (0, 20, 40, 60, 80):
        x = x0 + span * tick / 80
        line, = ax.plot([x, x], [0.3, 0.909], color='#E6ECF0', linewidth=0.5, zorder=2)
        info_artists.append(line)
        txt(x, 0.278, str(tick), 7.8, ha='center')
    labels = ['DM', 'T2D-related', 'IHD', 'Angina', 'Coronary\natherosclerosis', 'Peripheral nerve']
    for k, (row, y) in enumerate(zip(counts.itertuples(), ys)):
        txt(-0.985, y, f'{k + 1}  {labels[k]}', 9.6, linespacing=1.04)
        for side, start, n in (('negative', 0, row.negative), ('positive', row.negative, row.positive)):
            if n:
                rect(x0 + span * start / 80, y - height / 2, span * n / 80, height, STATE_COLOURS['SUPPORTED_' + side.upper()], f'bar-{row.disease_code}-{side}-{n}')
                artist = txt(x0 + span * (start + n / 2) / 80, y, str(n), 8.0 if n >= 6 else 7.2, ha='center', weight='bold')
                artist.set_color('white')
        if row.total == 0:
            txt(x0 + 0.015, y, '0', 8.5)
        txt(-0.065, y, str(row.total), 11.1, ha='center', weight='bold')
    baseline, = ax.plot([x0, x0 + span], [0.3, 0.3], color=INK, linewidth=0.55, zorder=3)
    info_artists.append(baseline)
    txt(x0 + span / 2, 0.246, 'Number of supported pathways', 8.5, ha='center')
    divider, = ax.plot([-0.997, -0.053], [0.216, 0.216], color='#DAE4EB', linewidth=0.7)
    info_artists.append(divider)
    txt(-0.985, 0.192, 'a × b', 9.0, weight='bold')
    txt(-0.505, 0.192, 'Marker domains', 9.0, weight='bold')
    txt(-0.985, 0.166, 'Cox log-hazard coefficient product', 6.4)
    bar_x, bar_y, bar_w, bar_h = (-0.98, 0.115, 0.414, 0.026)
    for i in range(256):
        rect(bar_x + bar_w * i / 256, bar_y, bar_w / 256 + 2e-05, bar_h, to_hex(EFFECT_CMAP((i + 0.5) / 256)), f'effect-colorbar-{mode}-{i:03d}')
    for side, xs, colour in [('low', [bar_x, bar_x - 0.014, bar_x], EFFECT_ANCHORS[0]), ('high', [bar_x + bar_w, bar_x + bar_w + 0.014, bar_x + bar_w], EFFECT_ANCHORS[-1])]:
        p = Polygon(list(zip(xs, [bar_y, bar_y + bar_h / 2, bar_y + bar_h])), facecolor=colour, edgecolor='none', zorder=4)
        p.set_gid(f'effect-colorbar-extend-{mode}-{side}')
        ax.add_patch(p)
        info_artists.append(p)
    for x, value in [(bar_x, -scale['c']), (bar_x + bar_w / 2, 0), (bar_x + bar_w, scale['c'])]:
        txt(x, 0.094, '0' if value == 0 else f'{value:+.4f}', 6.7, ha='center')
    rect(-0.985, 0.049, 0.025, 0.022, NO_ESTIMATE_COLOUR, 'no-estimate-legend')
    txt(-0.944, 0.06, 'No finite estimate', 7.4)
    nrows = math.ceil(len(present) / 2)
    dy = 0.103 / max(nrows - 1, 1)
    for idx, domain in enumerate(present):
        _, label, _, colour = domain
        col, j = divmod(idx, nrows)
        x, y = (-0.505 + col * 0.238, 0.155 - j * dy)
        rect(x, y - 0.007, 0.014, 0.014, colour, f'domain-legend-{domain[0]}')
        txt(x + 0.023, y, label, 6.6)
    txt(-0.985, 0.028, '*: q < 0.05', 7.4)

def polygons_overlap(first, second, tolerance=1e-10):
    """Separating-axis test for actual rotated text quadrilaterals, not AABBs."""
    a, b = (np.asarray(first), np.asarray(second))
    for poly in (a, b):
        for i in range(len(poly)):
            edge = poly[(i + 1) % len(poly)] - poly[i]
            axis = np.array([-edge[1], edge[0]])
            aa, bb = (a @ axis, b @ axis)
            if aa.max() <= bb.min() + tolerance or bb.max() <= aa.min() + tolerance:
                return False
    return True

def radial_labels(ax, fig, cells, markers, label_map, mode):
    selected = markers.merge(label_map, on='marker_id', validate='one_to_one', suffixes=('', '_label'))
    assert len(selected) == len(markers)
    slots = cells.drop_duplicates('marker_id').set_index('marker_id')
    artists = []
    size = 9.6
    for i, row in enumerate(selected.to_dict('records'), 1):
        theta = float(slots.loc[row['marker_id'], 'theta_mid'])
        physical = theta % 360
        left = 90 < physical < 270
        rotation = physical + 180 if left else physical
        x, y = polar(1.015, theta)
        text = ax.text(x, y, row['display_label'], ha='right' if left else 'left', va='center', fontsize=size, weight='normal', color='#000000', rotation=rotation, rotation_mode='anchor', zorder=8)
        text.set_gid(f'marker-name-{mode}-{i:03d}')
        artists.append((row, theta, rotation, text))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    records, polygons = ([], [])
    for i, (row, theta, rotation, artist) in enumerate(artists, 1):
        anchor = artist.get_transform().transform(artist.get_position())
        artist.set_rotation(0)
        box = artist.get_window_extent(renderer)
        p = np.array([[box.x0, box.y0], [box.x1, box.y0], [box.x1, box.y1], [box.x0, box.y1]])
        rad = math.radians(rotation)
        mat = np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])
        p = (p - anchor) @ mat.T + anchor
        artist.set_rotation(rotation)
        norm = fig.transFigure.inverted().transform(p)
        corners_data = ax.transData.inverted().transform(p)
        minimum_radius = float('inf')
        for v, w in zip(corners_data, np.roll(corners_data, -1, axis=0)):
            edge = w - v
            projection = np.clip(-np.dot(v, edge) / np.dot(edge, edge), 0.0, 1.0)
            minimum_radius = min(minimum_radius, float(np.linalg.norm(v + projection * edge)))
        points_per_unit = np.linalg.norm(ax.transData.transform((1.0, 0.0)) - ax.transData.transform((0.0, 0.0))) * 72 / fig.dpi
        assert minimum_radius > 1.0, ('Outer label enters domain band', row['marker_id'], minimum_radius)
        polygons.append(norm)
        assert ((norm >= 0) & (norm <= 1)).all(), f"Radial label outside page: {row['display_label']}"
        rec = dict(marker_id=row['marker_id'], domain_id=row['display_domain_id'], source_domain=row['source_domain'], index=i, full_label=row['full_label'], display_label=row['display_label'], theta=theta, rotation=rotation, radius=1.015, fontsize_pt=artist.get_fontsize(), minimum_bbox_radius=minimum_radius, band_clearance_pt=(minimum_radius - 1) * points_per_unit, label_x=artist.get_position()[0], label_y=artist.get_position()[1], svg_id=artist.get_gid(), ha=artist.get_ha(), within_page=True)
        for j, corner in enumerate(norm):
            rec[f'corner{j + 1}_fx'], rec[f'corner{j + 1}_fy'] = corner
        records.append(rec)
    collisions = []
    for i in range(len(polygons)):
        for j in range(i + 1, len(polygons)):
            if polygons_overlap(polygons[i], polygons[j]):
                collisions.append((records[i]['marker_id'], records[j]['marker_id']))
    assert not collisions, f'Overlapping radial labels: {collisions}'
    card_data = np.array([[-1.035, 0.015], [-0.015, 0.015], [-0.015, 1.035], [-1.035, 1.035]])
    card_poly = fig.transFigure.inverted().transform(ax.transData.transform(card_data))
    assert not any((polygons_overlap(poly, card_poly) for poly in polygons)), 'Radial label overlaps card'
    assert len(records) == len(markers)
    return (pd.DataFrame(records), len(collisions))

def curved_domain_label(ax, spec, rec, point_per_unit, mode):
    """Editable glyphs on a circular baseline; every glyph's BOTTOM faces out.

    Counterclockwise character placement and rotation theta+90 give an inward
    up-vector and outward down-vector. Reverse traversal as well as rotation
    so the words remain correctly ordered rather than reading backwards.
    """
    label, radius, size = (spec[2], 0.974, 9.2)
    font = FontProperties(family='Arial', weight='bold', size=size)
    measure = TextToPath()
    prefix = [0.0] + [measure.get_text_width_height_descent(label[:i], font, False)[0] for i in range(1, len(label) + 1)]
    width = prefix[-1]
    arc_length = radius * math.radians(rec['theta_end'] - rec['theta_start']) * point_per_unit
    if arc_length < width * 1.18:
        return ([], dict(domain_id=spec[0], label=label, drawn=False, reason='Full name retained in legend; sector too short', text_width_pt=width, available_arc_pt=arc_length))
    mid = (rec['theta_start'] + rec['theta_end']) / 2
    glyphs = []
    for i, char in enumerate(label):
        along = (prefix[i] + prefix[i + 1]) / 2 - width / 2
        theta = mid + math.degrees(along / (radius * point_per_unit))
        x, y = polar(radius, theta)
        rotation = theta + 90
        if not char.isspace():
            artist = ax.text(x, y, char, ha='center', va='center', rotation=rotation, rotation_mode='anchor', fontproperties=font, color='#FFFFFF', zorder=4)
            artist.set_gid(f'domain-curved-{mode}-{spec[0]}-{i + 1:02d}')
        up = (-math.sin(math.radians(rotation)), math.cos(math.radians(rotation)))
        radial = (math.cos(math.radians(theta)), math.sin(math.radians(theta)))
        dot = up[0] * radial[0] + up[1] * radial[1]
        assert dot < -0.999999999999
        glyphs.append(dict(mode=mode, domain_id=spec[0], label=label, glyph_index=i + 1, character=char, x=x, y=y, radius=radius, theta=theta, rotation=rotation, fontsize_pt=size, up_dot_outward=dot, bottom_dot_outward=-dot))
    assert all((glyphs[i + 1]['theta'] > glyphs[i]['theta'] for i in range(len(glyphs) - 1)))
    return (glyphs, dict(domain_id=spec[0], label=label, drawn=True, reason='Curved counterclockwise baseline; all glyph bottoms face outward', text_width_pt=width, available_arc_pt=arc_length))

def render(d, markers, counts, out, mode, dpi, label_map, scale):
    if mode == 'main':
        mids = set(d.loc[d.supported, 'marker_id'])
        markers = markers[markers.marker_id.isin(mids)].copy()
        d = d[d.marker_id.isin(mids)].copy()
        stem = 'Figure_AB_Main80_UnifiedDomains_BottomOutward'
    else:
        stem = 'Figure_AB_Supp140_UnifiedDomains_BottomOutward'
    cells, domains, slotwidth = geometry(d, markers)
    assert len(cells) == (480 if mode == 'main' else 840)
    assert len(domains) == (8 if mode == 'main' else 12)
    assert int(cells.finite_ab.sum()) == (479 if mode == 'main' else 532)
    if mode == 'main':
        assert cells.loc[~cells.finite_ab, ['marker_id', '疾病编码']].values.tolist() == [['clinical__hba1c', '351']]
    inches, limit = (14.0, 1.6) if mode == 'main' else (18.0, 1.75)
    fig, ax = plt.subplots(figsize=(inches, inches), dpi=120)
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect='equal')
    ax.axis('off')
    for k in range(6):
        outer = 0.922 - k * 0.105
        ax.add_patch(Wedge((0, 0), outer, 180, 450, width=0.095, fc='#F3F5F5', ec='none', zorder=1))
    starrows = []
    for i, row in enumerate(cells.to_dict('records')):
        p = Wedge((0, 0), row['r_outer'], row['theta_start'], row['theta_end'], width=row['r_outer'] - row['r_inner'], facecolor=row['fill'], edgecolor='white', linewidth=0.25 if mode == 'main' else 0.16, zorder=2)
        p.set_gid(f'cell-{mode}-{i + 1}')
        ax.add_patch(p)
    point_per_unit = inches * 0.95 / (2 * limit) * 72
    min_radius = cells.loc[cells.supported, 'r_inner'].min()
    max_diameter = 0.62 * min_radius * math.radians(slotwidth) * point_per_unit
    diameter = min(4.3 if mode == 'main' else 3.0, max_diameter)
    star_radius = diameter / 2 / point_per_unit
    star_path = MPath.unit_regular_star(5)
    for i, row in enumerate(cells[cells.supported].to_dict('records')):
        artist = ax.scatter([row['x']], [row['y']], marker=star_path, s=diameter ** 2, facecolor=row['star_color'], edgecolor='none', zorder=5)
        artist.set_gid(f'supported-star-{mode}-{i + 1:03d}')
        assert row['r_inner'] < math.hypot(row['x'], row['y']) - star_radius
        assert math.hypot(row['x'], row['y']) + star_radius < row['r_outer']
        assert math.degrees(math.asin(star_radius / math.hypot(row['x'], row['y']))) < slotwidth / 2
        starrows.append({**row, 'star_diameter_pt': diameter, 'star_fits_cell': True})
    assert len(starrows) == 291
    curved_glyphs, curved_words = ([], [])
    for sid, *_ in DOMAINS:
        z = [r for r in domains if r['domain_id'] == sid]
        if not z:
            continue
        rec = z[0]
        spec = next((x for x in DOMAINS if x[0] == sid))
        p = Wedge((0, 0), 1.0, rec['theta_start'], rec['theta_end'], width=0.052, facecolor=spec[3], edgecolor='white', linewidth=0.5, zorder=3)
        p.set_gid(f'domain-band-{mode}-{sid}')
        ax.add_patch(p)
        glyphs, word = curved_domain_label(ax, spec, rec, point_per_unit, mode)
        curved_glyphs.extend(glyphs)
        curved_words.append(word)
    save_csv(pd.DataFrame(curved_glyphs), out / 'results' / f'CURVED_DOMAIN_GLYPHS_{mode}.csv')
    save_csv(pd.DataFrame(curved_words), out / 'results' / f'CURVED_DOMAIN_LABELS_{mode}.csv')
    ax.text(0, -0.045, 'BIRTH-ORDER', ha='center', va='center', fontsize=17, weight='bold', color=NAVY)
    ax.text(0, -0.109, 'PATHWAY ATLAS', ha='center', va='center', fontsize=17, weight='bold', color=NAVY)
    ax.text(0, -0.173, 'Exploratory a × b products', ha='center', va='center', fontsize=8.7, color=INK)
    present = [x for x in DOMAINS if x[0] in {r['domain_id'] for r in domains}]
    info_artists = []
    draw_info(ax, counts, present, mode, info_artists, scale)
    links = []
    for k, code in enumerate(DISEASES):
        z = cells[cells['疾病编码'].eq(code)]
        rc = float((z.r_inner.iloc[0] + z.r_outer.iloc[0]) / 2)
        line, = ax.plot([-0.041, -0.012], [rc, rc], color='#607A88', linewidth=0.65, zorder=7)
        line.set_gid(f'bar-ring-link-{mode}-{code}')
        links.append(dict(disease_code=code, ring_order=k + 1, bar_y=rc, ring_center_y=rc, link_x0=-0.041, link_x1=-0.012, right_ring_number_drawn=False))
    save_csv(pd.DataFrame(links), out / 'results' / f'DISEASE_RING_LINKS_{mode}.csv')
    label_geometry, label_collisions = radial_labels(ax, fig, cells, markers, label_map, mode)
    save_csv(label_geometry, out / 'results' / f'LABEL_GEOMETRY_{mode}.csv')
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    info_checks = []
    for i, artist in enumerate(info_artists):
        box = artist.get_window_extent(renderer).transformed(ax.transData.inverted())
        corners = [(x, y) for x in (box.x0, box.x1) for y in (box.y0, box.y1)]
        inside = all((-1.035 <= x <= -0.015 and 0.015 <= y <= 1.035 for x, y in corners))
        label = artist.get_text() if hasattr(artist, 'get_text') else artist.get_gid()
        info_checks.append(dict(artist=i, label=label, x0=box.x0, y0=box.y0, x1=box.x1, y1=box.y1, inside_information_card=inside))
    save_csv(pd.DataFrame(info_checks), out / 'results' / f'INFO_QUADRANT_BOUNDS_{mode}.csv')
    bad = [x for x in info_checks if not x['inside_information_card']]
    assert not bad, f'Visible information area exceeds square card: {bad}'
    text_bounds = []
    for artist in ax.texts:
        box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())
        ok = box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= 1) and (box.y1 <= 1)
        text_bounds.append(dict(label=artist.get_text(), x0=box.x0, y0=box.y0, x1=box.x1, y1=box.y1, within_page=ok))
    assert all((x['within_page'] for x in text_bounds))
    save_csv(pd.DataFrame(text_bounds), out / 'results' / f'TEXT_PAGE_BOUNDS_{mode}.csv')
    save_csv(cells, out / 'results' / f'CELL_GEOMETRY_{mode}.csv')
    save_csv(cells, out / ('FIG_MAIN80_PLOTTING_DATA.csv' if mode == 'main' else 'FIG_SUPP140_PLOTTING_DATA.csv'))
    save_csv(pd.DataFrame(domains), out / 'results' / f'DOMAIN_GEOMETRY_{mode}.csv')
    save_csv(pd.DataFrame(starrows), out / 'results' / f'STAR_QC_291_{mode}.csv')
    if mode == 'main':
        save_csv(pd.DataFrame(starrows), out / 'results' / 'STAR_QC_291.csv')
    for ext in ('png', 'pdf', 'svg'):
        fig.savefig(out / 'outputs' / f'{stem}.{ext}', dpi=dpi, facecolor='white')
    fig.savefig(out / 'outputs' / f'{stem}_preview.png', dpi=150, facecolor='white')
    plt.close(fig)
    return dict(mode=mode, markers=len(markers), cells=len(cells), stars=len(starrows), colored=int(cells.finite_ab.sum()), no_estimate=int((~cells.finite_ab).sum()), finite_unsupported=int((cells.finite_ab & ~cells.supported).sum()), black_stars=int((cells.supported & cells.star_color.eq('#000000')).sum()), white_stars=int((cells.supported & cells.star_color.eq('#FFFFFF')).sum()), color_c=scale['c'], negative_stars=int(cells.state.eq('SUPPORTED_NEGATIVE').sum()), positive_stars=int(cells.state.eq('SUPPORTED_POSITIVE').sum()), information_layout='Latest user reference: square upper-left card; legends below bar plot', information_card_bounds=[-1.035, 0.015, -0.015, 1.035], data_arc_degrees=270, domains=len(domains), domain_text_style='True curved circular baseline; all glyph BOTTOMS outward; no half-plane flip', curved_domain_labels=sum((x['drawn'] for x in curved_words)), all_domain_glyph_bottoms_face_outward=True, star_diameter_pt=diameter, all_info_within_card=True, all_star_bounds_inside_cells=True, labels=len(label_geometry), label_collisions=label_collisions, figure_inches=inches, axes_limits=[-limit, limit], axes_position=list(ax.get_position().bounds), png_dpi=dpi, png_pixels=round(inches * dpi), svg_page_points=inches * 72)

def main():
    global DOMAINS
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input', type=Path, default=DEFAULT_SOURCE)
    p.add_argument('--mapping', type=Path, default=DEFAULT_MAP)
    p.add_argument('--labels', type=Path, default=DEFAULT_LABELS)
    p.add_argument('--out', type=Path, default=DEFAULT_OUT)
    p.add_argument('--style', type=Path, default=DEFAULT_STYLE)
    p.add_argument('--display-map', type=Path, default=DEFAULT_DISPLAY_MAP)
    p.add_argument('--figures', choices=('main', 'both'), default='both')
    p.add_argument('--dpi', type=int, default=600)
    args = p.parse_args()
    source, mapping, out = (args.input.resolve(), args.mapping.resolve(), args.out.resolve())
    style_path, display_path = (args.style.resolve(), args.display_map.resolve())
    if out.exists():
        raise FileExistsError(f'Refusing to overwrite output directory: {out}')
    for path in (source, mapping, style_path, display_path, args.labels):
        assert path.is_file(), f'Missing frozen aggregate/mapping path: {path}'
    print(f'INPUT={source}\nSCRIPT_ROOT={ROOT}\nSTYLE={style_path}\nDISPLAY_MAP={display_path}\nOUTPUT={out}', flush=True)
    style = pd.read_csv(style_path, dtype=str, keep_default_na=False)
    assert len(style) == 12 and style.display_domain_id.is_unique and style.display_label.is_unique
    style['display_order'] = pd.to_numeric(style.display_order)
    style = style.sort_values('display_order')
    shared = [('Early life', '#4C78A8'), ('Education & SES', '#2A9D8F'), ('Occupation & work', '#59A14F'), ('Lifestyle', '#E39D3D'), ('Psychosocial', '#9C6BB5'), ('Geographical mobility', '#607D8B')]
    assert list(zip(style.display_label.iloc[:6], style.colour.iloc[:6])) == shared
    DOMAINS = [(r.display_domain_id, r.display_label, r.arc_label, r.colour) for r in style.itertuples()]
    assert args.dpi >= 400
    font = font_manager.findfont('Arial', fallback_to_default=False)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False})
    for subdir in ('outputs', 'results', 'code', 'source', 'qa'):
        (out / subdir).mkdir(parents=True, exist_ok=False)
    before = digest(source)
    d, markers, counts, scale = load_data(source, mapping, out / 'results', display_path)
    label_path = args.labels.resolve()
    label_map = pd.read_csv(label_path, dtype=str, keep_default_na=False)
    assert len(label_map) == 140 and label_map.marker_id.is_unique
    assert label_map.display_label.is_unique and label_map.display_label.str.strip().ne('').all()
    assert set(label_map.marker_id) == set(markers.marker_id)
    check_names = label_map.merge(markers, on='marker_id', validate='one_to_one')
    assert check_names.full_label.eq(check_names.marker_label).all()
    assert check_names.source_domain_x.eq(check_names.source_domain_y).all()
    save_csv(label_map, out / 'results' / 'RADIAL_LABEL_MAPPING.csv')
    order_table = markers.merge(label_map, on='marker_id', validate='one_to_one', suffixes=('', '_label'))
    order_table['in_main80'] = order_table.marker_id.isin(set(d.loc[d.supported, 'marker_id']))
    save_csv(order_table, out / 'DOMAIN_AND_MARKER_ORDER.csv')
    (out / 'FIG_COLOR_SCALE_MANIFEST.json').write_text(json.dumps(scale, indent=2, ensure_ascii=False))
    qc_rows = []
    for scope, frame in [('main80', d[d.in_main80]), ('supp140', d)]:
        for (code, status), group in frame.groupby(['疾病编码', '估计状态'], sort=False):
            qc_rows.append(dict(scope=scope, disease_code=code, estimate_status=status, cells=len(group), finite=int(group.finite_ab.sum()), no_estimate=int((~group.finite_ab).sum()), supported=int(group.supported.sum()), finite_unsupported=int((group.finite_ab & ~group.supported).sum()), clipped_for_color=int(group.clipped_for_color.sum())))
    save_csv(pd.DataFrame(qc_rows), out / 'QC_COUNTS_AND_STATUS.csv')
    print(f'Frozen source verified: {before}; 1080/840/240, 291=258+33; main80/supp140', flush=True)
    render_qc = [render(d, markers, counts, out, 'main', args.dpi, label_map, scale)]
    if args.figures == 'both':
        render_qc.append(render(d, markers, counts, out, 'supp', args.dpi, label_map, scale))
    assert digest(source) == before
    shutil.copy2(source, out / 'source' / source.name)
    shutil.copy2(mapping, out / 'code' / 'MARKER_REFERENCE_MAPPING.csv')
    shutil.copy2(__file__, out / 'code' / Path(__file__).name)
    shutil.copy2(label_path, out / 'code' / 'RADIAL_LABEL_MAPPING.csv')
    shutil.copy2(style_path, out / 'code' / style_path.name)
    shutil.copy2(display_path, out / 'code' / display_path.name)
    manifest = dict(status='RENDERED_NUMERIC_GEOMETRY_PASS_PENDING_EXPORT_VISUAL_QA', task_id='BORDER_CIRCOS_BOTTOM_OUTWARD_20261022_V3', run_utc=datetime.now(timezone.utc).isoformat(), input=str(source), input_sha256=before, source_unchanged=True, source_shape=[1080, 28], code_sha256=digest(__file__), mapping=str(mapping), mapping_sha256=digest(mapping), label_mapping=str(label_path), label_mapping_sha256=digest(label_path), mapping_reuse='Frozen source classification preserved; separate shared forest/Circos display taxonomy with construct-based ordering', display_style=str(style_path), display_style_sha256=digest(style_path), display_mapping=str(display_path), display_mapping_sha256=digest(display_path), source_domains_changed=False, statistical_fields_changed=False, display_membership_changes=[dict(marker_id='commute_miles', source_domain='OCCUPATION', display_domain_id='HEALTH_BEHAVIOUR', reason='Match established forest Lifestyle display domain')], support_family='Original frozen M774 + M66; not posthoc M840', categorical_excluded=240, supported_paths=291, models_run=0, no_new_P_q_CI=True, color_scale=scale, latest_instruction_override='Keep unified domains, actual radial names, square card, numeric fills and frozen stars; domain-name glyphs follow arc with BOTTOM outward and TOP inward', geometry=render_qc, fonts=dict(family='Arial', file=font), software=dict(python=sys.version, matplotlib=matplotlib.__version__, pandas=pd.__version__, numpy=np.__version__, platform=platform.platform()))
    (out / 'REPRODUCIBILITY_MANIFEST.json').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))
    print(json.dumps(render_qc, indent=2), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Frozen aggregate-data rendering only. No models or P/q/CI calculations."""\nimport argparse\nimport hashlib\nimport json\nimport math\nimport os\nimport platform\nimport shutil\nimport sys\nfrom datetime import datetime, timezone\nfrom pathlib import Path\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib import font_manager\nfrom matplotlib.font_manager import FontProperties\nfrom matplotlib.textpath import TextToPath\nfrom matplotlib.patches import Circle, Rectangle, Wedge, FancyBboxPatch, Polygon\nfrom matplotlib.path import Path as MPath\nfrom matplotlib.colors import LinearSegmentedColormap, to_hex, to_rgb\nimport numpy as np\nimport pandas as pd\nROOT = Path(__file__).resolve().parent\nDRIVE_G = Path(birth_private_setting(\'PRIVATE_G_\') if os.name == \'nt\' else \'/mnt/g\')\nDEFAULT_SOURCE = Path(os.environ.get(\'BORDER_FIGURE_SOURCE\', str(DRIVE_G / \'UKB_results/birth/border_pathway_supplement_phase2_20261021_v1/delivery/BORDER_AB_INTEGRATED_20261021_V1/RESULTS_MAIN_1080_ZH.csv\')))\nDEFAULT_MAP = ROOT / \'MARKER_REFERENCE_MAPPING.csv\'\nDEFAULT_LABELS = ROOT / \'RADIAL_LABEL_MAPPING.csv\'\nDEFAULT_STYLE = ROOT / \'DOMAIN_STYLE_UNIFIED.csv\'\nDEFAULT_DISPLAY_MAP = ROOT / \'MARKER_DISPLAY_DOMAIN_UNIFIED.csv\'\nDEFAULT_OUT = Path(os.environ.get(\'BORDER_FIGURE_OUT\', str(DRIVE_G / \'UKB_results/birth/border_circos_bottom_out_20261022_v3\')))\nDISEASES = [\'250\', \'250.2\', \'411\', \'411.3\', \'411.4\', \'351\']\nSHORT = [\'DM\', \'T2D\', \'IHD\', \'Angina\', \'Cor. ather.\', \'351\']\nFULL = [\'DM\', \'T2D-related PheCode\', \'IHD\', \'Angina\', \'Coronary atherosclerosis-related PheCode\', \'Other peripheral nervous-system disorders\']\nEXPECTED = [(72, 6), (72, 6), (61, 3), (8, 0), (45, 18), (0, 0)]\nSTATE_COLOURS = {\'SUPPORTED_NEGATIVE\': \'#2C6E96\', \'SUPPORTED_POSITIVE\': \'#CE6C63\', \'ESTIMABLE_NOT_SUPPORTED\': \'#DEE7EB\', \'LIMITED_ESTIMATE\': \'#EAE5DC\', \'OTHER_REVIEW\': \'#D6CFD9\'}\nDOMAINS = []\nNAVY = \'#243E4C\'\nINK = \'#34434A\'\nEFFECT_ANCHORS = [\'#1F5D88\', \'#9CBFD6\', \'#FFFFFF\', \'#E7B1A8\', \'#AF4540\']\nEFFECT_CMAP = LinearSegmentedColormap.from_list(\'frozen_ab_blue_white_coral\', EFFECT_ANCHORS, N=1025)\nNO_ESTIMATE_COLOUR = \'#E3E5E7\'\n\ndef apply_effect_palette(d):\n    """Colour derives only from finite original product; never from CI/status/P."""\n    d[\'finite_ab\'] = np.isfinite(d.ab)\n    assert int(d.finite_ab.sum()) == 532, \'Finite count differs from frozen 532: stop, do not fabricate.\'\n    assert not (d.supported & ~d.finite_ab).any(), \'Original supported flag has nonfinite product.\'\n    c = float(np.quantile(np.abs(d.loc[d.finite_ab, \'ab\']), 0.98, method=\'linear\'))\n    assert math.isfinite(c) and c > 0\n    d[\'color_ab\'] = np.where(d.finite_ab, np.clip(d.ab, -c, c), np.nan)\n    d[\'clipped_for_color\'] = d.finite_ab & d.ab.abs().gt(c)\n    d[\'fill\'] = [to_hex(EFFECT_CMAP((value + c) / (2 * c))).upper() if finite else NO_ESTIMATE_COLOUR for value, finite in zip(d.color_ab, d.finite_ab)]\n    rgb = np.array([to_rgb(colour) for colour in d.fill])\n    linear = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)\n    luminance = linear @ np.array([0.2126, 0.7152, 0.0722])\n    d[\'relative_luminance\'] = luminance\n    d[\'contrast_black\'] = (luminance + 0.05) / 0.05\n    d[\'contrast_white\'] = 1.05 / (luminance + 0.05)\n    d[\'star_color\'] = [\'#FFFFFF\' if support and white > black else \'#000000\' if support else \'\' for white, black, support in zip(d.contrast_white, d.contrast_black, d.supported)]\n    d.loc[~d.supported, \'star_color\'] = \'\'\n    d[\'point_estimate_warning\'] = \'\'\n    review = d.finite_ab & d[\'估计状态\'].isin([\'BOOT_UNSTABLE\', \'FAILURE_REVIEW\', \'MODEL_LIMITATION\'])\n    d.loc[review, \'point_estimate_warning\'] = \'Finite exploratory point estimate; diagnostic/bootstrap limitation retained; colour is not reliable inference.\'\n    d.loc[~d.finite_ab, \'point_estimate_warning\'] = \'No finite product; not filled with zero.\'\n    scale = dict(quantile=0.98, quantile_method=\'NumPy linear; R type 7\', finite_source_n=int(d.finite_ab.sum()), center=0, c=c, vmin=-c, vmax=c, palette_anchors=EFFECT_ANCHORS, palette_steps=1025, missing_colour=NO_ESTIMATE_COLOUR, raw_min=float(d.loc[d.finite_ab, \'ab\'].min()), raw_max=float(d.loc[d.finite_ab, \'ab\'].max()), clipped_negative=int((d.finite_ab & d.ab.lt(-c)).sum()), clipped_positive=int((d.finite_ab & d.ab.gt(c)).sum()), original_ab_values_changed=False, applies_to=\'Both main80 and supplement140\', label=\'a × b (Cox log-hazard coefficient product)\', star_rule=\'Original-batch evidence flag only; supported + nonfinite is a stop condition\', contrast_rule=\'Choose black/white with greater contrast against actual sRGB hex fill; linearized relative luminance\')\n    return (d, scale)\n\ndef digest(path):\n    return hashlib.sha256(Path(path).read_bytes()).hexdigest()\n\ndef save_csv(df, path):\n    if path.exists():\n        raise FileExistsError(path)\n    df.to_csv(path, index=False, encoding=\'utf-8-sig\')\n\ndef load_data(source, mapping, results, display_mapping):\n    raw = pd.read_csv(source, dtype=str, keep_default_na=False)\n    assert raw.shape == (1080, 28), raw.shape\n    assert raw[\'候选ID\'].is_unique\n    assert (raw[\'类型\'] == \'分类关联_非连续乘积\').sum() == 240\n    d = raw.loc[raw[\'类型\'].eq(\'连续路径\')].copy()\n    d[\'source_row\'] = d.index + 1\n    assert len(d) == 840 and d[\'原指标ID\'].nunique() == 140\n    assert not d.duplicated([\'原指标ID\', \'疾病编码\']).any()\n    assert set(d[\'疾病编码\']) == set(DISEASES)\n    assert d.groupby(\'原指标ID\')[\'疾病编码\'].apply(lambda x: set(x) == set(DISEASES)).all()\n    for source_col, col in [(\'a乘b\', \'ab\'), (\'乘积CI下限\', \'ci_low\'), (\'乘积CI上限\', \'ci_high\'), (\'原批次联合BH_q\', \'frozen_q\')]:\n        d[col] = pd.to_numeric(d[source_col], errors=\'coerce\')\n    valid = np.isfinite(d[[\'ab\', \'ci_low\', \'ci_high\', \'frozen_q\']]).all(axis=1)\n    same_direction = d.ab.gt(0) & d.ci_low.gt(0) | d.ab.lt(0) & d.ci_high.lt(0)\n    flag = d[\'原批次证据\'].eq(\'EXPLORATORY_PATH_SUPPORTED\')\n    supported = flag\n    assert not (supported & ~np.isfinite(d.ab)).any(), \'Frozen supported position has no finite product.\'\n    assert (flag == flag & d[\'估计状态\'].eq(\'CI_AVAILABLE\') & d.frozen_q.lt(0.05) & valid & same_direction).all(), \'Frozen source audit conflict.\'\n    d[\'supported\'] = supported\n    d[\'state\'] = \'OTHER_REVIEW\'\n    d.loc[d[\'估计状态\'].isin([\'MODEL_LIMITATION\', \'BOOT_UNSTABLE\', \'FAILURE_REVIEW\']), \'state\'] = \'LIMITED_ESTIMATE\'\n    d.loc[d[\'估计状态\'].eq(\'CI_AVAILABLE\') & ~supported, \'state\'] = \'ESTIMABLE_NOT_SUPPORTED\'\n    d.loc[supported & d.ab.lt(0), \'state\'] = \'SUPPORTED_NEGATIVE\'\n    d.loc[supported & d.ab.gt(0), \'state\'] = \'SUPPORTED_POSITIVE\'\n    assert (d.state == \'OTHER_REVIEW\').sum() == 0, \'Unmapped state: explicit review required.\'\n    assert supported.sum() == 291\n    mids = set(d.loc[supported, \'原指标ID\'])\n    assert len(mids) == 80\n    d[\'in_main80\'] = d[\'原指标ID\'].isin(mids)\n    assert d.in_main80.sum() == 480\n    m = pd.read_csv(mapping, dtype=str, keep_default_na=False)\n    assert m.marker_id.is_unique and len(m) == 140\n    assert set(m.marker_id) == set(d[\'原指标ID\'])\n    source_domains = d.groupby(\'原指标ID\')[\'领域\'].agg(lambda x: list(set(x)))\n    assert source_domains.apply(len).eq(1).all()\n    assert all((source_domains[row.marker_id][0] == row.source_domain for row in m.itertuples()))\n    assert set(m.source_domain) == {x[0] for x in DOMAINS}\n    dm = pd.read_csv(display_mapping, dtype=str, keep_default_na=False)\n    assert len(dm) == 140 and dm.marker_id.is_unique\n    assert set(dm.marker_id) == set(m.marker_id)\n    assert set(dm.display_domain_id) == {row[0] for row in DOMAINS}\n    assert m.set_index(\'marker_id\').source_domain.eq(dm.set_index(\'marker_id\').source_domain).all()\n    changed = dm.loc[dm.source_domain.ne(dm.display_domain_id), \'marker_id\'].tolist()\n    assert changed == [\'commute_miles\'], changed\n    m[\'source_marker_order\'] = pd.to_numeric(m.marker_order)\n    m[\'legacy_domain_label\'] = m.domain\n    m[\'source_domain_order\'] = pd.to_numeric(m.domain_order)\n    add = dm[[\'marker_id\', \'source_domain\', \'display_domain_id\', \'display_label\', \'colour\', \'display_order\', \'change_reason\']].rename(columns={\'display_label\': \'display_domain\', \'colour\': \'display_domain_colour\', \'display_order\': \'display_domain_order\'})\n    m = m.merge(add, on=[\'marker_id\', \'source_domain\'], how=\'left\', validate=\'one_to_one\')\n    assert m.display_domain_id.notna().all()\n    m[\'display_domain_order\'] = pd.to_numeric(m.display_domain_order)\n    m = m.sort_values([\'display_domain_order\', \'source_marker_order\']).reset_index(drop=True)\n    m[\'marker_order\'] = np.arange(1, len(m) + 1)\n    m[\'domain\'] = m.display_domain\n    m[\'domain_order\'] = m.display_domain_order\n    m[\'within_domain_order\'] = m.groupby(\'display_domain_id\', sort=False).cumcount() + 1\n    d[\'marker_id\'] = d[\'原指标ID\']\n    d[\'source_domain\'] = d[\'领域\']\n    d[\'display_domain_id\'] = d.marker_id.map(m.set_index(\'marker_id\').display_domain_id)\n    d[\'display_domain\'] = d.marker_id.map(m.set_index(\'marker_id\').display_domain)\n    d[\'display_domain_colour\'] = d.marker_id.map(m.set_index(\'marker_id\').display_domain_colour)\n    d[\'display_change_reason\'] = d.marker_id.map(m.set_index(\'marker_id\').change_reason)\n    d[\'domain_id\'] = d.display_domain_id\n    order = dict(zip(m.marker_id, m.marker_order))\n    d[\'marker_order\'] = d.marker_id.map(order)\n    d[\'ring_order\'] = d[\'疾病编码\'].map({v: k + 1 for k, v in enumerate(DISEASES)})\n    d, scale = apply_effect_palette(d)\n    d = d.sort_values([\'marker_order\', \'ring_order\'])\n    counts = []\n    for code, name, exp in zip(DISEASES, FULL, EXPECTED):\n        z = d[d[\'疾病编码\'].eq(code)]\n        neg = int(z.state.eq(\'SUPPORTED_NEGATIVE\').sum())\n        pos = int(z.state.eq(\'SUPPORTED_POSITIVE\').sum())\n        assert (neg, pos) == exp, (code, neg, pos)\n        counts.append(dict(disease_code=code, disease=name, negative=neg, positive=pos, total=neg + pos))\n    counts = pd.DataFrame(counts)\n    assert counts.negative.sum() == 258 and counts.positive.sum() == 33\n    assert d.loc[d[\'疾病编码\'].eq(\'351\'), \'state\'].eq(\'LIMITED_ESTIMATE\').all()\n    domain_rows = []\n    main_expected = dict(NMR_METABOLOMICS=53, HEALTH_BEHAVIOUR=8, ADULT_BODY_COMPOSITION=7, CLINICAL_BIOCHEMISTRY=5, BLOOD_PRESSURE=3, EDUCATION_SOCIOECONOMIC=2, OCCUPATION=1, PSYCHOSOCIAL_SOCIAL=1)\n    for sid, label, arc, colour in DOMAINS:\n        n = int(m.display_domain_id.eq(sid).sum())\n        nm = int((m.display_domain_id.eq(sid) & m.marker_id.isin(mids)).sum())\n        assert nm == main_expected.get(sid, 0)\n        domain_rows.append(dict(domain_id=sid, display_label=label, arc_label=arc, colour=colour, main80_markers=nm, supp140_markers=n, source_domain_markers=int(m.source_domain.eq(sid).sum())))\n    assert m.display_domain_id.eq(\'HEALTH_BEHAVIOUR\').sum() == 13\n    assert m.display_domain_id.eq(\'OCCUPATION\').sum() == 2\n    original = raw.set_index(\'候选ID\').loc[d[\'候选ID\'], raw.columns.drop(\'候选ID\')]\n    current = d.set_index(\'候选ID\')[raw.columns.drop(\'候选ID\')]\n    pd.testing.assert_frame_equal(current, original, check_dtype=True)\n    save_csv(d, results / \'CONTINUOUS_840_REAL_DATA_PLOT.csv\')\n    save_csv(raw[raw[\'类型\'].ne(\'连续路径\')], results / \'CATEGORICAL_240_EXCLUDED.csv\')\n    save_csv(m[m.marker_id.isin(mids)], results / \'MAIN_80_MARKER_ORDER.csv\')\n    save_csv(m, results / \'SUPP_140_MARKER_ORDER.csv\')\n    save_csv(counts, results / \'DISEASE_SUPPORTED_COUNTS.csv\')\n    save_csv(pd.DataFrame(domain_rows), results / \'DOMAIN_MAPPING_AND_COUNTS.csv\')\n    save_csv(pd.DataFrame([dict(state=k, count=int(d.state.eq(k).sum()), role=\'Frozen support/status audit only; not used to choose cell fill\') for k in STATE_COLOURS]), results / \'SUPPORT_STATE_AUDIT_ONLY.csv\')\n    sample = pd.concat([d.sample(n=10, random_state=20261022), d.groupby(\'state\', sort=False).sample(n=2, random_state=20261022)])\n    sample = sample.drop_duplicates(\'候选ID\').copy()\n    sample[\'source_values_exact\'] = [all((str(row[c]) == str(raw.loc[int(row[\'source_row\']) - 1, c]) for c in raw.columns)) for row in sample.to_dict(\'records\')]\n    sample[\'star_expected_from_original_flag\'] = sample[\'原批次证据\'].eq(\'EXPLORATORY_PATH_SUPPORTED\')\n    assert sample.source_values_exact.all()\n    assert sample.star_expected_from_original_flag.eq(sample.supported).all()\n    save_csv(sample, results / \'RANDOM_PATH_SPOTCHECK.csv\')\n    return (d, m, counts, scale)\n\ndef polar(r, theta):\n    a = math.radians(theta)\n    return (r * math.cos(a), r * math.sin(a))\n\ndef geometry(d, markers):\n    """Budget every degree: 270 = cells + inter-domain gaps + boundary margins."""\n    present = [sid for sid, *_ in DOMAINS if markers.display_domain_id.eq(sid).any()]\n    gap, margin = (1.05, 0.65)\n    width = (270 - gap * (len(present) - 1) - 2 * margin) / len(markers)\n    theta = 180 + margin\n    slots, domains = ({}, [])\n    for sid in present:\n        start = theta\n        block = markers[markers.display_domain_id.eq(sid)]\n        for row in block.itertuples():\n            slots[row.marker_id] = (theta, theta + width)\n            theta += width\n        domains.append(dict(domain_id=sid, theta_start=start, theta_end=theta, markers=len(block)))\n        theta += gap\n    assert abs(theta - gap - (450 - margin)) < 1e-08\n    rows = []\n    for row in d.to_dict(\'records\'):\n        lo, hi = slots[row[\'marker_id\']]\n        ring = int(row[\'ring_order\'])\n        outer = 0.922 - (ring - 1) * 0.105\n        inner = outer - 0.095\n        radius = (outer + inner) / 2\n        x, y = polar(radius, (lo + hi) / 2)\n        rows.append({**row, \'theta_start\': lo, \'theta_end\': hi, \'theta_mid\': (lo + hi) / 2, \'r_inner\': inner, \'r_outer\': outer, \'x\': x, \'y\': y})\n    return (pd.DataFrame(rows), domains, width)\n\ndef draw_info(ax, counts, present, mode, info_artists, scale):\n    """Latest user reference overrides circular clipping: a square upper-left card.\n\n    The square extends outside the circle only in the unused quadrant. Its lower\n    and right bounds remain >0 and <0 respectively, so no 270-degree data cell\n    or supported star can be covered. Legends are below the bars, not beside them.\n    """\n    frame = FancyBboxPatch((-1.035, 0.015), 1.02, 1.02, boxstyle=\'round,pad=0,rounding_size=.007\', facecolor=\'#FBFCFE\', edgecolor=\'#B9CEDB\', linewidth=1.1, zorder=1)\n    ax.add_patch(frame)\n\n    def txt(x, y, s, size=8, **kwargs):\n        artist = ax.text(x, y, s, fontsize=size, color=INK, va=\'center\', zorder=6, **kwargs)\n        artist.set_gid(f\'information-{mode}-text-{len(info_artists):03d}\')\n        info_artists.append(artist)\n        return artist\n\n    def rect(x, y, w, h, colour, gid):\n        p = Rectangle((x, y), w, h, facecolor=colour, edgecolor=\'none\', zorder=4)\n        p.set_gid(gid)\n        ax.add_patch(p)\n        info_artists.append(p)\n        return p\n    txt(-0.985, 0.992, \'Supported pathways by disease\', 15, weight=\'bold\')\n    for x, state, label in [(-0.985, \'SUPPORTED_NEGATIVE\', \'Supported · negative\'), (-0.525, \'SUPPORTED_POSITIVE\', \'Supported · positive\')]:\n        rect(x, 0.919, 0.035, 0.035, STATE_COLOURS[state], \'bar-legend-\' + state)\n        txt(x + 0.057, 0.936, label, 9.0)\n    txt(-0.065, 0.936, \'Total\', 8.7, ha=\'center\', weight=\'bold\')\n    x0, span, height = (-0.675, 0.535, 0.052)\n    ys = np.array([0.922 - k * 0.105 - 0.095 / 2 for k in range(6)])\n    for tick in (0, 20, 40, 60, 80):\n        x = x0 + span * tick / 80\n        line, = ax.plot([x, x], [0.3, 0.909], color=\'#E6ECF0\', linewidth=0.5, zorder=2)\n        info_artists.append(line)\n        txt(x, 0.278, str(tick), 7.8, ha=\'center\')\n    labels = [\'DM\', \'T2D-related\', \'IHD\', \'Angina\', \'Coronary\\natherosclerosis\', \'Peripheral nerve\']\n    for k, (row, y) in enumerate(zip(counts.itertuples(), ys)):\n        txt(-0.985, y, f\'{k + 1}  {labels[k]}\', 9.6, linespacing=1.04)\n        for side, start, n in ((\'negative\', 0, row.negative), (\'positive\', row.negative, row.positive)):\n            if n:\n                rect(x0 + span * start / 80, y - height / 2, span * n / 80, height, STATE_COLOURS[\'SUPPORTED_\' + side.upper()], f\'bar-{row.disease_code}-{side}-{n}\')\n                artist = txt(x0 + span * (start + n / 2) / 80, y, str(n), 8.0 if n >= 6 else 7.2, ha=\'center\', weight=\'bold\')\n                artist.set_color(\'white\')\n        if row.total == 0:\n            txt(x0 + 0.015, y, \'0\', 8.5)\n        txt(-0.065, y, str(row.total), 11.1, ha=\'center\', weight=\'bold\')\n    baseline, = ax.plot([x0, x0 + span], [0.3, 0.3], color=INK, linewidth=0.55, zorder=3)\n    info_artists.append(baseline)\n    txt(x0 + span / 2, 0.246, \'Number of supported pathways\', 8.5, ha=\'center\')\n    divider, = ax.plot([-0.997, -0.053], [0.216, 0.216], color=\'#DAE4EB\', linewidth=0.7)\n    info_artists.append(divider)\n    txt(-0.985, 0.192, \'a × b\', 9.0, weight=\'bold\')\n    txt(-0.505, 0.192, \'Marker domains\', 9.0, weight=\'bold\')\n    txt(-0.985, 0.166, \'Cox log-hazard coefficient product\', 6.4)\n    bar_x, bar_y, bar_w, bar_h = (-0.98, 0.115, 0.414, 0.026)\n    for i in range(256):\n        rect(bar_x + bar_w * i / 256, bar_y, bar_w / 256 + 2e-05, bar_h, to_hex(EFFECT_CMAP((i + 0.5) / 256)), f\'effect-colorbar-{mode}-{i:03d}\')\n    for side, xs, colour in [(\'low\', [bar_x, bar_x - 0.014, bar_x], EFFECT_ANCHORS[0]), (\'high\', [bar_x + bar_w, bar_x + bar_w + 0.014, bar_x + bar_w], EFFECT_ANCHORS[-1])]:\n        p = Polygon(list(zip(xs, [bar_y, bar_y + bar_h / 2, bar_y + bar_h])), facecolor=colour, edgecolor=\'none\', zorder=4)\n        p.set_gid(f\'effect-colorbar-extend-{mode}-{side}\')\n        ax.add_patch(p)\n        info_artists.append(p)\n    for x, value in [(bar_x, -scale[\'c\']), (bar_x + bar_w / 2, 0), (bar_x + bar_w, scale[\'c\'])]:\n        txt(x, 0.094, \'0\' if value == 0 else f\'{value:+.4f}\', 6.7, ha=\'center\')\n    rect(-0.985, 0.049, 0.025, 0.022, NO_ESTIMATE_COLOUR, \'no-estimate-legend\')\n    txt(-0.944, 0.06, \'No finite estimate\', 7.4)\n    nrows = math.ceil(len(present) / 2)\n    dy = 0.103 / max(nrows - 1, 1)\n    for idx, domain in enumerate(present):\n        _, label, _, colour = domain\n        col, j = divmod(idx, nrows)\n        x, y = (-0.505 + col * 0.238, 0.155 - j * dy)\n        rect(x, y - 0.007, 0.014, 0.014, colour, f\'domain-legend-{domain[0]}\')\n        txt(x + 0.023, y, label, 6.6)\n    txt(-0.985, 0.028, \'*: q < 0.05\', 7.4)\n\ndef polygons_overlap(first, second, tolerance=1e-10):\n    """Separating-axis test for actual rotated text quadrilaterals, not AABBs."""\n    a, b = (np.asarray(first), np.asarray(second))\n    for poly in (a, b):\n        for i in range(len(poly)):\n            edge = poly[(i + 1) % len(poly)] - poly[i]\n            axis = np.array([-edge[1], edge[0]])\n            aa, bb = (a @ axis, b @ axis)\n            if aa.max() <= bb.min() + tolerance or bb.max() <= aa.min() + tolerance:\n                return False\n    return True\n\ndef radial_labels(ax, fig, cells, markers, label_map, mode):\n    selected = markers.merge(label_map, on=\'marker_id\', validate=\'one_to_one\', suffixes=(\'\', \'_label\'))\n    assert len(selected) == len(markers)\n    slots = cells.drop_duplicates(\'marker_id\').set_index(\'marker_id\')\n    artists = []\n    size = 9.6\n    for i, row in enumerate(selected.to_dict(\'records\'), 1):\n        theta = float(slots.loc[row[\'marker_id\'], \'theta_mid\'])\n        physical = theta % 360\n        left = 90 < physical < 270\n        rotation = physical + 180 if left else physical\n        x, y = polar(1.015, theta)\n        text = ax.text(x, y, row[\'display_label\'], ha=\'right\' if left else \'left\', va=\'center\', fontsize=size, weight=\'normal\', color=\'#000000\', rotation=rotation, rotation_mode=\'anchor\', zorder=8)\n        text.set_gid(f\'marker-name-{mode}-{i:03d}\')\n        artists.append((row, theta, rotation, text))\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    records, polygons = ([], [])\n    for i, (row, theta, rotation, artist) in enumerate(artists, 1):\n        anchor = artist.get_transform().transform(artist.get_position())\n        artist.set_rotation(0)\n        box = artist.get_window_extent(renderer)\n        p = np.array([[box.x0, box.y0], [box.x1, box.y0], [box.x1, box.y1], [box.x0, box.y1]])\n        rad = math.radians(rotation)\n        mat = np.array([[math.cos(rad), -math.sin(rad)], [math.sin(rad), math.cos(rad)]])\n        p = (p - anchor) @ mat.T + anchor\n        artist.set_rotation(rotation)\n        norm = fig.transFigure.inverted().transform(p)\n        corners_data = ax.transData.inverted().transform(p)\n        minimum_radius = float(\'inf\')\n        for v, w in zip(corners_data, np.roll(corners_data, -1, axis=0)):\n            edge = w - v\n            projection = np.clip(-np.dot(v, edge) / np.dot(edge, edge), 0.0, 1.0)\n            minimum_radius = min(minimum_radius, float(np.linalg.norm(v + projection * edge)))\n        points_per_unit = np.linalg.norm(ax.transData.transform((1.0, 0.0)) - ax.transData.transform((0.0, 0.0))) * 72 / fig.dpi\n        assert minimum_radius > 1.0, (\'Outer label enters domain band\', row[\'marker_id\'], minimum_radius)\n        polygons.append(norm)\n        assert ((norm >= 0) & (norm <= 1)).all(), f"Radial label outside page: {row[\'display_label\']}"\n        rec = dict(marker_id=row[\'marker_id\'], domain_id=row[\'display_domain_id\'], source_domain=row[\'source_domain\'], index=i, full_label=row[\'full_label\'], display_label=row[\'display_label\'], theta=theta, rotation=rotation, radius=1.015, fontsize_pt=artist.get_fontsize(), minimum_bbox_radius=minimum_radius, band_clearance_pt=(minimum_radius - 1) * points_per_unit, label_x=artist.get_position()[0], label_y=artist.get_position()[1], svg_id=artist.get_gid(), ha=artist.get_ha(), within_page=True)\n        for j, corner in enumerate(norm):\n            rec[f\'corner{j + 1}_fx\'], rec[f\'corner{j + 1}_fy\'] = corner\n        records.append(rec)\n    collisions = []\n    for i in range(len(polygons)):\n        for j in range(i + 1, len(polygons)):\n            if polygons_overlap(polygons[i], polygons[j]):\n                collisions.append((records[i][\'marker_id\'], records[j][\'marker_id\']))\n    assert not collisions, f\'Overlapping radial labels: {collisions}\'\n    card_data = np.array([[-1.035, 0.015], [-0.015, 0.015], [-0.015, 1.035], [-1.035, 1.035]])\n    card_poly = fig.transFigure.inverted().transform(ax.transData.transform(card_data))\n    assert not any((polygons_overlap(poly, card_poly) for poly in polygons)), \'Radial label overlaps card\'\n    assert len(records) == len(markers)\n    return (pd.DataFrame(records), len(collisions))\n\ndef curved_domain_label(ax, spec, rec, point_per_unit, mode):\n    """Editable glyphs on a circular baseline; every glyph\'s BOTTOM faces out.\n\n    Counterclockwise character placement and rotation theta+90 give an inward\n    up-vector and outward down-vector. Reverse traversal as well as rotation\n    so the words remain correctly ordered rather than reading backwards.\n    """\n    label, radius, size = (spec[2], 0.974, 9.2)\n    font = FontProperties(family=\'Arial\', weight=\'bold\', size=size)\n    measure = TextToPath()\n    prefix = [0.0] + [measure.get_text_width_height_descent(label[:i], font, False)[0] for i in range(1, len(label) + 1)]\n    width = prefix[-1]\n    arc_length = radius * math.radians(rec[\'theta_end\'] - rec[\'theta_start\']) * point_per_unit\n    if arc_length < width * 1.18:\n        return ([], dict(domain_id=spec[0], label=label, drawn=False, reason=\'Full name retained in legend; sector too short\', text_width_pt=width, available_arc_pt=arc_length))\n    mid = (rec[\'theta_start\'] + rec[\'theta_end\']) / 2\n    glyphs = []\n    for i, char in enumerate(label):\n        along = (prefix[i] + prefix[i + 1]) / 2 - width / 2\n        theta = mid + math.degrees(along / (radius * point_per_unit))\n        x, y = polar(radius, theta)\n        rotation = theta + 90\n        if not char.isspace():\n            artist = ax.text(x, y, char, ha=\'center\', va=\'center\', rotation=rotation, rotation_mode=\'anchor\', fontproperties=font, color=\'#FFFFFF\', zorder=4)\n            artist.set_gid(f\'domain-curved-{mode}-{spec[0]}-{i + 1:02d}\')\n        up = (-math.sin(math.radians(rotation)), math.cos(math.radians(rotation)))\n        radial = (math.cos(math.radians(theta)), math.sin(math.radians(theta)))\n        dot = up[0] * radial[0] + up[1] * radial[1]\n        assert dot < -0.999999999999\n        glyphs.append(dict(mode=mode, domain_id=spec[0], label=label, glyph_index=i + 1, character=char, x=x, y=y, radius=radius, theta=theta, rotation=rotation, fontsize_pt=size, up_dot_outward=dot, bottom_dot_outward=-dot))\n    assert all((glyphs[i + 1][\'theta\'] > glyphs[i][\'theta\'] for i in range(len(glyphs) - 1)))\n    return (glyphs, dict(domain_id=spec[0], label=label, drawn=True, reason=\'Curved counterclockwise baseline; all glyph bottoms face outward\', text_width_pt=width, available_arc_pt=arc_length))\n\ndef render(d, markers, counts, out, mode, dpi, label_map, scale):\n    if mode == \'main\':\n        mids = set(d.loc[d.supported, \'marker_id\'])\n        markers = markers[markers.marker_id.isin(mids)].copy()\n        d = d[d.marker_id.isin(mids)].copy()\n        stem = \'Figure_AB_Main80_UnifiedDomains_BottomOutward\'\n    else:\n        stem = \'Figure_AB_Supp140_UnifiedDomains_BottomOutward\'\n    cells, domains, slotwidth = geometry(d, markers)\n    assert len(cells) == (480 if mode == \'main\' else 840)\n    assert len(domains) == (8 if mode == \'main\' else 12)\n    assert int(cells.finite_ab.sum()) == (479 if mode == \'main\' else 532)\n    if mode == \'main\':\n        assert cells.loc[~cells.finite_ab, [\'marker_id\', \'疾病编码\']].values.tolist() == [[\'clinical__hba1c\', \'351\']]\n    inches, limit = (14.0, 1.6) if mode == \'main\' else (18.0, 1.75)\n    fig, ax = plt.subplots(figsize=(inches, inches), dpi=120)\n    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)\n    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect=\'equal\')\n    ax.axis(\'off\')\n    for k in range(6):\n        outer = 0.922 - k * 0.105\n        ax.add_patch(Wedge((0, 0), outer, 180, 450, width=0.095, fc=\'#F3F5F5\', ec=\'none\', zorder=1))\n    starrows = []\n    for i, row in enumerate(cells.to_dict(\'records\')):\n        p = Wedge((0, 0), row[\'r_outer\'], row[\'theta_start\'], row[\'theta_end\'], width=row[\'r_outer\'] - row[\'r_inner\'], facecolor=row[\'fill\'], edgecolor=\'white\', linewidth=0.25 if mode == \'main\' else 0.16, zorder=2)\n        p.set_gid(f\'cell-{mode}-{i + 1}\')\n        ax.add_patch(p)\n    point_per_unit = inches * 0.95 / (2 * limit) * 72\n    min_radius = cells.loc[cells.supported, \'r_inner\'].min()\n    max_diameter = 0.62 * min_radius * math.radians(slotwidth) * point_per_unit\n    diameter = min(4.3 if mode == \'main\' else 3.0, max_diameter)\n    star_radius = diameter / 2 / point_per_unit\n    star_path = MPath.unit_regular_star(5)\n    for i, row in enumerate(cells[cells.supported].to_dict(\'records\')):\n        artist = ax.scatter([row[\'x\']], [row[\'y\']], marker=star_path, s=diameter ** 2, facecolor=row[\'star_color\'], edgecolor=\'none\', zorder=5)\n        artist.set_gid(f\'supported-star-{mode}-{i + 1:03d}\')\n        assert row[\'r_inner\'] < math.hypot(row[\'x\'], row[\'y\']) - star_radius\n        assert math.hypot(row[\'x\'], row[\'y\']) + star_radius < row[\'r_outer\']\n        assert math.degrees(math.asin(star_radius / math.hypot(row[\'x\'], row[\'y\']))) < slotwidth / 2\n        starrows.append({**row, \'star_diameter_pt\': diameter, \'star_fits_cell\': True})\n    assert len(starrows) == 291\n    curved_glyphs, curved_words = ([], [])\n    for sid, *_ in DOMAINS:\n        z = [r for r in domains if r[\'domain_id\'] == sid]\n        if not z:\n            continue\n        rec = z[0]\n        spec = next((x for x in DOMAINS if x[0] == sid))\n        p = Wedge((0, 0), 1.0, rec[\'theta_start\'], rec[\'theta_end\'], width=0.052, facecolor=spec[3], edgecolor=\'white\', linewidth=0.5, zorder=3)\n        p.set_gid(f\'domain-band-{mode}-{sid}\')\n        ax.add_patch(p)\n        glyphs, word = curved_domain_label(ax, spec, rec, point_per_unit, mode)\n        curved_glyphs.extend(glyphs)\n        curved_words.append(word)\n    save_csv(pd.DataFrame(curved_glyphs), out / \'results\' / f\'CURVED_DOMAIN_GLYPHS_{mode}.csv\')\n    save_csv(pd.DataFrame(curved_words), out / \'results\' / f\'CURVED_DOMAIN_LABELS_{mode}.csv\')\n    ax.text(0, -0.045, \'BIRTH-ORDER\', ha=\'center\', va=\'center\', fontsize=17, weight=\'bold\', color=NAVY)\n    ax.text(0, -0.109, \'PATHWAY ATLAS\', ha=\'center\', va=\'center\', fontsize=17, weight=\'bold\', color=NAVY)\n    ax.text(0, -0.173, \'Exploratory a × b products\', ha=\'center\', va=\'center\', fontsize=8.7, color=INK)\n    present = [x for x in DOMAINS if x[0] in {r[\'domain_id\'] for r in domains}]\n    info_artists = []\n    draw_info(ax, counts, present, mode, info_artists, scale)\n    links = []\n    for k, code in enumerate(DISEASES):\n        z = cells[cells[\'疾病编码\'].eq(code)]\n        rc = float((z.r_inner.iloc[0] + z.r_outer.iloc[0]) / 2)\n        line, = ax.plot([-0.041, -0.012], [rc, rc], color=\'#607A88\', linewidth=0.65, zorder=7)\n        line.set_gid(f\'bar-ring-link-{mode}-{code}\')\n        links.append(dict(disease_code=code, ring_order=k + 1, bar_y=rc, ring_center_y=rc, link_x0=-0.041, link_x1=-0.012, right_ring_number_drawn=False))\n    save_csv(pd.DataFrame(links), out / \'results\' / f\'DISEASE_RING_LINKS_{mode}.csv\')\n    label_geometry, label_collisions = radial_labels(ax, fig, cells, markers, label_map, mode)\n    save_csv(label_geometry, out / \'results\' / f\'LABEL_GEOMETRY_{mode}.csv\')\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    info_checks = []\n    for i, artist in enumerate(info_artists):\n        box = artist.get_window_extent(renderer).transformed(ax.transData.inverted())\n        corners = [(x, y) for x in (box.x0, box.x1) for y in (box.y0, box.y1)]\n        inside = all((-1.035 <= x <= -0.015 and 0.015 <= y <= 1.035 for x, y in corners))\n        label = artist.get_text() if hasattr(artist, \'get_text\') else artist.get_gid()\n        info_checks.append(dict(artist=i, label=label, x0=box.x0, y0=box.y0, x1=box.x1, y1=box.y1, inside_information_card=inside))\n    save_csv(pd.DataFrame(info_checks), out / \'results\' / f\'INFO_QUADRANT_BOUNDS_{mode}.csv\')\n    bad = [x for x in info_checks if not x[\'inside_information_card\']]\n    assert not bad, f\'Visible information area exceeds square card: {bad}\'\n    text_bounds = []\n    for artist in ax.texts:\n        box = artist.get_window_extent(renderer).transformed(fig.transFigure.inverted())\n        ok = box.x0 >= 0 and box.y0 >= 0 and (box.x1 <= 1) and (box.y1 <= 1)\n        text_bounds.append(dict(label=artist.get_text(), x0=box.x0, y0=box.y0, x1=box.x1, y1=box.y1, within_page=ok))\n    assert all((x[\'within_page\'] for x in text_bounds))\n    save_csv(pd.DataFrame(text_bounds), out / \'results\' / f\'TEXT_PAGE_BOUNDS_{mode}.csv\')\n    save_csv(cells, out / \'results\' / f\'CELL_GEOMETRY_{mode}.csv\')\n    save_csv(cells, out / (\'FIG_MAIN80_PLOTTING_DATA.csv\' if mode == \'main\' else \'FIG_SUPP140_PLOTTING_DATA.csv\'))\n    save_csv(pd.DataFrame(domains), out / \'results\' / f\'DOMAIN_GEOMETRY_{mode}.csv\')\n    save_csv(pd.DataFrame(starrows), out / \'results\' / f\'STAR_QC_291_{mode}.csv\')\n    if mode == \'main\':\n        save_csv(pd.DataFrame(starrows), out / \'results\' / \'STAR_QC_291.csv\')\n    for ext in (\'png\', \'pdf\', \'svg\'):\n        fig.savefig(out / \'outputs\' / f\'{stem}.{ext}\', dpi=dpi, facecolor=\'white\')\n    fig.savefig(out / \'outputs\' / f\'{stem}_preview.png\', dpi=150, facecolor=\'white\')\n    plt.close(fig)\n    return dict(mode=mode, markers=len(markers), cells=len(cells), stars=len(starrows), colored=int(cells.finite_ab.sum()), no_estimate=int((~cells.finite_ab).sum()), finite_unsupported=int((cells.finite_ab & ~cells.supported).sum()), black_stars=int((cells.supported & cells.star_color.eq(\'#000000\')).sum()), white_stars=int((cells.supported & cells.star_color.eq(\'#FFFFFF\')).sum()), color_c=scale[\'c\'], negative_stars=int(cells.state.eq(\'SUPPORTED_NEGATIVE\').sum()), positive_stars=int(cells.state.eq(\'SUPPORTED_POSITIVE\').sum()), information_layout=\'Latest user reference: square upper-left card; legends below bar plot\', information_card_bounds=[-1.035, 0.015, -0.015, 1.035], data_arc_degrees=270, domains=len(domains), domain_text_style=\'True curved circular baseline; all glyph BOTTOMS outward; no half-plane flip\', curved_domain_labels=sum((x[\'drawn\'] for x in curved_words)), all_domain_glyph_bottoms_face_outward=True, star_diameter_pt=diameter, all_info_within_card=True, all_star_bounds_inside_cells=True, labels=len(label_geometry), label_collisions=label_collisions, figure_inches=inches, axes_limits=[-limit, limit], axes_position=list(ax.get_position().bounds), png_dpi=dpi, png_pixels=round(inches * dpi), svg_page_points=inches * 72)\n\ndef main():\n    global DOMAINS\n    p = argparse.ArgumentParser(description=__doc__)\n    p.add_argument(\'--input\', type=Path, default=DEFAULT_SOURCE)\n    p.add_argument(\'--mapping\', type=Path, default=DEFAULT_MAP)\n    p.add_argument(\'--labels\', type=Path, default=DEFAULT_LABELS)\n    p.add_argument(\'--out\', type=Path, default=DEFAULT_OUT)\n    p.add_argument(\'--style\', type=Path, default=DEFAULT_STYLE)\n    p.add_argument(\'--display-map\', type=Path, default=DEFAULT_DISPLAY_MAP)\n    p.add_argument(\'--figures\', choices=(\'main\', \'both\'), default=\'both\')\n    p.add_argument(\'--dpi\', type=int, default=600)\n    args = p.parse_args()\n    source, mapping, out = (args.input.resolve(), args.mapping.resolve(), args.out.resolve())\n    style_path, display_path = (args.style.resolve(), args.display_map.resolve())\n    if out.exists():\n        raise FileExistsError(f\'Refusing to overwrite output directory: {out}\')\n    for path in (source, mapping, style_path, display_path, args.labels):\n        assert path.is_file(), f\'Missing frozen aggregate/mapping path: {path}\'\n    print(f\'INPUT={source}\\nSCRIPT_ROOT={ROOT}\\nSTYLE={style_path}\\nDISPLAY_MAP={display_path}\\nOUTPUT={out}\', flush=True)\n    style = pd.read_csv(style_path, dtype=str, keep_default_na=False)\n    assert len(style) == 12 and style.display_domain_id.is_unique and style.display_label.is_unique\n    style[\'display_order\'] = pd.to_numeric(style.display_order)\n    style = style.sort_values(\'display_order\')\n    shared = [(\'Early life\', \'#4C78A8\'), (\'Education & SES\', \'#2A9D8F\'), (\'Occupation & work\', \'#59A14F\'), (\'Lifestyle\', \'#E39D3D\'), (\'Psychosocial\', \'#9C6BB5\'), (\'Geographical mobility\', \'#607D8B\')]\n    assert list(zip(style.display_label.iloc[:6], style.colour.iloc[:6])) == shared\n    DOMAINS = [(r.display_domain_id, r.display_label, r.arc_label, r.colour) for r in style.itertuples()]\n    assert args.dpi >= 400\n    font = font_manager.findfont(\'Arial\', fallback_to_default=False)\n    plt.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 8, \'svg.fonttype\': \'none\', \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'axes.unicode_minus\': False})\n    for subdir in (\'outputs\', \'results\', \'code\', \'source\', \'qa\'):\n        (out / subdir).mkdir(parents=True, exist_ok=False)\n    before = digest(source)\n    d, markers, counts, scale = load_data(source, mapping, out / \'results\', display_path)\n    label_path = args.labels.resolve()\n    label_map = pd.read_csv(label_path, dtype=str, keep_default_na=False)\n    assert len(label_map) == 140 and label_map.marker_id.is_unique\n    assert label_map.display_label.is_unique and label_map.display_label.str.strip().ne(\'\').all()\n    assert set(label_map.marker_id) == set(markers.marker_id)\n    check_names = label_map.merge(markers, on=\'marker_id\', validate=\'one_to_one\')\n    assert check_names.full_label.eq(check_names.marker_label).all()\n    assert check_names.source_domain_x.eq(check_names.source_domain_y).all()\n    save_csv(label_map, out / \'results\' / \'RADIAL_LABEL_MAPPING.csv\')\n    order_table = markers.merge(label_map, on=\'marker_id\', validate=\'one_to_one\', suffixes=(\'\', \'_label\'))\n    order_table[\'in_main80\'] = order_table.marker_id.isin(set(d.loc[d.supported, \'marker_id\']))\n    save_csv(order_table, out / \'DOMAIN_AND_MARKER_ORDER.csv\')\n    (out / \'FIG_COLOR_SCALE_MANIFEST.json\').write_text(json.dumps(scale, indent=2, ensure_ascii=False))\n    qc_rows = []\n    for scope, frame in [(\'main80\', d[d.in_main80]), (\'supp140\', d)]:\n        for (code, status), group in frame.groupby([\'疾病编码\', \'估计状态\'], sort=False):\n            qc_rows.append(dict(scope=scope, disease_code=code, estimate_status=status, cells=len(group), finite=int(group.finite_ab.sum()), no_estimate=int((~group.finite_ab).sum()), supported=int(group.supported.sum()), finite_unsupported=int((group.finite_ab & ~group.supported).sum()), clipped_for_color=int(group.clipped_for_color.sum())))\n    save_csv(pd.DataFrame(qc_rows), out / \'QC_COUNTS_AND_STATUS.csv\')\n    print(f\'Frozen source verified: {before}; 1080/840/240, 291=258+33; main80/supp140\', flush=True)\n    render_qc = [render(d, markers, counts, out, \'main\', args.dpi, label_map, scale)]\n    if args.figures == \'both\':\n        render_qc.append(render(d, markers, counts, out, \'supp\', args.dpi, label_map, scale))\n    assert digest(source) == before\n    shutil.copy2(source, out / \'source\' / source.name)\n    shutil.copy2(mapping, out / \'code\' / \'MARKER_REFERENCE_MAPPING.csv\')\n    shutil.copy2(__file__, out / \'code\' / Path(__file__).name)\n    shutil.copy2(label_path, out / \'code\' / \'RADIAL_LABEL_MAPPING.csv\')\n    shutil.copy2(style_path, out / \'code\' / style_path.name)\n    shutil.copy2(display_path, out / \'code\' / display_path.name)\n    manifest = dict(status=\'RENDERED_NUMERIC_GEOMETRY_PASS_PENDING_EXPORT_VISUAL_QA\', task_id=\'BORDER_CIRCOS_BOTTOM_OUTWARD_20261022_V3\', run_utc=datetime.now(timezone.utc).isoformat(), input=str(source), input_sha256=before, source_unchanged=True, source_shape=[1080, 28], code_sha256=digest(__file__), mapping=str(mapping), mapping_sha256=digest(mapping), label_mapping=str(label_path), label_mapping_sha256=digest(label_path), mapping_reuse=\'Frozen source classification preserved; separate shared forest/Circos display taxonomy with construct-based ordering\', display_style=str(style_path), display_style_sha256=digest(style_path), display_mapping=str(display_path), display_mapping_sha256=digest(display_path), source_domains_changed=False, statistical_fields_changed=False, display_membership_changes=[dict(marker_id=\'commute_miles\', source_domain=\'OCCUPATION\', display_domain_id=\'HEALTH_BEHAVIOUR\', reason=\'Match established forest Lifestyle display domain\')], support_family=\'Original frozen M774 + M66; not posthoc M840\', categorical_excluded=240, supported_paths=291, models_run=0, no_new_P_q_CI=True, color_scale=scale, latest_instruction_override=\'Keep unified domains, actual radial names, square card, numeric fills and frozen stars; domain-name glyphs follow arc with BOTTOM outward and TOP inward\', geometry=render_qc, fonts=dict(family=\'Arial\', file=font), software=dict(python=sys.version, matplotlib=matplotlib.__version__, pandas=pd.__version__, numpy=np.__version__, platform=platform.platform()))\n    (out / \'REPRODUCIBILITY_MANIFEST.json\').write_text(json.dumps(manifest, ensure_ascii=False, indent=2))\n    print(json.dumps(render_qc, indent=2), flush=True)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
