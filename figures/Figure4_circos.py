"""Function library; no source-time data loading, downloads or analysis dispatch.
Private inputs and contexts must be provided explicitly.
"""
import argparse
import hashlib
import json
import math
import os
import shutil
import sys
from datetime import datetime
from pathlib import Path
from zoneinfo import ZoneInfo
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.pyplot as plt
from matplotlib.colors import to_hex, to_rgb
from matplotlib.patches import Wedge, Rectangle, PathPatch, Polygon
from matplotlib.path import Path as MPath
from matplotlib.text import Text
from matplotlib.transforms import Bbox
import pymupdf
import legacy_circos_geometry as legacy
ENDPOINTS = ['250', '250.2', '411', '411.3', '411.4']
SHORT = ['DM', 'T2D', 'IHD', 'Angina', 'Cor. atheroscl.']
Q_COLUMN = 'q_BH_global605_exploratory'
NEW_NAMES = {'alcohol': 'Alcohol status', 'body_size10': 'Childhood body size', 'breastfed': 'Breastfed', 'chronotype': 'Chronotype', 'clinical__ldl': 'LDL-C (clinical)', 'clinical__tg': 'Triglycerides', 'commute_car': 'Car commuting', 'commute_public_transport': 'Public transport', 'commute_walking': 'Walking to work', 'confiding': 'Confiding frequency', 'employment_employed': 'Paid/self-employment', 'employment_student': 'Student status', 'employment_voluntary': 'Voluntary work', 'height10': 'Childhood height', 'housing_tenure': 'Housing tenure', 'income': 'Household income', 'insomnia': 'Insomnia', 'living_alone': 'Living alone', 'loneliness': 'Loneliness', 'manual_work': 'Manual work', 'maternal_smoking': 'Maternal smoking', 'nap': 'Daytime naps', 'oily_fish': 'Oily fish', 'participation_adult_education': 'Adult-education group', 'participation_other': 'Other social group', 'participation_religious': 'Religious group', 'qualification': 'Education', 'smoking': 'Smoking status', 'vehicles': 'Vehicles', 'visits': 'Visits with friends/family'}
DOMAIN_MAP = {'Early life': 'EARLY_LIFE', 'Education & SES': 'EDUCATION_SOCIOECONOMIC', 'Employment & occupation': 'OCCUPATION', 'Lifestyle (excluding smoking and alcohol)': 'HEALTH_BEHAVIOUR', 'Psychosocial & social connections': 'PSYCHOSOCIAL_SOCIAL', 'Geographical mobility': 'GEOGRAPHY_ENVIRONMENT'}

def sha(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()

def jwrite(p, d):
    with p.open('x') as f:
        json.dump(d, f, ensure_ascii=False, indent=2, default=str, allow_nan=False)

def mapping(d, style_root, social_path):
    style = pd.read_csv(style_root / 'DOMAIN_STYLE_UNIFIED.csv').sort_values('display_order')
    old = pd.read_csv(style_root / 'MARKER_DISPLAY_DOMAIN_UNIFIED.csv').set_index('marker_id')
    names = pd.read_csv(style_root / 'RADIAL_LABEL_MAPPING.csv').set_index('marker_id')
    social = pd.read_csv(social_path).drop_duplicates('canonical_id').set_index('canonical_id')
    records = []
    for i, row in enumerate(d.drop_duplicates('indicator_id').itertuples()):
        ident = row.indicator_id
        if ident in old.index:
            domain = old.loc[ident, 'display_domain_id']
            basis = 'Reused approved forest/Circos display map'
        elif ident in social.index:
            domain = DOMAIN_MAP[social.loc[ident, 'domain']]
            basis = 'Existing nontraditional A-forest domain, not regrouped by effect'
        elif ident.startswith('clinical__'):
            domain = 'CLINICAL_BIOCHEMISTRY'
            basis = 'Existing clinical-assay domain'
        elif ident in {'smoking', 'alcohol'}:
            domain = 'HEALTH_BEHAVIOUR'
            basis = 'Existing lifestyle domain'
        else:
            raise ValueError(f'Unmapped indicator: {ident}')
        label = NEW_NAMES.get(ident)
        if label is None:
            assert ident in names.index, ident
            label = str(names.loc[ident, 'display_label'])
        if ident == 'clinical__hdl':
            label = 'HDL-C (clinical)'
        if ident == 'clinical__apoa':
            label = 'ApoA1 (clinical)'
        if ident == 'clinical__tc':
            label = 'Total cholesterol (clinical)'
        if ident == 'mobility':
            label = 'Birth–home distance'
        if ident == 'nmr__nmr_23404':
            label = 'LDL-C (NMR)'
        if ident == 'nmr__nmr_23406':
            label = 'HDL-C (NMR)'
        if ident in social.index:
            order = float(social.loc[ident, 'marker_order'])
        elif ident in old.index:
            order = float(old.loc[ident, 'source_marker_order'])
        else:
            order = 1000 + i
        records.append(dict(marker_id=ident, source_domain=row.domain, display_domain_id=domain, marker_order=order, source_indicator_name=row.indicator_name, full_label=str(social.loc[ident, 'label']) if ident in social.index else str(names.loc[ident, 'full_label']) if ident in names.index else label, display_label=label, domain_mapping_basis=basis, label_mapping_basis='Existing radial label or audited construct translation'))
    m = pd.DataFrame(records)
    meta = style.set_index('display_domain_id')
    m['domain_order'] = m.display_domain_id.map(meta.display_order)
    m['display_domain'] = m.display_domain_id.map(meta.display_label)
    m['display_colour'] = m.display_domain_id.map(meta.colour)
    assert m.domain_order.notna().all() and m.marker_id.is_unique and (len(m) == 121)
    assert m.display_label.is_unique, m.loc[m.display_label.duplicated(False)].to_dict('records')
    m = m.sort_values(['domain_order', 'marker_order', 'marker_id'], kind='stable').reset_index(drop=True)
    m['marker_order'] = np.arange(1, len(m) + 1)
    m['within_domain_order'] = m.groupby('display_domain_id', sort=False).cumcount() + 1
    return (m, style)

def colorize(d):
    d = d.copy()
    d['finite_ie'] = np.isfinite(d.IE_per_1000)
    assert d.finite_ie.sum() == 550
    d['supported'] = d[Q_COLUMN].lt(0.05)
    assert d.supported.sum() == 429 and (not (d.supported & ~d.finite_ie).any())
    assert (d.loc[d.supported, 'IE_CI_low'] * d.loc[d.supported, 'IE_CI_high']).gt(0).all()
    c = float(np.quantile(np.abs(d.loc[d.finite_ie, 'IE_per_1000']), 0.98, method='linear'))
    d['display_ie_clipped'] = d.IE_per_1000.clip(-c, c)
    d['display_color_clipped'] = d.finite_ie & d.IE_per_1000.abs().gt(c)
    d['fill'] = [to_hex(legacy.EFFECT_CMAP((x + c) / (2 * c))) if finite else legacy.NO_ESTIMATE_COLOUR for x, finite in zip(d.display_ie_clipped, d.finite_ie)]
    rgb = np.array([to_rgb(x) for x in d.fill])
    rgb_lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)
    lum = rgb_lin @ np.array([0.2126, 0.7152, 0.0722])
    d['relative_luminance'] = lum
    d['star_color'] = np.where(d.supported, np.where(1.05 / (lum + 0.05) > (lum + 0.05) / 0.05, '#FFFFFF', '#000000'), '')
    d['q_column_used'] = Q_COLUMN
    d['q_family_M'] = 605
    d['conditional_or_warned'] = d.result_status.isin(['PARTIAL_BOOTSTRAP_INFERENCE', 'WARNING_AWARE_PARTIAL_BOOTSTRAP_INFERENCE', 'WARNING_AWARE_BOOTSTRAP_INFERENCE'])
    return (d, dict(c=c, vmin=-c, vmax=c, center=0, finite_values=550, finite_min=float(d.IE_per_1000.min()), finite_max=float(d.IE_per_1000.max()), unit='11-year indirect risk difference per 1,000 participants', percentile=98, percentile_method='Type 7; display only', clipping_changes_source_values=False, missing_color=legacy.NO_ESTIMATE_COLOUR, palette=legacy.EFFECT_ANCHORS, star_rule=Q_COLUMN + ' < 0.05; fixed M=605', star_contrast='Higher black/white relative-luminance contrast against actual fill'))

def info_card(ax, counts, domains, scale, frame):
    """Outward information region with a true circular cut-out, not a mask."""
    left, top, right, bottom = (frame['left'], frame['top'], frame['right'], frame['bottom'])
    radius = frame['notch_radius']
    theta0 = math.degrees(math.acos(right / radius))
    theta1 = 180 - math.degrees(math.asin(bottom / radius))
    arc = MPath.arc(theta0, theta1, n=32)
    av = arc.vertices * radius
    vertices = np.vstack([[(left, top), (right, top), av[0]], av[1:], [(left, bottom), (left, top)]])
    codes = [MPath.MOVETO, MPath.LINETO, MPath.LINETO] + list(arc.codes[1:]) + [MPath.LINETO, MPath.CLOSEPOLY]
    box = PathPatch(MPath(vertices, codes), fc='#FBFCFE', ec='#B9CEDB', lw=1.1, zorder=1)
    box.set_gid('outward-information-card-with-circular-notch')
    ax.add_patch(box)
    artists = []

    def txt(x, y, text, size=8, **kw):
        weight = kw.pop('weight', 'normal')
        a = ax.text(x, y, text, fontsize=size * 1.3, color='#000000', weight=weight, va='center', zorder=7, **kw)
        artists.append(a)
        return a

    def rect(x, y, w, h, color):
        a = Rectangle((x, y), w, h, fc=color, ec='none', zorder=4)
        ax.add_patch(a)
        artists.append(a)
        return a
    text_left = left + 0.048
    title_y = top - 0.057
    legend_y = 0.98
    txt(text_left, title_y, 'Indirect effects by disease', 15, weight='bold')
    for x, label, color in [(text_left, 'Negative', '#2C6E96'), (text_left + 0.41, 'Positive', '#CE6C63')]:
        rect(x, legend_y - 0.015, 0.029, 0.029, color)
        txt(x + 0.047, legend_y, label, 15 / 1.3)
    txt(-0.065, legend_y, 'Total', 9.2, ha='center', weight='bold')
    x0, span, height = (-0.87, 0.735, 0.053)
    names = ['Diabetes mellitus', 'Type 2 diabetes', 'Ischemic heart\ndisease', 'Angina pectoris', 'Coronary\natherosclerosis']
    for t in range(0, 101, 20):
        x = x0 + span * t / 100
        ax.plot([x, x], [0.402, 0.908], color='#E5EBEF', lw=0.5, zorder=2)
        txt(x, 0.414, str(t), 8, ha='center')
    for i, row in enumerate(counts.itertuples()):
        y = 0.922 - i * 0.105 - 0.095 / 2
        txt(text_left, y, names[i], 18.5 / 1.3, linespacing=1.0)
        for start, n, color in [(0, row.negative, '#2C6E96'), (row.negative, row.positive, '#CE6C63')]:
            rect(x0 + span * start / 100, y - height / 2, span * n / 100, height, color)
            if n:
                a = txt(x0 + span * (start + n / 2) / 100, y, str(n), 9, ha='center')
                a.set_color('#FFFFFF')
                a.set_gid(f'bar-count-white-{i}-{start}')
        txt(-0.065, y, str(row.total), 11.2, ha='center')
    ax.plot([x0, x0 + span], [0.402, 0.402], color=legacy.INK, lw=0.6, zorder=4)
    txt(x0 + span / 2, 0.35, 'Number with BH q < 0.05', 8.8, ha='center')
    ax.plot([left + 0.025, -0.29], [0.325, 0.325], color='#DAE4EB', lw=0.7)
    txt(text_left, 0.294, 'IE per 1,000', 10, weight='bold')
    txt(text_left, 0.255, 'Risk difference', 8)
    bx, by, bw, bh = (text_left + 0.005, 0.208, 0.29, 0.027)
    for k in range(256):
        rect(bx + bw * k / 256, by, bw / 256 + 2e-05, bh, to_hex(legacy.EFFECT_CMAP((k + 0.5) / 256)))
    for xs, colour in [([bx, bx - 0.014, bx], legacy.EFFECT_ANCHORS[0]), ([bx + bw, bx + bw + 0.014, bx + bw], legacy.EFFECT_ANCHORS[-1])]:
        ax.add_patch(Polygon(list(zip(xs, [by, by + bh / 2, by + bh])), fc=colour, ec='none', zorder=4))
    for x, value in [(bx, -scale['c']), (bx + bw / 2, 0), (bx + bw, scale['c'])]:
        txt(x, 0.182, '0' if value == 0 else f'{value:+.2f}', 8, ha='center')
    txt(text_left, 0.139, '*: BH q < 0.05 (M = 605)', 8)
    rect(text_left, 0.082, 0.024, 0.023, legacy.NO_ESTIMATE_COLOUR)
    txt(text_left + 0.041, 0.093, 'No estimate', 8)
    txt(text_left, 0.048, '550 estimates; 55 unavailable', 8)
    txt(-0.955, 0.294, 'Domains', 10, weight='bold')
    nrows = math.ceil(len(domains) / 2)
    for k, (_, label, _, color) in enumerate(domains):
        col, j = divmod(k, nrows)
        x = -0.955 + col * 0.295
        y = 0.253 - j * 0.044
        rect(x, y - 0.008, 0.015, 0.015, color)
        txt(x + 0.025, y, label, 7.2)
    return (artists, dict(**frame, label_left=text_left, bar_start=x0, bar_width=span, previous_bar_width=0.535, bar_width_multiplier=span / 0.535, domain_legend_columns=[-0.955, -0.66], no_white_mask_used=True, card_path_vertices=vertices.tolist(), card_path_codes=[int(x) for x in codes]))

def render(d, m, style, counts, scale, out, dpi):
    present = set(m.display_domain_id)
    legacy.DOMAINS = [(r.display_domain_id, r.display_label, r.arc_label, r.colour) for r in style.itertuples() if r.display_domain_id in present]
    cells, domains, width = legacy.geometry(d, m)
    assert len(cells) == 605 and len(domains) == 11
    inches, limit = (16.0 * 1.55 / 1.43, 1.55)
    fig, ax = plt.subplots(figsize=(inches, inches), dpi=120)
    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)
    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect='equal')
    ax.axis('off')
    for k in range(5):
        ax.add_patch(Wedge((0, 0), 0.922 - k * 0.105, 180, 450, width=0.095, fc='#F3F5F5', ec='none', zorder=1))
    for i, row in enumerate(cells.itertuples()):
        a = Wedge((0, 0), row.r_outer, row.theta_start, row.theta_end, width=row.r_outer - row.r_inner, fc=row.fill, ec='white', lw=0.2, zorder=2)
        a.set_gid(f'cell-B100-{i + 1:03d}')
        ax.add_patch(a)
    ppu = inches * 0.95 / (2 * limit) * 72
    diameter = min(5.6, 0.84 * cells.loc[cells.supported, 'r_inner'].min() * math.radians(width) * ppu)
    star_radius = diameter / 2 / ppu
    starrows = []
    for i, row in enumerate(cells[cells.supported].to_dict('records')):
        a = ax.scatter(row['x'], row['y'], marker=MPath.unit_regular_star(5), s=diameter ** 2, c=row['star_color'], edgecolors='none', zorder=5)
        a.set_gid(f'B100-q-star-{i + 1:03d}')
        assert row['r_inner'] < math.hypot(row['x'], row['y']) - star_radius
        assert math.hypot(row['x'], row['y']) + star_radius < row['r_outer']
        assert math.degrees(math.asin(star_radius / math.hypot(row['x'], row['y']))) < width / 2
        starrows.append({**row, 'star_diameter_pt': diameter, 'star_fits_cell': True})
    assert len(starrows) == 429
    glyphs, words = ([], [])
    for rec in domains:
        spec = next((x for x in legacy.DOMAINS if x[0] == rec['domain_id']))
        a = Wedge((0, 0), 1, rec['theta_start'], rec['theta_end'], width=0.052, fc=spec[3], ec='white', lw=0.5, zorder=3)
        a.set_gid('domain-band-' + spec[0])
        ax.add_patch(a)
        gs, word = legacy.curved_domain_label(ax, spec, rec, ppu, 'B100')
        glyphs.extend(gs)
        words.append(word)
    center_texts = []
    for y, label, size, weight in [(0.125, 'BIRTH ORDER', 31.2, 'bold'), (0.035, 'INDIRECT EFFECTS', 29.9, 'bold'), (-0.07, 'Relative birth order (0 → 1)', 16.9, 'normal')]:
        center_texts.append(ax.text(0, y, label, ha='center', va='center', fontsize=size, weight=weight, color='#000000'))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    boxes = [a.get_window_extent(renderer).transformed(ax.transData.inverted()) for a in center_texts]
    dx = -(min((b.x0 for b in boxes)) + max((b.x1 for b in boxes))) / 2
    dy = -(min((b.y0 for b in boxes)) + max((b.y1 for b in boxes))) / 2
    for a in center_texts:
        x, y = a.get_position()
        a.set_position((x + dx, y + dy))
    label_map = m[['marker_id', 'full_label', 'display_label']]
    geometry, collisions = legacy.radial_labels(ax, fig, cells, m, label_map, 'B100')
    fx = geometry[[f'corner{i}_fx' for i in range(1, 5)]].to_numpy()
    fy = geometry[[f'corner{i}_fy' for i in range(1, 5)]].to_numpy()
    label_left = float((fx.min() - 0.5) * (2 * limit) / 0.95)
    label_top = float((fy.max() - 0.5) * (2 * limit) / 0.95)
    frame = dict(left=-1.358555753374844, top=1.2137663378025365, right=0.0, bottom=0.0, center_hole_radius=float(cells.r_inner.min()), notch_radius=float(cells.r_inner.min()), alignment_source='Frozen accepted fullnames_clean_20261105_v6 frame; not current label extrema', radial_label_left=label_left, radial_label_top=label_top)
    info, card_geometry = info_card(ax, counts, legacy.DOMAINS, scale, frame)
    links = []
    for i, ep in enumerate(ENDPOINTS):
        radius = 0.922 - i * 0.105 - 0.095 / 2
        ax.plot([-0.041, -0.012], [radius, radius], color='#607A88', lw=0.65, zorder=7)
        links.append(dict(endpoint_id=ep, ring_order_outer_to_inner=i + 1, bar_y=radius, ring_y=radius))
    fig.canvas.draw()
    renderer = fig.canvas.get_renderer()
    bounds = []
    for a in fig.findobj(Text):
        if a.get_visible() and a.get_text():
            b = a.get_window_extent(renderer).transformed(fig.transFigure.inverted())
            assert b.x0 >= 0 and b.y0 >= 0 and (b.x1 <= 1) and (b.y1 <= 1), (a.get_text(), b)
            expected_color = '#ffffff' if (a.get_gid() or '').startswith(('bar-count-white-', 'domain-curved-')) else '#000000'
            assert to_hex(a.get_color()) == expected_color, (a.get_text(), a.get_color())
            bounds.append(dict(text=a.get_text(), left=b.x0, bottom=b.y0, right=b.x1, top=b.y1, color=to_hex(a.get_color()), font_size_pt=a.get_fontsize(), weight=a.get_fontweight()))
    info_checks = []
    for a in info:
        b = a.get_window_extent(renderer).transformed(ax.transData.inverted())
        label = getattr(a, 'get_text', lambda: 'shape')()
        assert b.x0 >= frame['left'] and b.x1 <= frame['right'] and (b.y0 >= frame['bottom']) and (b.y1 <= frame['top']), (label, b)
        r_min = math.hypot(b.x1, b.y0)
        assert r_min > frame['notch_radius'], ('Information overlaps central-circle cutout', label, b, r_min)
        info_checks.append(dict(label=label, left=b.x0, right=b.x1, bottom=b.y0, top=b.y1, minimum_radius=r_min, outside_center_hole=True, inside_outward_frame=True))
    assert frame['notch_radius'] == frame['center_hole_radius']
    center_records = []
    for a in center_texts:
        b = a.get_window_extent(renderer).transformed(ax.transData.inverted())
        max_radius = max((math.hypot(x, y) for x in [b.x0, b.x1] for y in [b.y0, b.y1]))
        assert max_radius < frame['center_hole_radius'] - 0.01, (a.get_text(), b)
        center_records.append(dict(text=a.get_text(), font_size=a.get_fontsize(), weight=a.get_fontweight(), x=a.get_position()[0], y=a.get_position()[1], left=b.x0, right=b.x1, bottom=b.y0, top=b.y1, maximum_corner_radius=max_radius))
    center_x = (min((r['left'] for r in center_records)) + max((r['right'] for r in center_records))) / 2
    center_y = (min((r['bottom'] for r in center_records)) + max((r['top'] for r in center_records))) / 2
    assert abs(center_x) < 1e-12 and abs(center_y) < 1e-12, (center_x, center_y)
    card_geometry.update(inner_circle_center=[0, 0], notch_center=[0, 0], notch_exactly_matches_heatmap_inner_radius=True, previous_clearance_removed=0.01, central_text_block_center=[center_x, center_y], central_text_sizes=[31.2, 29.9, 16.9], five_disease_names_bold=False)
    jwrite(out / 'data/OUTWARD_CARD_GEOMETRY.json', card_geometry)
    for name, frame in [('B100_CIRCOS_605_CELLS', cells), ('STAR_QC_429', pd.DataFrame(starrows)), ('DOMAIN_GEOMETRY', pd.DataFrame(domains)), ('CURVED_DOMAIN_GLYPHS', pd.DataFrame(glyphs)), ('CURVED_DOMAIN_LABELS', pd.DataFrame(words)), ('RADIAL_LABEL_GEOMETRY', geometry), ('BAR_RING_ALIGNMENT', pd.DataFrame(links)), ('TEXT_PAGE_BOUNDS', pd.DataFrame(bounds)), ('INFORMATION_CONTENT_BOUNDS', pd.DataFrame(info_checks)), ('CENTER_TEXT_GEOMETRY', pd.DataFrame(center_records))]:
        frame.to_csv(out / 'data' / (name + '.csv'), index=False)
    crop = Bbox.from_extents(1.3260323596861439, 1.4415376014360621, 15.653637506740257, 15.247114802307188)
    for ext in ['pdf', 'svg', 'png']:
        fig.savefig(out / 'figures' / ('B100_all121_five_disease_IE_circos.' + ext), dpi=dpi, facecolor='white', bbox_inches=crop)
    fig.savefig(out / 'figures/B100_all121_five_disease_IE_circos_preview.png', dpi=150, facecolor='white', bbox_inches=crop)
    plt.close(fig)
    with pymupdf.open(out / 'figures/B100_all121_five_disease_IE_circos.pdf') as pdf:
        assert len(pdf) == 1 and (not pdf[0].get_images(full=True))
        fonts = pdf[0].get_fonts(full=True)
        assert fonts and all(('Arial' in f[3] and pdf.extract_font(f[0])[3] for f in fonts)), fonts
        spans = [s for b in pdf[0].get_text('dict')['blocks'] for line in b.get('lines', []) for s in line.get('spans', [])]
        assert all((pdf[0].rect.contains(pymupdf.Rect(s['bbox'])) for s in spans))
        pdf[0].get_pixmap(matrix=pymupdf.Matrix(1.5, 1.5), alpha=False).save(out / 'qa/ACTUAL_PDF_RENDER.png')
    return dict(indicators=121, diseases=5, cells=605, finite_IE=550, NA=55, stars=429, black_stars=int((cells.supported & cells.star_color.eq('#000000')).sum()), white_stars=int((cells.supported & cells.star_color.eq('#FFFFFF')).sum()), label_count=len(geometry), radial_label_collisions=collisions, domains=11, curved_labels=int(sum((x['drawn'] for x in words))), curved_glyph_bottoms_face_outward=True, bars_align_with_rings=True, card_left_and_top_align_with_radial_text=False, card_frame_frozen_from_accepted_v6=True, export_crop_frozen_from_accepted_v6=True, center_circle_upper_left_quadrant_unobstructed=True, information_content_outside_center_hole=True, circular_corner_cutout_not_white_overlay=True, notch_radius_equals_heatmap_inner_radius=True, central_text_block_centered=True, central_text_enlarged=True, five_disease_names_bold=False, card_geometry=card_geometry, no_text_clipping=True, all_stars_fit_cells=True, Arial_embedded=True, vector_PDF_SVG=True, figure_canvas_inches=[inches, inches], figure_inches=[crop.width, crop.height], PNG_dpi=dpi, all_nonbar_nondomain_text_pure_black=True, bar_count_text_white=True, domain_arc_titles_white=True, outer_label_leaders_drawn=0, bar_ring_links_retained=5, info_card_font_multiplier=1.3, disease_label_font_size_pt=18.5, disease_names_full_English=True, legend_y=0.98, negative_positive_font_size_pt=15.0, nonheading_info_text_weight='normal', radial_font_size_pt=9.6, radial_weight='normal', radial_label_radius=1.015, outer_label_positions_frozen_from_v7=True, only_star_diameter_changed_vs_v8=True, star_diameter_multiplier_vs_v8=1.4, star_area_multiplier_vs_v8=1.96, tightly_cropped_without_reducing_ring_scale=True, all_121_included_no_significance_filter=True, models_fitted=0, new_P_q_CI=0)

def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', type=Path, required=True)
    ap.add_argument('--style-root', type=Path, required=True)
    ap.add_argument('--social-map', type=Path, required=True)
    ap.add_argument('--out', type=Path, required=True)
    ap.add_argument('--dpi', type=int, default=600)
    args = ap.parse_args()
    source = args.input.resolve(strict=True)
    style_root = args.style_root.resolve(strict=True)
    social = args.social_map.resolve(strict=True)
    out = args.out.resolve()
    code = Path(__file__).resolve().parent
    print(f'INPUT={source}\nSTYLE_ROOT={style_root}\nSCRIPT_ROOT={code}\nHELPER_DIR={HELPER_DIR}\nOUTPUT={out}', flush=True)
    assert not out.exists() and args.dpi >= 300
    assert sha(source) == PINNED_SOURCE, 'Wrong B100 table: latest completion 550 CI required'
    manifest = json.loads((source.parent / 'TABLE_MANIFEST.json').read_text())
    assert manifest[source.name]['sha256'] == sha(source)
    d = pd.read_csv(source, dtype={'endpoint_id': str})
    assert len(d) == 605 and d.indicator_id.nunique() == 121
    assert set(d.endpoint_id) == set(ENDPOINTS)
    assert not d.duplicated(['indicator_id', 'endpoint_id']).any()
    assert d.horizon_years.eq(11).all() and d.bootstrap_planned.eq(100).all()
    assert d.model_followup.eq('FULL_FOLLOWUP').all()
    assert d.version.value_counts().to_dict() == {'B100_ORIGINAL_ESTIMATOR_WARNING_AWARE_RELEASE_V1': 565, 'B100_MISSING_INFERENCE_COMPLETION_20261030_V1_ORIGINAL_ESTIMATOR': 40}
    assert d.source_method_hash.eq(birth_private_setting('PRIVATE_7AA11EFF31702D72C96A44823457AB674B93D8CC8A1AD7965B72932F60DA')).all()
    assert np.isfinite(d.IE).sum() == np.isfinite(d.IE_CI_low).sum() == np.isfinite(d[Q_COLUMN]).sum() == 550
    valid = np.isfinite(d.IE)
    assert np.allclose(d.loc[valid, 'TE'], d.loc[valid, 'DE'] + d.loc[valid, 'IE'], atol=1e-10)
    assert np.allclose(d.loc[valid, 'IE_per_1000'], d.loc[valid, 'IE'] * 1000, atol=1e-10)
    assert d.loc[~valid, 'result_status'].eq('NOT_ESTIMABLE').all()
    m, style = mapping(d, style_root, social)
    d['source_csv_row'] = np.arange(1, len(d) + 1)
    d['marker_id'] = d.indicator_id
    d['display_domain_id'] = d.marker_id.map(m.set_index('marker_id').display_domain_id)
    d['marker_order'] = d.marker_id.map(m.set_index('marker_id').marker_order)
    d['ring_order'] = d.endpoint_id.map({v: i + 1 for i, v in enumerate(ENDPOINTS)})
    d, scale = colorize(d)
    counts = []
    for ep, (neg, pos) in zip(ENDPOINTS, EXPECTED_COUNTS):
        z = d[d.endpoint_id.eq(ep)]
        n = int((z.supported & z.IE.lt(0)).sum())
        p = int((z.supported & z.IE.gt(0)).sum())
        assert (n, p) == (neg, pos)
        counts.append(dict(endpoint_id=ep, negative=n, positive=p, total=n + p, estimated=int(z.finite_ie.sum()), missing=int((~z.finite_ie).sum())))
    counts = pd.DataFrame(counts)
    for folder in ['figures', 'data', 'source', 'code', 'qa', 'reports']:
        (out / folder).mkdir(parents=True)
    sources = [source, source.parent / 'TABLE_MANIFEST.json', source.parent / 'COMPLETION_REPORT_ZH.md', source.parent / 'COMPLETION_QA.json', social, style_root / 'DOMAIN_STYLE_UNIFIED.csv', style_root / 'MARKER_DISPLAY_DOMAIN_UNIFIED.csv', style_root / 'RADIAL_LABEL_MAPPING.csv']
    before = {str(p): sha(p) for p in sources}
    for p in sources:
        shutil.copy2(p, out / 'source' / p.name)
    for p in code.iterdir():
        if p.is_file():
            shutil.copy2(p, out / 'code' / p.name)
    counts.to_csv(out / 'data/DISEASE_SUPPORTED_COUNTS.csv', index=False)
    m.to_csv(out / 'data/INDICATOR_DOMAIN_LABEL_MAPPING.csv', index=False)
    d.groupby(['endpoint_id', 'result_status'], dropna=False).agg(cells=('indicator_id', 'size'), finite=('finite_ie', 'sum'), stars=('supported', 'sum')).reset_index().to_csv(out / 'data/STATUS_COUNTS.csv', index=False)
    jwrite(out / 'COLOR_SCALE.json', scale)
    plt.rcParams.update({'font.family': 'Arial', 'font.size': 8, 'svg.fonttype': 'none', 'pdf.fonttype': 42, 'ps.fonttype': 42, 'axes.unicode_minus': False})
    qc = render(d.sort_values(['marker_order', 'ring_order']), m, style, counts, scale, out, args.dpi)
    for p in sources:
        assert sha(p) == before[str(p)]
    plotted = pd.read_csv(out / 'data/B100_CIRCOS_605_CELLS.csv', dtype={'endpoint_id': str}).set_index(['indicator_id', 'endpoint_id'])
    original = pd.read_csv(source, dtype={'endpoint_id': str}).set_index(['indicator_id', 'endpoint_id'])
    pd.testing.assert_frame_equal(plotted.loc[original.index, original.columns], original, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-14)
    qc.update(status='NUMERIC_GEOMETRY_PASS_PENDING_VISUAL_REVIEW', latest_source_sha256=sha(source), source_statistical_columns_unchanged=True, source_warning_and_missingness_preserved=True, created_at=datetime.now(ZoneInfo('Asia/Shanghai')).isoformat())
    jwrite(out / 'QC.json', qc)
    jwrite(out / 'SOURCE_MANIFEST.json', dict(input=str(source), input_sha256=sha(source), source_files=before, source_unchanged=True, old_source_not_used='Earlier release had 510 CIs; completion table has 550 CIs', reused_functions=['legacy polar/geometry', 'radial labels and polygon collision checks', 'curved domain glyphs'], legacy_renderer=str(style_root / 'render_circos_bottom_out.py'), legacy_renderer_sha256=sha(style_root / 'render_circos_bottom_out.py'), legacy_changes='Only radial font size and label-count assertion adapted; legacy data loader and statistical flags not used', plotting_estimand='Saved per-indicator IE, not sum across indicators, TE, a*b, HR or mediation proportion', fixed_multiplicity_family=605, q_column=Q_COLUMN, models_fitted=0, new_P_q_CI=0))
    (out / 'reports/CAPTION_EN.md').write_text('Exploratory indirect-effect atlas of relative birth order (0 to 1) among participants from sibships of 2–8 children. All 121 registered B100 indicators and five disease outcomes are shown (605 combinations). Rings run from outer to inner: diabetes mellitus, type 2 diabetes, ischaemic heart disease, angina pectoris, and coronary atherosclerosis. Each disease bar is aligned with its corresponding ring at the upper seam. Colour encodes the saved 11-year indirect risk difference per 1,000 participants (IE); the symmetric colour limits use the 98th percentile of absolute finite IE values for display only. All 550 finite estimates are coloured regardless of significance; 55 unavailable estimates are grey, not zero. Stars denote saved exploratory BH q<0.05 across the fixed 605-combination family; black or white is chosen solely for contrast. Blue/red stacked counts denote negative/positive starred combinations. The latest completion release includes intervals for 40 previously point-only combinations. Partial bootstrap intervals and P values remain conditional on successful replicates; model-warning limitations are preserved. These are per-indicator exploratory effects, not additive contributions or independent mechanisms. No mediation models, bootstrap replicates, P values or q values were recomputed for this figure.\n', encoding='utf-8')
    (out / 'reports/README_ZH.md').write_text(f'# 最新B100五病全指标间接效应环图\n\n使用2026-10-30补齐版B100_121_COMPLETE_RESULTS.csv（SHA256见清单），不是较早只有510格区间的release。121指标×5病=605格，550格有IE与CI，55格不可估计原样灰色保留；全部指标均纳入，不按显著性挑子集。TE为总效应，IE为间接效应；本图按中介结果总览用途绘制IE，不把TE、旧a×b或PM混用，亦不相加121项IE。\n\n五病外至内为DM、T2D、IHD、Angina、冠脉粥样硬化。左上方形柱图逐行对齐对应疾病环；颜色使用每千人11年风险差IE，蓝负红正，大小由色阶表达。星号只用原表q_BH_global605_exploratory<0.05，固定M=605，不采用更小的病内121家族，不沿用旧291星。429星分别94、94、78、81、82；负/正计数见CSV。色标以全部有限IE绝对值98百分位对称截色，仅影响颜色，不改数值；浅色格黑星、深色格白星。\n\n分类沿用已统一的森林图/Circos领域字典，分类表和名称/缩写逐项可追溯。本批B100不含CMR或351，不补造这些结局/指标。NMR高度相关，不能把429格当作429个独立机制。部分B100格成功数少于100，新增蛋白等仍有警告接受策略；CI与P保留原条件性、探索性限制，星号不代表无警告、完整B100或因果识别已获证明。\n\n本轮只读聚合数值，未读取个体数据，未重拟合、重抽样或重算P/q/CI；旧结果及旧图不覆盖。输出可编辑PDF/SVG、{args.dpi}dpi PNG、每格源数据、映射、几何验收与英文图注。\n', encoding='utf-8')
    (out / 'RUN_COMMAND.txt').write_text(' '.join([sys.executable, str(code / 'build.py'), '--input', str(source), '--style-root', str(style_root), '--social-map', str(social), '--out', str(out), '--dpi', str(args.dpi)]) + '\n')
    print(json.dumps({'out': str(out), 'QC': qc}, ensure_ascii=False), flush=True)

def run_in_private_context(context):
    if not isinstance(context, dict):
        raise TypeError('Explicit private stage namespace required')
    exec(compile('"""Outer labels close to bands; every other position/crop frozen to accepted v6."""\nimport argparse\nimport hashlib\nimport json\nimport math\nimport os\nimport shutil\nimport sys\nfrom datetime import datetime\nfrom pathlib import Path\nfrom zoneinfo import ZoneInfo\nimport numpy as np\nimport pandas as pd\nimport matplotlib\nmatplotlib.use(\'Agg\')\nimport matplotlib.pyplot as plt\nfrom matplotlib.colors import to_hex, to_rgb\nfrom matplotlib.patches import Wedge, Rectangle, PathPatch, Polygon\nfrom matplotlib.path import Path as MPath\nfrom matplotlib.text import Text\nfrom matplotlib.transforms import Bbox\nimport pymupdf\nHELPER_DIR = Path(os.environ.get(\'BIRTH_CIRCOS_HELPER_DIR\', str(Path(__file__).resolve().parent))).resolve(strict=True)\nsys.path.insert(0, str(HELPER_DIR))\nimport legacy_circos_geometry as legacy\nENDPOINTS = [\'250\', \'250.2\', \'411\', \'411.3\', \'411.4\']\nSHORT = [\'DM\', \'T2D\', \'IHD\', \'Angina\', \'Cor. atheroscl.\']\nEXPECTED_COUNTS = [(77, 17), (77, 17), (55, 23), (53, 28), (49, 33)]\nQ_COLUMN = \'q_BH_global605_exploratory\'\nPINNED_SOURCE = birth_private_setting(\'PRIVATE_5990C4F5E198AC8DF73A29ABE9CB6E969FC3708409EDC2E2BB1C26C66FB0\')\nNEW_NAMES = {\'alcohol\': \'Alcohol status\', \'body_size10\': \'Childhood body size\', \'breastfed\': \'Breastfed\', \'chronotype\': \'Chronotype\', \'clinical__ldl\': \'LDL-C (clinical)\', \'clinical__tg\': \'Triglycerides\', \'commute_car\': \'Car commuting\', \'commute_public_transport\': \'Public transport\', \'commute_walking\': \'Walking to work\', \'confiding\': \'Confiding frequency\', \'employment_employed\': \'Paid/self-employment\', \'employment_student\': \'Student status\', \'employment_voluntary\': \'Voluntary work\', \'height10\': \'Childhood height\', \'housing_tenure\': \'Housing tenure\', \'income\': \'Household income\', \'insomnia\': \'Insomnia\', \'living_alone\': \'Living alone\', \'loneliness\': \'Loneliness\', \'manual_work\': \'Manual work\', \'maternal_smoking\': \'Maternal smoking\', \'nap\': \'Daytime naps\', \'oily_fish\': \'Oily fish\', \'participation_adult_education\': \'Adult-education group\', \'participation_other\': \'Other social group\', \'participation_religious\': \'Religious group\', \'qualification\': \'Education\', \'smoking\': \'Smoking status\', \'vehicles\': \'Vehicles\', \'visits\': \'Visits with friends/family\'}\nDOMAIN_MAP = {\'Early life\': \'EARLY_LIFE\', \'Education & SES\': \'EDUCATION_SOCIOECONOMIC\', \'Employment & occupation\': \'OCCUPATION\', \'Lifestyle (excluding smoking and alcohol)\': \'HEALTH_BEHAVIOUR\', \'Psychosocial & social connections\': \'PSYCHOSOCIAL_SOCIAL\', \'Geographical mobility\': \'GEOGRAPHY_ENVIRONMENT\'}\n\ndef sha(p):\n    return hashlib.sha256(p.read_bytes()).hexdigest()\n\ndef jwrite(p, d):\n    with p.open(\'x\') as f:\n        json.dump(d, f, ensure_ascii=False, indent=2, default=str, allow_nan=False)\n\ndef mapping(d, style_root, social_path):\n    style = pd.read_csv(style_root / \'DOMAIN_STYLE_UNIFIED.csv\').sort_values(\'display_order\')\n    old = pd.read_csv(style_root / \'MARKER_DISPLAY_DOMAIN_UNIFIED.csv\').set_index(\'marker_id\')\n    names = pd.read_csv(style_root / \'RADIAL_LABEL_MAPPING.csv\').set_index(\'marker_id\')\n    social = pd.read_csv(social_path).drop_duplicates(\'canonical_id\').set_index(\'canonical_id\')\n    records = []\n    for i, row in enumerate(d.drop_duplicates(\'indicator_id\').itertuples()):\n        ident = row.indicator_id\n        if ident in old.index:\n            domain = old.loc[ident, \'display_domain_id\']\n            basis = \'Reused approved forest/Circos display map\'\n        elif ident in social.index:\n            domain = DOMAIN_MAP[social.loc[ident, \'domain\']]\n            basis = \'Existing nontraditional A-forest domain, not regrouped by effect\'\n        elif ident.startswith(\'clinical__\'):\n            domain = \'CLINICAL_BIOCHEMISTRY\'\n            basis = \'Existing clinical-assay domain\'\n        elif ident in {\'smoking\', \'alcohol\'}:\n            domain = \'HEALTH_BEHAVIOUR\'\n            basis = \'Existing lifestyle domain\'\n        else:\n            raise ValueError(f\'Unmapped indicator: {ident}\')\n        label = NEW_NAMES.get(ident)\n        if label is None:\n            assert ident in names.index, ident\n            label = str(names.loc[ident, \'display_label\'])\n        if ident == \'clinical__hdl\':\n            label = \'HDL-C (clinical)\'\n        if ident == \'clinical__apoa\':\n            label = \'ApoA1 (clinical)\'\n        if ident == \'clinical__tc\':\n            label = \'Total cholesterol (clinical)\'\n        if ident == \'mobility\':\n            label = \'Birth–home distance\'\n        if ident == \'nmr__nmr_23404\':\n            label = \'LDL-C (NMR)\'\n        if ident == \'nmr__nmr_23406\':\n            label = \'HDL-C (NMR)\'\n        if ident in social.index:\n            order = float(social.loc[ident, \'marker_order\'])\n        elif ident in old.index:\n            order = float(old.loc[ident, \'source_marker_order\'])\n        else:\n            order = 1000 + i\n        records.append(dict(marker_id=ident, source_domain=row.domain, display_domain_id=domain, marker_order=order, source_indicator_name=row.indicator_name, full_label=str(social.loc[ident, \'label\']) if ident in social.index else str(names.loc[ident, \'full_label\']) if ident in names.index else label, display_label=label, domain_mapping_basis=basis, label_mapping_basis=\'Existing radial label or audited construct translation\'))\n    m = pd.DataFrame(records)\n    meta = style.set_index(\'display_domain_id\')\n    m[\'domain_order\'] = m.display_domain_id.map(meta.display_order)\n    m[\'display_domain\'] = m.display_domain_id.map(meta.display_label)\n    m[\'display_colour\'] = m.display_domain_id.map(meta.colour)\n    assert m.domain_order.notna().all() and m.marker_id.is_unique and (len(m) == 121)\n    assert m.display_label.is_unique, m.loc[m.display_label.duplicated(False)].to_dict(\'records\')\n    m = m.sort_values([\'domain_order\', \'marker_order\', \'marker_id\'], kind=\'stable\').reset_index(drop=True)\n    m[\'marker_order\'] = np.arange(1, len(m) + 1)\n    m[\'within_domain_order\'] = m.groupby(\'display_domain_id\', sort=False).cumcount() + 1\n    return (m, style)\n\ndef colorize(d):\n    d = d.copy()\n    d[\'finite_ie\'] = np.isfinite(d.IE_per_1000)\n    assert d.finite_ie.sum() == 550\n    d[\'supported\'] = d[Q_COLUMN].lt(0.05)\n    assert d.supported.sum() == 429 and (not (d.supported & ~d.finite_ie).any())\n    assert (d.loc[d.supported, \'IE_CI_low\'] * d.loc[d.supported, \'IE_CI_high\']).gt(0).all()\n    c = float(np.quantile(np.abs(d.loc[d.finite_ie, \'IE_per_1000\']), 0.98, method=\'linear\'))\n    d[\'display_ie_clipped\'] = d.IE_per_1000.clip(-c, c)\n    d[\'display_color_clipped\'] = d.finite_ie & d.IE_per_1000.abs().gt(c)\n    d[\'fill\'] = [to_hex(legacy.EFFECT_CMAP((x + c) / (2 * c))) if finite else legacy.NO_ESTIMATE_COLOUR for x, finite in zip(d.display_ie_clipped, d.finite_ie)]\n    rgb = np.array([to_rgb(x) for x in d.fill])\n    rgb_lin = np.where(rgb <= 0.04045, rgb / 12.92, ((rgb + 0.055) / 1.055) ** 2.4)\n    lum = rgb_lin @ np.array([0.2126, 0.7152, 0.0722])\n    d[\'relative_luminance\'] = lum\n    d[\'star_color\'] = np.where(d.supported, np.where(1.05 / (lum + 0.05) > (lum + 0.05) / 0.05, \'#FFFFFF\', \'#000000\'), \'\')\n    d[\'q_column_used\'] = Q_COLUMN\n    d[\'q_family_M\'] = 605\n    d[\'conditional_or_warned\'] = d.result_status.isin([\'PARTIAL_BOOTSTRAP_INFERENCE\', \'WARNING_AWARE_PARTIAL_BOOTSTRAP_INFERENCE\', \'WARNING_AWARE_BOOTSTRAP_INFERENCE\'])\n    return (d, dict(c=c, vmin=-c, vmax=c, center=0, finite_values=550, finite_min=float(d.IE_per_1000.min()), finite_max=float(d.IE_per_1000.max()), unit=\'11-year indirect risk difference per 1,000 participants\', percentile=98, percentile_method=\'Type 7; display only\', clipping_changes_source_values=False, missing_color=legacy.NO_ESTIMATE_COLOUR, palette=legacy.EFFECT_ANCHORS, star_rule=Q_COLUMN + \' < 0.05; fixed M=605\', star_contrast=\'Higher black/white relative-luminance contrast against actual fill\'))\n\ndef info_card(ax, counts, domains, scale, frame):\n    """Outward information region with a true circular cut-out, not a mask."""\n    left, top, right, bottom = (frame[\'left\'], frame[\'top\'], frame[\'right\'], frame[\'bottom\'])\n    radius = frame[\'notch_radius\']\n    theta0 = math.degrees(math.acos(right / radius))\n    theta1 = 180 - math.degrees(math.asin(bottom / radius))\n    arc = MPath.arc(theta0, theta1, n=32)\n    av = arc.vertices * radius\n    vertices = np.vstack([[(left, top), (right, top), av[0]], av[1:], [(left, bottom), (left, top)]])\n    codes = [MPath.MOVETO, MPath.LINETO, MPath.LINETO] + list(arc.codes[1:]) + [MPath.LINETO, MPath.CLOSEPOLY]\n    box = PathPatch(MPath(vertices, codes), fc=\'#FBFCFE\', ec=\'#B9CEDB\', lw=1.1, zorder=1)\n    box.set_gid(\'outward-information-card-with-circular-notch\')\n    ax.add_patch(box)\n    artists = []\n\n    def txt(x, y, text, size=8, **kw):\n        weight = kw.pop(\'weight\', \'normal\')\n        a = ax.text(x, y, text, fontsize=size * 1.3, color=\'#000000\', weight=weight, va=\'center\', zorder=7, **kw)\n        artists.append(a)\n        return a\n\n    def rect(x, y, w, h, color):\n        a = Rectangle((x, y), w, h, fc=color, ec=\'none\', zorder=4)\n        ax.add_patch(a)\n        artists.append(a)\n        return a\n    text_left = left + 0.048\n    title_y = top - 0.057\n    legend_y = 0.98\n    txt(text_left, title_y, \'Indirect effects by disease\', 15, weight=\'bold\')\n    for x, label, color in [(text_left, \'Negative\', \'#2C6E96\'), (text_left + 0.41, \'Positive\', \'#CE6C63\')]:\n        rect(x, legend_y - 0.015, 0.029, 0.029, color)\n        txt(x + 0.047, legend_y, label, 15 / 1.3)\n    txt(-0.065, legend_y, \'Total\', 9.2, ha=\'center\', weight=\'bold\')\n    x0, span, height = (-0.87, 0.735, 0.053)\n    names = [\'Diabetes mellitus\', \'Type 2 diabetes\', \'Ischemic heart\\ndisease\', \'Angina pectoris\', \'Coronary\\natherosclerosis\']\n    for t in range(0, 101, 20):\n        x = x0 + span * t / 100\n        ax.plot([x, x], [0.402, 0.908], color=\'#E5EBEF\', lw=0.5, zorder=2)\n        txt(x, 0.414, str(t), 8, ha=\'center\')\n    for i, row in enumerate(counts.itertuples()):\n        y = 0.922 - i * 0.105 - 0.095 / 2\n        txt(text_left, y, names[i], 18.5 / 1.3, linespacing=1.0)\n        for start, n, color in [(0, row.negative, \'#2C6E96\'), (row.negative, row.positive, \'#CE6C63\')]:\n            rect(x0 + span * start / 100, y - height / 2, span * n / 100, height, color)\n            if n:\n                a = txt(x0 + span * (start + n / 2) / 100, y, str(n), 9, ha=\'center\')\n                a.set_color(\'#FFFFFF\')\n                a.set_gid(f\'bar-count-white-{i}-{start}\')\n        txt(-0.065, y, str(row.total), 11.2, ha=\'center\')\n    ax.plot([x0, x0 + span], [0.402, 0.402], color=legacy.INK, lw=0.6, zorder=4)\n    txt(x0 + span / 2, 0.35, \'Number with BH q < 0.05\', 8.8, ha=\'center\')\n    ax.plot([left + 0.025, -0.29], [0.325, 0.325], color=\'#DAE4EB\', lw=0.7)\n    txt(text_left, 0.294, \'IE per 1,000\', 10, weight=\'bold\')\n    txt(text_left, 0.255, \'Risk difference\', 8)\n    bx, by, bw, bh = (text_left + 0.005, 0.208, 0.29, 0.027)\n    for k in range(256):\n        rect(bx + bw * k / 256, by, bw / 256 + 2e-05, bh, to_hex(legacy.EFFECT_CMAP((k + 0.5) / 256)))\n    for xs, colour in [([bx, bx - 0.014, bx], legacy.EFFECT_ANCHORS[0]), ([bx + bw, bx + bw + 0.014, bx + bw], legacy.EFFECT_ANCHORS[-1])]:\n        ax.add_patch(Polygon(list(zip(xs, [by, by + bh / 2, by + bh])), fc=colour, ec=\'none\', zorder=4))\n    for x, value in [(bx, -scale[\'c\']), (bx + bw / 2, 0), (bx + bw, scale[\'c\'])]:\n        txt(x, 0.182, \'0\' if value == 0 else f\'{value:+.2f}\', 8, ha=\'center\')\n    txt(text_left, 0.139, \'*: BH q < 0.05 (M = 605)\', 8)\n    rect(text_left, 0.082, 0.024, 0.023, legacy.NO_ESTIMATE_COLOUR)\n    txt(text_left + 0.041, 0.093, \'No estimate\', 8)\n    txt(text_left, 0.048, \'550 estimates; 55 unavailable\', 8)\n    txt(-0.955, 0.294, \'Domains\', 10, weight=\'bold\')\n    nrows = math.ceil(len(domains) / 2)\n    for k, (_, label, _, color) in enumerate(domains):\n        col, j = divmod(k, nrows)\n        x = -0.955 + col * 0.295\n        y = 0.253 - j * 0.044\n        rect(x, y - 0.008, 0.015, 0.015, color)\n        txt(x + 0.025, y, label, 7.2)\n    return (artists, dict(**frame, label_left=text_left, bar_start=x0, bar_width=span, previous_bar_width=0.535, bar_width_multiplier=span / 0.535, domain_legend_columns=[-0.955, -0.66], no_white_mask_used=True, card_path_vertices=vertices.tolist(), card_path_codes=[int(x) for x in codes]))\n\ndef render(d, m, style, counts, scale, out, dpi):\n    present = set(m.display_domain_id)\n    legacy.DOMAINS = [(r.display_domain_id, r.display_label, r.arc_label, r.colour) for r in style.itertuples() if r.display_domain_id in present]\n    cells, domains, width = legacy.geometry(d, m)\n    assert len(cells) == 605 and len(domains) == 11\n    inches, limit = (16.0 * 1.55 / 1.43, 1.55)\n    fig, ax = plt.subplots(figsize=(inches, inches), dpi=120)\n    fig.subplots_adjust(left=0.025, right=0.975, bottom=0.025, top=0.975)\n    ax.set(xlim=(-limit, limit), ylim=(-limit, limit), aspect=\'equal\')\n    ax.axis(\'off\')\n    for k in range(5):\n        ax.add_patch(Wedge((0, 0), 0.922 - k * 0.105, 180, 450, width=0.095, fc=\'#F3F5F5\', ec=\'none\', zorder=1))\n    for i, row in enumerate(cells.itertuples()):\n        a = Wedge((0, 0), row.r_outer, row.theta_start, row.theta_end, width=row.r_outer - row.r_inner, fc=row.fill, ec=\'white\', lw=0.2, zorder=2)\n        a.set_gid(f\'cell-B100-{i + 1:03d}\')\n        ax.add_patch(a)\n    ppu = inches * 0.95 / (2 * limit) * 72\n    diameter = min(5.6, 0.84 * cells.loc[cells.supported, \'r_inner\'].min() * math.radians(width) * ppu)\n    star_radius = diameter / 2 / ppu\n    starrows = []\n    for i, row in enumerate(cells[cells.supported].to_dict(\'records\')):\n        a = ax.scatter(row[\'x\'], row[\'y\'], marker=MPath.unit_regular_star(5), s=diameter ** 2, c=row[\'star_color\'], edgecolors=\'none\', zorder=5)\n        a.set_gid(f\'B100-q-star-{i + 1:03d}\')\n        assert row[\'r_inner\'] < math.hypot(row[\'x\'], row[\'y\']) - star_radius\n        assert math.hypot(row[\'x\'], row[\'y\']) + star_radius < row[\'r_outer\']\n        assert math.degrees(math.asin(star_radius / math.hypot(row[\'x\'], row[\'y\']))) < width / 2\n        starrows.append({**row, \'star_diameter_pt\': diameter, \'star_fits_cell\': True})\n    assert len(starrows) == 429\n    glyphs, words = ([], [])\n    for rec in domains:\n        spec = next((x for x in legacy.DOMAINS if x[0] == rec[\'domain_id\']))\n        a = Wedge((0, 0), 1, rec[\'theta_start\'], rec[\'theta_end\'], width=0.052, fc=spec[3], ec=\'white\', lw=0.5, zorder=3)\n        a.set_gid(\'domain-band-\' + spec[0])\n        ax.add_patch(a)\n        gs, word = legacy.curved_domain_label(ax, spec, rec, ppu, \'B100\')\n        glyphs.extend(gs)\n        words.append(word)\n    center_texts = []\n    for y, label, size, weight in [(0.125, \'BIRTH ORDER\', 31.2, \'bold\'), (0.035, \'INDIRECT EFFECTS\', 29.9, \'bold\'), (-0.07, \'Relative birth order (0 → 1)\', 16.9, \'normal\')]:\n        center_texts.append(ax.text(0, y, label, ha=\'center\', va=\'center\', fontsize=size, weight=weight, color=\'#000000\'))\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    boxes = [a.get_window_extent(renderer).transformed(ax.transData.inverted()) for a in center_texts]\n    dx = -(min((b.x0 for b in boxes)) + max((b.x1 for b in boxes))) / 2\n    dy = -(min((b.y0 for b in boxes)) + max((b.y1 for b in boxes))) / 2\n    for a in center_texts:\n        x, y = a.get_position()\n        a.set_position((x + dx, y + dy))\n    label_map = m[[\'marker_id\', \'full_label\', \'display_label\']]\n    geometry, collisions = legacy.radial_labels(ax, fig, cells, m, label_map, \'B100\')\n    fx = geometry[[f\'corner{i}_fx\' for i in range(1, 5)]].to_numpy()\n    fy = geometry[[f\'corner{i}_fy\' for i in range(1, 5)]].to_numpy()\n    label_left = float((fx.min() - 0.5) * (2 * limit) / 0.95)\n    label_top = float((fy.max() - 0.5) * (2 * limit) / 0.95)\n    frame = dict(left=-1.358555753374844, top=1.2137663378025365, right=0.0, bottom=0.0, center_hole_radius=float(cells.r_inner.min()), notch_radius=float(cells.r_inner.min()), alignment_source=\'Frozen accepted fullnames_clean_20261105_v6 frame; not current label extrema\', radial_label_left=label_left, radial_label_top=label_top)\n    info, card_geometry = info_card(ax, counts, legacy.DOMAINS, scale, frame)\n    links = []\n    for i, ep in enumerate(ENDPOINTS):\n        radius = 0.922 - i * 0.105 - 0.095 / 2\n        ax.plot([-0.041, -0.012], [radius, radius], color=\'#607A88\', lw=0.65, zorder=7)\n        links.append(dict(endpoint_id=ep, ring_order_outer_to_inner=i + 1, bar_y=radius, ring_y=radius))\n    fig.canvas.draw()\n    renderer = fig.canvas.get_renderer()\n    bounds = []\n    for a in fig.findobj(Text):\n        if a.get_visible() and a.get_text():\n            b = a.get_window_extent(renderer).transformed(fig.transFigure.inverted())\n            assert b.x0 >= 0 and b.y0 >= 0 and (b.x1 <= 1) and (b.y1 <= 1), (a.get_text(), b)\n            expected_color = \'#ffffff\' if (a.get_gid() or \'\').startswith((\'bar-count-white-\', \'domain-curved-\')) else \'#000000\'\n            assert to_hex(a.get_color()) == expected_color, (a.get_text(), a.get_color())\n            bounds.append(dict(text=a.get_text(), left=b.x0, bottom=b.y0, right=b.x1, top=b.y1, color=to_hex(a.get_color()), font_size_pt=a.get_fontsize(), weight=a.get_fontweight()))\n    info_checks = []\n    for a in info:\n        b = a.get_window_extent(renderer).transformed(ax.transData.inverted())\n        label = getattr(a, \'get_text\', lambda: \'shape\')()\n        assert b.x0 >= frame[\'left\'] and b.x1 <= frame[\'right\'] and (b.y0 >= frame[\'bottom\']) and (b.y1 <= frame[\'top\']), (label, b)\n        r_min = math.hypot(b.x1, b.y0)\n        assert r_min > frame[\'notch_radius\'], (\'Information overlaps central-circle cutout\', label, b, r_min)\n        info_checks.append(dict(label=label, left=b.x0, right=b.x1, bottom=b.y0, top=b.y1, minimum_radius=r_min, outside_center_hole=True, inside_outward_frame=True))\n    assert frame[\'notch_radius\'] == frame[\'center_hole_radius\']\n    center_records = []\n    for a in center_texts:\n        b = a.get_window_extent(renderer).transformed(ax.transData.inverted())\n        max_radius = max((math.hypot(x, y) for x in [b.x0, b.x1] for y in [b.y0, b.y1]))\n        assert max_radius < frame[\'center_hole_radius\'] - 0.01, (a.get_text(), b)\n        center_records.append(dict(text=a.get_text(), font_size=a.get_fontsize(), weight=a.get_fontweight(), x=a.get_position()[0], y=a.get_position()[1], left=b.x0, right=b.x1, bottom=b.y0, top=b.y1, maximum_corner_radius=max_radius))\n    center_x = (min((r[\'left\'] for r in center_records)) + max((r[\'right\'] for r in center_records))) / 2\n    center_y = (min((r[\'bottom\'] for r in center_records)) + max((r[\'top\'] for r in center_records))) / 2\n    assert abs(center_x) < 1e-12 and abs(center_y) < 1e-12, (center_x, center_y)\n    card_geometry.update(inner_circle_center=[0, 0], notch_center=[0, 0], notch_exactly_matches_heatmap_inner_radius=True, previous_clearance_removed=0.01, central_text_block_center=[center_x, center_y], central_text_sizes=[31.2, 29.9, 16.9], five_disease_names_bold=False)\n    jwrite(out / \'data/OUTWARD_CARD_GEOMETRY.json\', card_geometry)\n    for name, frame in [(\'B100_CIRCOS_605_CELLS\', cells), (\'STAR_QC_429\', pd.DataFrame(starrows)), (\'DOMAIN_GEOMETRY\', pd.DataFrame(domains)), (\'CURVED_DOMAIN_GLYPHS\', pd.DataFrame(glyphs)), (\'CURVED_DOMAIN_LABELS\', pd.DataFrame(words)), (\'RADIAL_LABEL_GEOMETRY\', geometry), (\'BAR_RING_ALIGNMENT\', pd.DataFrame(links)), (\'TEXT_PAGE_BOUNDS\', pd.DataFrame(bounds)), (\'INFORMATION_CONTENT_BOUNDS\', pd.DataFrame(info_checks)), (\'CENTER_TEXT_GEOMETRY\', pd.DataFrame(center_records))]:\n        frame.to_csv(out / \'data\' / (name + \'.csv\'), index=False)\n    crop = Bbox.from_extents(1.3260323596861439, 1.4415376014360621, 15.653637506740257, 15.247114802307188)\n    for ext in [\'pdf\', \'svg\', \'png\']:\n        fig.savefig(out / \'figures\' / (\'B100_all121_five_disease_IE_circos.\' + ext), dpi=dpi, facecolor=\'white\', bbox_inches=crop)\n    fig.savefig(out / \'figures/B100_all121_five_disease_IE_circos_preview.png\', dpi=150, facecolor=\'white\', bbox_inches=crop)\n    plt.close(fig)\n    with pymupdf.open(out / \'figures/B100_all121_five_disease_IE_circos.pdf\') as pdf:\n        assert len(pdf) == 1 and (not pdf[0].get_images(full=True))\n        fonts = pdf[0].get_fonts(full=True)\n        assert fonts and all((\'Arial\' in f[3] and pdf.extract_font(f[0])[3] for f in fonts)), fonts\n        spans = [s for b in pdf[0].get_text(\'dict\')[\'blocks\'] for line in b.get(\'lines\', []) for s in line.get(\'spans\', [])]\n        assert all((pdf[0].rect.contains(pymupdf.Rect(s[\'bbox\'])) for s in spans))\n        pdf[0].get_pixmap(matrix=pymupdf.Matrix(1.5, 1.5), alpha=False).save(out / \'qa/ACTUAL_PDF_RENDER.png\')\n    return dict(indicators=121, diseases=5, cells=605, finite_IE=550, NA=55, stars=429, black_stars=int((cells.supported & cells.star_color.eq(\'#000000\')).sum()), white_stars=int((cells.supported & cells.star_color.eq(\'#FFFFFF\')).sum()), label_count=len(geometry), radial_label_collisions=collisions, domains=11, curved_labels=int(sum((x[\'drawn\'] for x in words))), curved_glyph_bottoms_face_outward=True, bars_align_with_rings=True, card_left_and_top_align_with_radial_text=False, card_frame_frozen_from_accepted_v6=True, export_crop_frozen_from_accepted_v6=True, center_circle_upper_left_quadrant_unobstructed=True, information_content_outside_center_hole=True, circular_corner_cutout_not_white_overlay=True, notch_radius_equals_heatmap_inner_radius=True, central_text_block_centered=True, central_text_enlarged=True, five_disease_names_bold=False, card_geometry=card_geometry, no_text_clipping=True, all_stars_fit_cells=True, Arial_embedded=True, vector_PDF_SVG=True, figure_canvas_inches=[inches, inches], figure_inches=[crop.width, crop.height], PNG_dpi=dpi, all_nonbar_nondomain_text_pure_black=True, bar_count_text_white=True, domain_arc_titles_white=True, outer_label_leaders_drawn=0, bar_ring_links_retained=5, info_card_font_multiplier=1.3, disease_label_font_size_pt=18.5, disease_names_full_English=True, legend_y=0.98, negative_positive_font_size_pt=15.0, nonheading_info_text_weight=\'normal\', radial_font_size_pt=9.6, radial_weight=\'normal\', radial_label_radius=1.015, outer_label_positions_frozen_from_v7=True, only_star_diameter_changed_vs_v8=True, star_diameter_multiplier_vs_v8=1.4, star_area_multiplier_vs_v8=1.96, tightly_cropped_without_reducing_ring_scale=True, all_121_included_no_significance_filter=True, models_fitted=0, new_P_q_CI=0)\n\ndef main():\n    ap = argparse.ArgumentParser(description=__doc__)\n    ap.add_argument(\'--input\', type=Path, required=True)\n    ap.add_argument(\'--style-root\', type=Path, required=True)\n    ap.add_argument(\'--social-map\', type=Path, required=True)\n    ap.add_argument(\'--out\', type=Path, required=True)\n    ap.add_argument(\'--dpi\', type=int, default=600)\n    args = ap.parse_args()\n    source = args.input.resolve(strict=True)\n    style_root = args.style_root.resolve(strict=True)\n    social = args.social_map.resolve(strict=True)\n    out = args.out.resolve()\n    code = Path(__file__).resolve().parent\n    print(f\'INPUT={source}\\nSTYLE_ROOT={style_root}\\nSCRIPT_ROOT={code}\\nHELPER_DIR={HELPER_DIR}\\nOUTPUT={out}\', flush=True)\n    assert not out.exists() and args.dpi >= 300\n    assert sha(source) == PINNED_SOURCE, \'Wrong B100 table: latest completion 550 CI required\'\n    manifest = json.loads((source.parent / \'TABLE_MANIFEST.json\').read_text())\n    assert manifest[source.name][\'sha256\'] == sha(source)\n    d = pd.read_csv(source, dtype={\'endpoint_id\': str})\n    assert len(d) == 605 and d.indicator_id.nunique() == 121\n    assert set(d.endpoint_id) == set(ENDPOINTS)\n    assert not d.duplicated([\'indicator_id\', \'endpoint_id\']).any()\n    assert d.horizon_years.eq(11).all() and d.bootstrap_planned.eq(100).all()\n    assert d.model_followup.eq(\'FULL_FOLLOWUP\').all()\n    assert d.version.value_counts().to_dict() == {\'B100_ORIGINAL_ESTIMATOR_WARNING_AWARE_RELEASE_V1\': 565, \'B100_MISSING_INFERENCE_COMPLETION_20261030_V1_ORIGINAL_ESTIMATOR\': 40}\n    assert d.source_method_hash.eq(birth_private_setting(\'PRIVATE_7AA11EFF31702D72C96A44823457AB674B93D8CC8A1AD7965B72932F60DA\')).all()\n    assert np.isfinite(d.IE).sum() == np.isfinite(d.IE_CI_low).sum() == np.isfinite(d[Q_COLUMN]).sum() == 550\n    valid = np.isfinite(d.IE)\n    assert np.allclose(d.loc[valid, \'TE\'], d.loc[valid, \'DE\'] + d.loc[valid, \'IE\'], atol=1e-10)\n    assert np.allclose(d.loc[valid, \'IE_per_1000\'], d.loc[valid, \'IE\'] * 1000, atol=1e-10)\n    assert d.loc[~valid, \'result_status\'].eq(\'NOT_ESTIMABLE\').all()\n    m, style = mapping(d, style_root, social)\n    d[\'source_csv_row\'] = np.arange(1, len(d) + 1)\n    d[\'marker_id\'] = d.indicator_id\n    d[\'display_domain_id\'] = d.marker_id.map(m.set_index(\'marker_id\').display_domain_id)\n    d[\'marker_order\'] = d.marker_id.map(m.set_index(\'marker_id\').marker_order)\n    d[\'ring_order\'] = d.endpoint_id.map({v: i + 1 for i, v in enumerate(ENDPOINTS)})\n    d, scale = colorize(d)\n    counts = []\n    for ep, (neg, pos) in zip(ENDPOINTS, EXPECTED_COUNTS):\n        z = d[d.endpoint_id.eq(ep)]\n        n = int((z.supported & z.IE.lt(0)).sum())\n        p = int((z.supported & z.IE.gt(0)).sum())\n        assert (n, p) == (neg, pos)\n        counts.append(dict(endpoint_id=ep, negative=n, positive=p, total=n + p, estimated=int(z.finite_ie.sum()), missing=int((~z.finite_ie).sum())))\n    counts = pd.DataFrame(counts)\n    for folder in [\'figures\', \'data\', \'source\', \'code\', \'qa\', \'reports\']:\n        (out / folder).mkdir(parents=True)\n    sources = [source, source.parent / \'TABLE_MANIFEST.json\', source.parent / \'COMPLETION_REPORT_ZH.md\', source.parent / \'COMPLETION_QA.json\', social, style_root / \'DOMAIN_STYLE_UNIFIED.csv\', style_root / \'MARKER_DISPLAY_DOMAIN_UNIFIED.csv\', style_root / \'RADIAL_LABEL_MAPPING.csv\']\n    before = {str(p): sha(p) for p in sources}\n    for p in sources:\n        shutil.copy2(p, out / \'source\' / p.name)\n    for p in code.iterdir():\n        if p.is_file():\n            shutil.copy2(p, out / \'code\' / p.name)\n    counts.to_csv(out / \'data/DISEASE_SUPPORTED_COUNTS.csv\', index=False)\n    m.to_csv(out / \'data/INDICATOR_DOMAIN_LABEL_MAPPING.csv\', index=False)\n    d.groupby([\'endpoint_id\', \'result_status\'], dropna=False).agg(cells=(\'indicator_id\', \'size\'), finite=(\'finite_ie\', \'sum\'), stars=(\'supported\', \'sum\')).reset_index().to_csv(out / \'data/STATUS_COUNTS.csv\', index=False)\n    jwrite(out / \'COLOR_SCALE.json\', scale)\n    plt.rcParams.update({\'font.family\': \'Arial\', \'font.size\': 8, \'svg.fonttype\': \'none\', \'pdf.fonttype\': 42, \'ps.fonttype\': 42, \'axes.unicode_minus\': False})\n    qc = render(d.sort_values([\'marker_order\', \'ring_order\']), m, style, counts, scale, out, args.dpi)\n    for p in sources:\n        assert sha(p) == before[str(p)]\n    plotted = pd.read_csv(out / \'data/B100_CIRCOS_605_CELLS.csv\', dtype={\'endpoint_id\': str}).set_index([\'indicator_id\', \'endpoint_id\'])\n    original = pd.read_csv(source, dtype={\'endpoint_id\': str}).set_index([\'indicator_id\', \'endpoint_id\'])\n    pd.testing.assert_frame_equal(plotted.loc[original.index, original.columns], original, check_dtype=False, check_exact=False, rtol=1e-12, atol=1e-14)\n    qc.update(status=\'NUMERIC_GEOMETRY_PASS_PENDING_VISUAL_REVIEW\', latest_source_sha256=sha(source), source_statistical_columns_unchanged=True, source_warning_and_missingness_preserved=True, created_at=datetime.now(ZoneInfo(\'Asia/Shanghai\')).isoformat())\n    jwrite(out / \'QC.json\', qc)\n    jwrite(out / \'SOURCE_MANIFEST.json\', dict(input=str(source), input_sha256=sha(source), source_files=before, source_unchanged=True, old_source_not_used=\'Earlier release had 510 CIs; completion table has 550 CIs\', reused_functions=[\'legacy polar/geometry\', \'radial labels and polygon collision checks\', \'curved domain glyphs\'], legacy_renderer=str(style_root / \'render_circos_bottom_out.py\'), legacy_renderer_sha256=sha(style_root / \'render_circos_bottom_out.py\'), legacy_changes=\'Only radial font size and label-count assertion adapted; legacy data loader and statistical flags not used\', plotting_estimand=\'Saved per-indicator IE, not sum across indicators, TE, a*b, HR or mediation proportion\', fixed_multiplicity_family=605, q_column=Q_COLUMN, models_fitted=0, new_P_q_CI=0))\n    (out / \'reports/CAPTION_EN.md\').write_text(\'Exploratory indirect-effect atlas of relative birth order (0 to 1) among participants from sibships of 2–8 children. All 121 registered B100 indicators and five disease outcomes are shown (605 combinations). Rings run from outer to inner: diabetes mellitus, type 2 diabetes, ischaemic heart disease, angina pectoris, and coronary atherosclerosis. Each disease bar is aligned with its corresponding ring at the upper seam. Colour encodes the saved 11-year indirect risk difference per 1,000 participants (IE); the symmetric colour limits use the 98th percentile of absolute finite IE values for display only. All 550 finite estimates are coloured regardless of significance; 55 unavailable estimates are grey, not zero. Stars denote saved exploratory BH q<0.05 across the fixed 605-combination family; black or white is chosen solely for contrast. Blue/red stacked counts denote negative/positive starred combinations. The latest completion release includes intervals for 40 previously point-only combinations. Partial bootstrap intervals and P values remain conditional on successful replicates; model-warning limitations are preserved. These are per-indicator exploratory effects, not additive contributions or independent mechanisms. No mediation models, bootstrap replicates, P values or q values were recomputed for this figure.\\n\', encoding=\'utf-8\')\n    (out / \'reports/README_ZH.md\').write_text(f\'# 最新B100五病全指标间接效应环图\\n\\n使用2026-10-30补齐版B100_121_COMPLETE_RESULTS.csv（SHA256见清单），不是较早只有510格区间的release。121指标×5病=605格，550格有IE与CI，55格不可估计原样灰色保留；全部指标均纳入，不按显著性挑子集。TE为总效应，IE为间接效应；本图按中介结果总览用途绘制IE，不把TE、旧a×b或PM混用，亦不相加121项IE。\\n\\n五病外至内为DM、T2D、IHD、Angina、冠脉粥样硬化。左上方形柱图逐行对齐对应疾病环；颜色使用每千人11年风险差IE，蓝负红正，大小由色阶表达。星号只用原表q_BH_global605_exploratory<0.05，固定M=605，不采用更小的病内121家族，不沿用旧291星。429星分别94、94、78、81、82；负/正计数见CSV。色标以全部有限IE绝对值98百分位对称截色，仅影响颜色，不改数值；浅色格黑星、深色格白星。\\n\\n分类沿用已统一的森林图/Circos领域字典，分类表和名称/缩写逐项可追溯。本批B100不含CMR或351，不补造这些结局/指标。NMR高度相关，不能把429格当作429个独立机制。部分B100格成功数少于100，新增蛋白等仍有警告接受策略；CI与P保留原条件性、探索性限制，星号不代表无警告、完整B100或因果识别已获证明。\\n\\n本轮只读聚合数值，未读取个体数据，未重拟合、重抽样或重算P/q/CI；旧结果及旧图不覆盖。输出可编辑PDF/SVG、{args.dpi}dpi PNG、每格源数据、映射、几何验收与英文图注。\\n\', encoding=\'utf-8\')\n    (out / \'RUN_COMMAND.txt\').write_text(\' \'.join([sys.executable, str(code / \'build.py\'), \'--input\', str(source), \'--style-root\', str(style_root), \'--social-map\', str(social), \'--out\', str(out), \'--dpi\', str(args.dpi)]) + \'\\n\')\n    print(json.dumps({\'out\': str(out), \'QC\': qc}, ensure_ascii=False), flush=True)\nif __name__ == \'__main__\':\n    main()', '<private-stage>', 'exec'), context)
