#!/usr/bin/env python3
"""Render core-four gene expression on supplied UMAP and t-SNE coordinates.

All panels use the same expression mapping and publication geometry.
No embedding, normalization, statistical testing or model fitting is performed.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import sys
sys.dont_write_bytecode = True


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as handle:
        for chunk in iter(lambda: handle.read(1048576), b''):
            h.update(chunk)
    return h.hexdigest()


def write_json(path, value):
    with path.open('x') as handle:
        json.dump(value, handle, indent=2, ensure_ascii=False)
        handle.write('\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for key in ['input-root', 'output-root']:
        default = os.getenv('SCRNA_SHARED_SCALE_' + key.replace('-', '_').upper())
        parser.add_argument('--' + key, default=default, required=default is None)
    args = parser.parse_args()
    inp, out = [Path(getattr(args, key)).expanduser().resolve()
                for key in ['input_root', 'output_root']]
    if out.exists():
        raise FileExistsError(f'Refusing existing output: {out}')
    if out.is_relative_to(inp):
        raise ValueError('Output root must be outside the input directory')
    dslist = ['GSE289422', 'GSE255566']
    genes = ['FCRL6', 'GZMA', 'GZMH', 'KLRD1']
    sources = [inp/f'results/{ds}/coordinates.csv.gz' for ds in dslist]
    for path in sources:
        if not path.is_file():
            raise FileNotFoundError(path)
    hashes = {str(path):sha(path) for path in sources}
    out.mkdir(parents=True)
    figures, metadata = out/'figures', out/'metadata'
    figures.mkdir(); metadata.mkdir()
    os.environ.setdefault('MPLCONFIGDIR', str(metadata/'matplotlib_cache'))
    import numpy as np
    import pandas as pd
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    from matplotlib.colors import LinearSegmentedColormap, Normalize
    from matplotlib.cm import ScalarMappable
    from matplotlib.font_manager import FontProperties, findfont

    font_family = 'DejaVu Sans'
    for candidate in ['Arial', 'Helvetica']:
        try:
            findfont(FontProperties(family=candidate), fallback_to_default=False)
            font_family = candidate
            break
        except ValueError:
            pass
    plt.rcParams.update({'font.family':'sans-serif',
        'font.sans-serif':[font_family,'Arial','Helvetica','DejaVu Sans'], 'font.size':7.3,
        'axes.labelsize':7, 'axes.titlesize':8.5, 'xtick.labelsize':7, 'ytick.labelsize':7,
        'pdf.fonttype':42, 'ps.fonttype':42, 'svg.fonttype':'none',
        'text.color':'#202020','axes.labelcolor':'#202020',
        'figure.facecolor':'white','savefig.facecolor':'white',
        'svg.hashsalt':'core4-ocean-gold'})
    dat = {ds:pd.read_csv(path) for ds,path in zip(dslist, sources[:2])}
    for ds,d in dat.items():
        assert d.cell_id.is_unique
        assert len(d) == (44000 if ds == 'GSE289422' else 6000)
        assert set(d.group) == {'Control', 'CAD' if ds == 'GSE289422' else 'T2D'}
        values = d[[g+'_logCP10k' for g in genes]].to_numpy()
        assert np.isfinite(values).all() and (values >= 0).all() and (values <= 5).all()

    # Identical numerical mapping for all genes, cohorts and algorithms.
    norm = Normalize(vmin=0, vmax=5, clip=False)
    cmap = LinearSegmentedColormap.from_list('Ocean_gold', list(zip(
        [0,.14,.40,.66,.84,1], ['#E9EEF3','#C2D6E1','#215B7B','#169E9A','#A7C77D','#F3D46A'])), N=256)
    colour_limits = []
    for ds,d in dat.items():
        for gene in genes:
            values = d[gene+'_logCP10k'].to_numpy()
            colour_limits.append({'dataset':ds, 'gene':gene, 'vmin':0, 'vmax':5,
                'observed_min':float(values.min()), 'observed_max':float(values.max()),
                'clipped_cells':0, 'units':'ln(1 + 10000 * gene UMI / all-source-gene total UMI)',
                'mapping':'One shared Normalize(0,5), same Ocean_gold colormap, no cohort/gene rescaling'})
    pd.DataFrame(colour_limits).to_csv(metadata/'SHARED_COLOUR_SCALE.csv', index=False, mode='x')
    cohort = {'GSE289422':'CAD cohort', 'GSE255566':'T2D cohort'}
    cohortcolor = {'GSE289422':'#155B8A', 'GSE255566':'#C55432'}
    ledger, alignment, files = [], [], []
    for mode,cols in [('UMAP',['UMAP1','UMAP2']), ('tSNE',['TSNE1','TSNE2'])]:
        centers = {ds:((d[cols].max()+d[cols].min())/2).to_numpy() for ds,d in dat.items()}
        span = max(float(np.ptp(d[col])) for d in dat.values() for col in cols)*1.12
        # Equal coordinate units and shared bounds retain every cell.
        xspan = max(float(np.ptp(d[cols[0]])) for d in dat.values())*1.045
        yspan = max(float(np.ptp(d[cols[1]])) for d in dat.values())*1.045
        width_mm, left_mm, right_mm, gap_mm = 180., 11., 2., 3.
        panel_width_mm = (width_mm-left_mm-right_mm-3*gap_mm)/4
        panel_height_mm = panel_width_mm*yspan/xspan
        legend_mm, rowgap_mm, header_mm = 12., 2., 7.
        height_mm = legend_mm+2*panel_height_mm+rowgap_mm+header_mm
        fw,fh = width_mm/25.4,height_mm/25.4
        fig = plt.figure(figsize=(fw,fh))
        left,pw,gap = left_mm/width_mm,panel_width_mm/width_mm,gap_mm/width_mm
        ph = panel_height_mm/height_mm
        rowbottom = [(legend_mm+panel_height_mm+rowgap_mm)/height_mm,legend_mm/height_mm]
        axes = []
        for row,ds in enumerate(dslist):
            d = dat[ds]
            xy = d[cols].to_numpy()-centers[ds]
            fig.text(.030, rowbottom[row]+ph/2, cohort[ds]+'\n'+ds,
                fontsize=7.5,fontweight='bold',color=cohortcolor[ds],rotation=90,
                ha='center',va='center')
            for col,gene in enumerate(genes):
                ax = fig.add_axes([left+col*(pw+gap),rowbottom[row],pw,ph])
                axes.append((row,col,ax))
                ax.set_aspect('equal',adjustable='box')
                ax.set_xlim(-xspan/2,xspan/2); ax.set_ylim(-yspan/2,yspan/2)
                ax.set_xticks([]); ax.set_yticks([])
                for spine in ax.spines.values():
                    spine.set_visible(False)
                if row == 0:
                    ax.set_title(gene,fontsize=8.5,fontstyle='italic',fontweight='bold',pad=4)
                values = d[gene+'_logCP10k'].to_numpy()
                ix = np.argsort(values,kind='stable')
                ax.scatter(xy[ix,0],xy[ix,1],c=values[ix],s=.68,alpha=1,linewidths=0,
                    rasterized=True,cmap=cmap,norm=norm,zorder=1)
                # Alpha=1 makes the numeric shared colorbar exactly match each dot.
                ledger.append({'algorithm':mode,'dataset':ds,'gene':gene,
                    'display_cells':len(d),'dx':float(-centers[ds][0]),'dy':float(-centers[ds][1]),
                    'scale_x':1,'scale_y':1,'rotation_degrees':0,
                    'xlim':[-xspan/2,xspan/2], 'ylim':[-yspan/2,yspan/2],
                    'native_coordinates_translated_only':True})

        # A single orientation key replaces all repeated panel axis labels.
        key = fig.add_axes([11./width_mm,3.2/height_mm,12./width_mm,6.8/height_mm]); key.set_axis_off()
        arrow = dict(arrowstyle='-|>',color='#444444',lw=.65,mutation_scale=6)
        key.annotate('',xy=(.98,.18),xytext=(.16,.18),arrowprops=arrow)
        key.annotate('',xy=(.16,.96),xytext=(.16,.18),arrowprops=arrow)
        key.text(.57,-.015,('UMAP' if mode=='UMAP' else 't-SNE')+' 1',
                 fontsize=6.7,ha='center',va='top')
        key.text(.015,.60,('UMAP' if mode=='UMAP' else 't-SNE')+' 2',
                 fontsize=6.7,ha='right',va='center',rotation=90)
        # One colorbar for all eight panels; tick labels correspond to data values.
        cax = fig.add_axes([.405,7./height_mm,.54,1.4/height_mm])
        cb = fig.colorbar(ScalarMappable(norm=norm,cmap=cmap),cax=cax,
                         orientation='horizontal',ticks=[0,1,2,3,4,5])
        cb.outline.set_linewidth(.4)
        cb.ax.tick_params(labelsize=7,length=2,width=.4,pad=1.5)
        fig.text(.675,.6/height_mm,'Expression [ln(1 + CP10k)]',ha='center',va='bottom',fontsize=7.3)

        fig.canvas.draw()
        for row in [0,1]:
            bounds=[a.get_position().bounds for r,c,a in axes if r==row]
            assert np.ptp([b[1] for b in bounds]) < 1e-10
            assert np.ptp([b[2] for b in bounds]) < 1e-10
            assert np.ptp([b[3] for b in bounds]) < 1e-10
        for col in range(4):
            bounds=[a.get_position().bounds for r,c,a in axes if c==col]
            assert np.ptp([b[0] for b in bounds]) < 1e-10
        panel_bounds=[list(map(float,a.get_position().bounds)) for r,c,a in axes]
        alignment.append({'algorithm':mode,'canvas_mm':[width_mm,height_mm],
            'eight_equal_rectangular_panels':True,'row_column_edges_aligned':True,
            'panel_bounds':panel_bounds,'common_coordinate_xlim':[-xspan/2,xspan/2], 'common_coordinate_ylim':[-yspan/2,yspan/2],
            'coordinate_scale_x_y_equal':True,'shared_colour_scale':[0,5],
            'colourbar_count':1,'axis_key_count':1,'panel_letter_count':0})
        stem='Figure5_features_'+mode
        for ext,dpi in [('.pdf',600),('.svg',600),('_300dpi.png',300),('_600dpi.png',600)]:
            path=figures/(stem+ext)
            fig.savefig(path,dpi=dpi)
            files.append({'path':str(path),'sha256':sha(path),'bytes':path.stat().st_size})
        plt.close(fig)
        print('Saved',stem,flush=True)

    pd.DataFrame(ledger).to_csv(metadata/'COORDINATE_LEDGER.csv',index=False,mode='x')
    for path in sources:
        assert sha(path)==hashes[str(path)],f'Input changed: {path}'
    write_json(metadata/'RENDER_MANIFEST.json',{'status':'COMPLETE',
        'input_sha256':hashes,'script':str(Path(__file__).resolve()),'script_sha256':sha(Path(__file__).resolve()),
        'expression_palette':'Ocean_gold',
        'expression_values_unchanged':True,'cell_selection_unchanged':True,
        'source_coordinates_unchanged':True,'no_model_refits':True,'no_new_tests':True,
        'new_normalization':False,'cross_dataset_disease_comparison_estimated':False,
        'display_cell_counts':{ds:len(d) for ds,d in dat.items()},
        'maximum_expression':max(r['observed_max'] for r in colour_limits),
        'clipped_cells_all_panels':0,'alignment':alignment,'outputs':files})


if __name__=='__main__':
    main()
