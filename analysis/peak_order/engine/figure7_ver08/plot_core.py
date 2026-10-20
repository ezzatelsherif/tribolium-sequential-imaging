"""Shared peak-position plotting data and editable vector primitives for version 08."""
from pathlib import Path
import html, json, math, subprocess, re, os
from functools import lru_cache
import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parent
OLD = ROOT.parent / 'figure7_ver04_analysis'
NEW = ROOT.parent / 'figure7_ver05_analysis'
COL = {"kr":"#ED2124","gt":"#C89C0A","mlpt":"#0D8140","svb":"#A62A2A",
       "eve":"#7D277D","run":"#F06CA8","odd":"#6FCDDD"}
LABELS = [f'{c}.{s}' for c in range(1,9) for s in range(1,4)]
STAGES = ['Stage3_3','Stage4_1','Stage4_3','Stage5_2']

def normalize(y):
    y = np.asarray(y, float)
    lo, hi = float(np.nanmin(y)), float(np.nanmax(y))
    return (y-lo)/(hi-lo) if hi > lo else np.zeros_like(y)

@lru_cache(None)
def spatial_curve(stage, embryo, gene):
    if gene in ['eve','run','odd']:
        src = OLD/'pairrule/spatial_profiles.csv'; value = 'normalized_intensity'
    elif gene in ['gt','mlpt']:
        src = OLD/'gt_mlpt/spatial_profiles.csv'; value = 'normalized_intensity'
    else:
        src = OLD/f'{gene}_eve/spatial_profiles.csv'; value = 'smoothed_intensity'
    d = pd.read_csv(src)
    d = d[(d.stage==stage)&(d.embryo==embryo)&(d.gene.str.lower()==gene)].sort_values('x_ap')
    return d.x_ap.to_numpy(), normalize(d[value].to_numpy()), str(src)

@lru_cache(None)
def temporal_curve(cohort, gene):
    if cohort == 'pairrule':
        src = OLD/'pairrule/temporal_profiles.csv'; value = 'normalized_mean'
    else:
        src = NEW/'gap_peaks/temporal_profiles.csv'; value = 'normalized_intensity'
    d = pd.read_csv(src)
    d = d[d.gene.str.lower()==gene]
    if 'cohort' in d and (d.cohort.str.lower()==cohort).any():
        d=d[d.cohort.str.lower()==cohort]
    d=d.sort_values('substage_index')
    return d.substage_index.to_numpy(), normalize(d[value].to_numpy()), str(src)

def crop(x,y,lo,hi):
    """Crop without masks, including true interpolated crop endpoints."""
    keep=(x>lo)&(x<hi)
    xx=np.r_[max(lo,x[0]),x[keep],min(hi,x[-1])]
    xx=np.unique(xx)
    return xx,np.interp(xx,x,y)

class VectorFigure:
    def __init__(self,w,h):
        self.doc={'width':w,'height':h,'layer_name':'Figure 7','objects':[]}
        self.o=self.doc['objects']
    def path(self,pts,stroke='#222222',width=.65,fill=None,closed=False):
        self.o.append(dict(type='path',points=[[float(x),float(y)] for x,y in pts],stroke=stroke,width=width,fill=fill,closed=closed,dash=[]))
    def line(self,x,y,x2,y2,c='#222222',w=.6):self.path([(x,y),(x2,y2)],c,w)
    def text(self,x,y,s,size=8,color='#222222',anchor='start',bold=False,italic=False,rotation=0):
        """Colour each gene name individually, including names inside comparisons."""
        s=str(s)
        pattern=re.compile(r'\b(mlpt|eve|run|odd|svb|gt|kr)([0-9]*)\b',re.I)
        matches=list(pattern.finditer(s))
        if not matches:
            self.o.append(dict(type='text',x=float(x),y=float(y),text=s,font='Arial',size=size,fill=color,anchor=anchor,bold=bold,italic=italic,rotation=rotation))
            return
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        suffix='BoldItalic' if bold and italic else 'Bold' if bold else 'Italic' if italic else 'Regular'
        font='FigureSans'+suffix
        if font not in pdfmetrics.getRegisteredFontNames():
            fontdir=ROOT.parent/'fonts'
            pdfmetrics.registerFont(TTFont(font,str(fontdir/f'LiberationSans-{suffix}.ttf')))
        width=pdfmetrics.stringWidth(s,font,size)
        offset=-width/2 if anchor=='middle' else -width if anchor=='end' else 0
        theta=math.radians(rotation)
        segments=[];end=0
        for m in matches:
            if m.start()>end:segments.append((s[end:m.start()],color))
            segments.append((m.group(0),COL[m.group(1).lower()]))
            end=m.end()
        if end<len(s):segments.append((s[end:],color))
        for text,fill in segments:
            # Leading/trailing XML whitespace is collapsed by some SVG readers.
            # Encode its width in position, so all readers keep word spacing.
            leading=text[:len(text)-len(text.lstrip(' '))]
            drawn=text.strip(' ')
            left=offset+pdfmetrics.stringWidth(leading,font,size)
            if drawn:
                self.o.append(dict(type='text',x=float(x+left*math.cos(theta)),y=float(y+left*math.sin(theta)),
                                   text=drawn,font='Arial',size=size,fill=fill,anchor='start',bold=bold,italic=italic,rotation=rotation))
            offset+=pdfmetrics.stringWidth(text,font,size)
    def marker(self,x,y,color,kind='circle',r=2.35,filled=True):
        if kind=='diamond': pts=[(x,y-r*1.2),(x+r,y),(x,y+r*1.2),(x-r,y)]
        else: pts=[(x+r*math.cos(t),y+r*math.sin(t)) for t in np.linspace(0,2*math.pi,24,endpoint=False)]
        self.path(pts,color,.7,color if filled else '#ffffff',True)
    def axes(self,x,y,w,h,xlim,ticks,xlabel='',embryo=None,ylabel=False):
        fx=lambda v:x+(v-xlim[0])/(xlim[1]-xlim[0])*w
        fy=lambda v:y+h-v/1.12*h
        self.line(x,y+h,x+w,y+h);self.line(x,y,x,y+h)
        for v,l in ticks:
            xx=fx(v);self.line(xx,y+h,xx,y+h+2.4)
            self.text(xx,y+h+12,l,8,anchor='middle')
        for v in [0,1]:
            yy=fy(v);self.line(x-2.4,yy,x,yy)
            self.text(x-5,yy+2.5,str(v),8,anchor='end',color='#666666')
        if xlabel:self.text(x+w/2,y+h+26,xlabel,8,anchor='middle')
        if embryo:self.text(x+w-3,y-4,f'Embryo {embryo}',8,color='#666666',anchor='end')
        if ylabel:self.text(x-24,y+h/2,'Relative expression',8,anchor='middle',rotation=-90)
        return fx,fy
    def curve(self,fx,fy,x,y,color):self.path([(fx(a),fy(b)) for a,b in zip(x,y)],color,1.05)
    def save(self,stem):
        """Write standard SVG with live text and a direct vector PDF.

        PDF creation does not use the earlier legacy Illustrator exporter.
        SVG holds all paths and text in one layer, without nested groups or masks.
        """
        from reportlab.pdfgen import canvas
        from reportlab.pdfbase import pdfmetrics
        from reportlab.pdfbase.ttfonts import TTFont
        from reportlab.lib.colors import HexColor
        import os
        name=ROOT/stem
        name.with_suffix('.json').write_text(json.dumps(self.doc,indent=2))
        w,h=self.doc['width'],self.doc['height']
        svg=['<?xml version="1.0" encoding="UTF-8"?>',
             f'<svg xmlns="http://www.w3.org/2000/svg" xmlns:inkscape="http://www.inkscape.org/namespaces/inkscape" version="1.1" width="{w}pt" height="{h}pt" viewBox="0 0 {w} {h}">',
             '<title>Figure 7: temporal and spatial expression peak order</title>',
             '<desc>Version 08. Gene colours follow submitted Figure 6 with darker gt yellow; every gene name follows its own colour. Editable vector paths and live text, one layer. Analysis and layout retained from version 06.</desc>',
             '<g id="Figure7" inkscape:groupmode="layer" inkscape:label="Figure 7">']
        for ob in self.o:
            if ob['type']=='path':
                d='M '+' L '.join(f'{x:.5f} {y:.5f}' for x,y in ob['points'])+(' Z' if ob['closed'] else '')
                svg.append(f'<path d="{d}" fill="{ob["fill"] or "none"}" stroke="{ob["stroke"] or "none"}" stroke-width="{ob["width"]}" stroke-linejoin="round" stroke-linecap="round"/>')
            else:
                tr=f' transform="rotate({ob["rotation"]} {ob["x"]} {ob["y"]})"' if ob['rotation'] else ''
                svg.append(f'<text x="{ob["x"]}" y="{ob["y"]}" font-family="Arial, Liberation Sans, Nimbus Sans, sans-serif" font-size="{ob["size"]}" font-weight="{"bold" if ob["bold"] else "normal"}" font-style="{"italic" if ob["italic"] else "normal"}" text-anchor="{ob["anchor"]}" fill="{ob["fill"]}"{tr}>{html.escape(ob["text"])}</text>')
        svg+=['</g>','</svg>']
        name.with_suffix('.svg').write_text('\n'.join(svg),encoding='utf-8')
        # Liberation Sans is metrically compatible with Arial and embedded in PDF.
        fontdir=ROOT.parent/'fonts'
        for suffix in ['Regular','Bold','Italic','BoldItalic']:
            pdfmetrics.registerFont(TTFont('FigureSans'+suffix,str(fontdir/f'LiberationSans-{suffix}.ttf')))
        c=canvas.Canvas(str(name.with_suffix('.pdf')),pagesize=(w,h),pageCompression=1,
                        initialFontName='FigureSansRegular',initialFontSize=8)
        c.setTitle('Figure 7 - temporal and spatial expression peak order (version 08)')
        c.setAuthor('Ezzat El-Sherif')
        c.setSubject('Gene names colour-matched individually, gt yellow darkened; curves, peaks and timing intervals retained from version 06')
        for ob in self.o:
            c.saveState()
            if ob['type']=='path':
                pts=ob['points'];p=c.beginPath();p.moveTo(pts[0][0],h-pts[0][1])
                for x,y in pts[1:]:p.lineTo(x,h-y)
                if ob['closed']:p.close()
                if ob['stroke']:c.setStrokeColor(HexColor(ob['stroke']))
                if ob['fill']:c.setFillColor(HexColor(ob['fill']))
                c.setLineWidth(ob['width']);c.setLineCap(1);c.setLineJoin(1)
                c.drawPath(p,stroke=int(bool(ob['stroke'])),fill=int(bool(ob['fill'])))
            else:
                suffix='BoldItalic' if ob['bold'] and ob['italic'] else 'Bold' if ob['bold'] else 'Italic' if ob['italic'] else 'Regular'
                c.setFont('FigureSans'+suffix,ob['size']);c.setFillColor(HexColor(ob['fill']))
                c.translate(ob['x'],h-ob['y']);c.rotate(-ob['rotation'])
                draw=c.drawCentredString if ob['anchor']=='middle' else c.drawRightString if ob['anchor']=='end' else c.drawString
                draw(0,0,ob['text'])
            c.restoreState()
        c.showPage();c.save()
        import fitz
        rendered=fitz.open(name.with_suffix('.pdf'))
        rendered[0].get_pixmap(matrix=fitz.Matrix(150/72,150/72),alpha=False).save(name.with_suffix('.png'))
        stats=dict(paths=sum(o['type']=='path' for o in self.o),text_objects=sum(o['type']=='text' for o in self.o),
                   svg_layers=1,nested_groups=0,clipping_paths=0,masks=0,raster_images=0,
                   width_pt=w,height_pt=h,minimum_text_size=min(o['size'] for o in self.o if o['type']=='text'),
                   pdf_direct_vector_export=True,svg_live_text=True,native_illustrator_reopening_tested=False,
                   peak_position_markers=True,temporal_peak_ci95_bars=True,direct_temporal_difference_ci95=True,
                   halfheight_markers=0,spatial_error_bars=False,gene_palette=COL)
        return stats
