#!/usr/bin/env python3
"""Build Fig. 7 with landmark schematic A, pooled summary B, and cases C-F.

The original four examples are copied as SVG objects and arranged in one row
without scaling. Only their panel letters change from ver10 A-D to C-F. Counts are
independently checked against every eligible row of S2 Data.
"""
from pathlib import Path
from io import BytesIO
import copy
import csv
import hashlib
import json
import math
import os
import subprocess
import unicodedata

from lxml import etree as ET
import fitz
from PIL import Image

ROOT = Path(__file__).resolve().parent
SOURCE = ROOT / 'source'
SVG_NS = 'http://www.w3.org/2000/svg'
INK_NS = 'http://www.inkscape.org/namespaces/inkscape'
NS = {'s': SVG_NS}
CATS = ['preserved_all3', 'reversed_all3', 'spatial_consistent_temporal_unresolved', 'mixed', 'other_unresolved']
COLORS = ['#326F99', '#AD483B', '#6D6299', '#E4BE62', '#DDE1E5']
TEXT_COLORS = ['#FFFFFF', '#FFFFFF', '#FFFFFF', '#202B34', '#202B34']
GROUPS = ['pair-rule', 'gap-gap', 'gap-pair-rule']
LABELS = ['Pair-rule pairs', 'Gap-gene pairs', 'Gap / pair-rule pairs']
EXPECTED = {'pair-rule': [32,0,0,1,0], 'gap-gap': [28,3,5,4,1],
            'gap-pair-rule': [117,13,9,15,5]}
WIDTH = 1080
TOP_HEIGHT = 196
HEIGHT = 640
FONT = 'Arial, Liberation Sans, Nimbus Sans, sans-serif'


def atomic_write(path, payload):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_name(path.name + '.writing')
    with tmp.open('wb') as stream:
        stream.write(payload)
        stream.flush()
        os.fsync(stream.fileno())
    os.replace(tmp, path)
    assert path.stat().st_size == len(payload)


def add(parent, name, attrs=None, text=None):
    node = ET.SubElement(parent, '{'+SVG_NS+'}'+name,
                         {k: str(v) for k,v in (attrs or {}).items()})
    if text is not None:
        node.text = text
    return node


def text(parent, x, y, value, size=8, weight='normal', anchor='start', color='#222222'):
    return add(parent, 'text', {'x':x, 'y':y, 'font-family':FONT,
        'font-size':size, 'font-weight':weight, 'font-style':'normal',
        'text-anchor':anchor, 'fill':color}, value)


def schematic(layer):
    panel=add(layer,'g',{'id':'landmark_schematic'})
    text(panel,9,16,'A',12,'bold')
    text(panel,26,16,'Expression landmarks',10,'bold')
    x0,x1,top,bottom=52,251,42,126
    center=151.5
    sigma=26
    half_y=(top+bottom)/2
    half=math.sqrt(2*math.log(2))*sigma
    # An illustrative Gaussian pulse is used only to define the three landmarks.
    points=[]
    for i in range(251):
        x=x0+(x1-x0)*i/250
        value=math.exp(-0.5*((x-center)/sigma)**2)
        points.append((x,bottom-(bottom-top)*value))
    for y,label in [(bottom,'0'),((top+bottom)/2,'0.5'),(top,'1')]:
        add(panel,'line',{'x1':x0-3,'x2':x0,'y1':y,'y2':y,'stroke':'#555555','stroke-width':.65})
        text(panel,x0-5,y+2.6,label,7,anchor='end',color='#555555')
    add(panel,'line',{'x1':x0,'x2':x0,'y1':top-2,'y2':bottom,'stroke':'#555555','stroke-width':.7})
    add(panel,'line',{'x1':x0,'x2':x1,'y1':bottom,'y2':bottom,'stroke':'#555555','stroke-width':.7})
    add(panel,'line',{'x1':x0,'x2':x1,'y1':half_y,'y2':half_y,'stroke':'#A5ADB4','stroke-width':.7,'stroke-dasharray':'2 2'})
    path='M '+' L '.join(f'{x:.5f},{y:.5f}' for x,y in points)
    add(panel,'path',{'d':path,'fill':'none','stroke':'#505962','stroke-width':1.6})
    for x,y in [(center,top),(center-half,half_y),(center+half,half_y)]:
        add(panel,'circle',{'cx':x,'cy':y,'r':2.4,'fill':'#222222','stroke':'#FFFFFF','stroke-width':.6})
    text(panel,center,33,'Peak',8,'bold','middle')
    text(panel,91,62,'Rising',8,anchor='middle')
    text(panel,91,72,'half-height',8,anchor='middle')
    text(panel,212,62,'Falling',8,anchor='middle')
    text(panel,212,72,'half-height',8,anchor='middle')
    add(panel,'line',{'x1':108,'x2':center-half-2,'y1':76,'y2':half_y-2,'stroke':'#777777','stroke-width':.6})
    add(panel,'line',{'x1':195,'x2':center+half+2,'y1':76,'y2':half_y-2,'stroke':'#777777','stroke-width':.6})
    lab=text(panel,17,84,'Normalized expression',8,anchor='middle',color='#555555')
    lab.set('transform','rotate(-90 17 84)')
    text(panel,center,142,'Time or position',8,anchor='middle')
    text(panel,center,155,'Illustrative profile',7.5,anchor='middle',color='#666666')
    return panel


def main():
    baseline_path = SOURCE/'Fig7_ver10.svg'
    baseline = ET.parse(str(baseline_path)).getroot()
    old_layer = baseline.find('s:g', NS)
    assert old_layer is not None
    with (SOURCE/'summary_counts_and_percentages.csv').open(newline='') as stream:
        rows=list(csv.DictReader(stream))
    with (SOURCE/'S2_Data_Landmark_Comparisons.csv').open(newline='') as stream:
        inventory=list(csv.DictReader(stream))
    eligible=[r for r in inventory if r['eligible'].lower()=='true']
    assert len(inventory)==859 and len(eligible)==233
    verified={}
    for group, expected in EXPECTED.items():
        subset=[r for r in rows if r['panel']=='all_landmarks' and r['group']==group]
        assert len(subset)==5
        counts=[int(next(r['count'] for r in subset if r['category']==cat)) for cat in CATS]
        independently_counted=[sum(r['group']==group and r['summary_category']==cat for r in eligible) for cat in CATS]
        assert counts==expected==independently_counted, (group,counts,expected,independently_counted)
        assert all(int(r['total'])==sum(counts) and r['figure_panel']=='B' for r in subset)
        assert next(r['category_label'] for r in subset if r['category']==CATS[2])=='Unresolved in time, separated in space'
        verified[group]=counts

    root=ET.Element('{'+SVG_NS+'}svg', nsmap={None:SVG_NS,'inkscape':INK_NS},
                    version='1.1', width=f'{WIDTH}pt',height=f'{HEIGHT}pt',
                    viewBox=f'0 0 {WIDTH} {HEIGHT}')
    add(root,'title',text='Figure 7: temporal-to-spatial landmark ordering')
    add(root,'desc',text='Version 12. Panel A illustrates the peak and rising/falling half-height landmarks. Panel B pools all 233 eligible landmark comparisons. Panels C-F retain the four version-10 examples, arranged horizontally without scaling. Live text and vector paths in one editable layer.')
    layer=add(root,'g',{'id':'Figure7', '{'+INK_NS+'}groupmode':'layer', '{'+INK_NS+'}label':'Figure 7'})
    add(layer,'rect',{'x':0,'y':0,'width':WIDTH,'height':HEIGHT,'fill':'#FFFFFF'})
    schematic(layer)
    summary=add(layer,'g',{'id':'landmark_summary'})
    text(summary,279,16,'B',12,'bold')
    text(summary,296,16,'Landmark ordering',10,'bold')
    text(summary,296,30,'Peaks and both half-height borders (n = 233)',8)
    x0=400; w=665; centers=[52,80,108]; bar_h=18
    for p in [0,25,50,75,100]:
        x=x0+w*p/100
        add(summary,'line',{'x1':x,'x2':x,'y1':42,'y2':120,'stroke':'#E7EAED','stroke-width':.5})
        text(summary,x,132,f'{p}%',7,anchor='middle',color='#555555')
    for group,label,y in zip(GROUPS,LABELS,centers):
        counts=verified[group]; total=sum(counts)
        text(summary,279,y-1,label,8)
        text(summary,279,y+8,'n = '+str(total),7.5,color='#555555')
        left=x0
        for count,color,txt_color in zip(counts,COLORS,TEXT_COLORS):
            if count:
                seg_w=w*count/total
                add(summary,'rect',{'x':left,'y':y-bar_h/2,'width':seg_w,'height':bar_h,
                                   'fill':color,'stroke':'#FFFFFF','stroke-width':.35})
                text(summary,left+seg_w/2,y+2.65,str(count),8.5,'bold','middle',txt_color)
                left += seg_w
    text(summary,732.5,148,'Bar widths show percentages; numbers show comparison counts.',8,anchor='middle',color='#555555')
    legend=[
        (9,170,0,['Preserved order in all 3 embryos']),
        (215,170,1,['Reversed order in all 3 embryos']),
        (421,170,2,['Unresolved in time,', 'separated in space']),
        (868,170,4,['Other unresolved comparisons']),
        (627,170,3,['Spatial order differs among embryos;', 'temporal order resolved']),
    ]
    for x,y,i,lines in legend:
        add(summary,'rect',{'x':x,'y':y-6,'width':7,'height':7,'fill':COLORS[i]})
        for j,line in enumerate(lines):
            text(summary,x+11,y+j*9,line,7.8)
    text(summary,540,190,'Each count is one landmark comparison at one stage, assessed across 3 embryos.',7.5,anchor='middle',color='#555555')

    preserved=add(layer,'g',{'id':'original_examples'})
    letter_map={'A':'C','B':'D','C':'E','D':'F'}
    changed=[]
    copied=[]
    # Case ranges are fixed by the source's panel headers. Each complete panel is
    # translated as one block; no original element coordinates or styles change.
    layout=[('C',0,130,0,196),('D',130,256,0,196),
            ('E',256,358,540,-219),('F',358,466,540,-219),
            ('case_legend',466,470,270,-216)]
    source_objects=list(old_layer)
    for label,start,end,tx,ty in layout:
        block=add(preserved,'g',{'id':'case_'+label,'transform':f'translate({tx} {ty})'})
        if label!='case_legend':
            assert source_objects[start].text==dict(C='A',D='B',E='C',F='D')[label]
            assert source_objects[start].get('font-size')=='12'
        for child in source_objects[start:end]:
            node=copy.deepcopy(child)
            if ET.QName(node).localname=='text' and node.text in letter_map and node.get('font-size')=='12':
                changed.append((node.text,letter_map[node.text]))
                node.text=letter_map[node.text]
            block.append(node)
            copied.append(node)
    assert changed==[('A','C'),('B','D'),('C','E'),('D','F')],changed
    original_paths=old_layer.findall('s:path',NS)
    assert [dict(x.attrib) for x in original_paths]==[dict(x.attrib) for x in copied if ET.QName(x).localname=='path']
    a=list(old_layer);b=copied
    assert len(a)==len(b)==470
    for original,new in zip(a,b):
        assert dict(original.attrib)==dict(new.attrib)
        assert original.tag==new.tag
        if original.text!=new.text:
            assert (original.text,new.text) in changed and original.get('font-size')=='12'
    svg_path=ROOT/'Fig7_ver12.svg'
    atomic_write(svg_path,ET.tostring(root,xml_declaration=True,encoding='UTF-8',pretty_print=True))
    temp_pdf=ROOT/'qa'/'Fig7_ver12_export.pdf'
    temp_pdf.parent.mkdir(exist_ok=True)
    subprocess.run(['inkscape',str(svg_path),'--export-type=pdf',f'--export-filename={temp_pdf}'],capture_output=True,text=True,check=True)
    payload=temp_pdf.read_bytes();assert payload.startswith(b'%PDF') and len(payload)>20000
    atomic_write(ROOT/'Fig7_ver12.pdf',payload)
    doc=fitz.open(stream=payload,filetype='pdf');assert len(doc)==1
    page=doc[0];assert abs(page.rect.width-WIDTH)<.01 and abs(page.rect.height-HEIGHT)<.01
    assert len(page.get_images())==0,'PDF must remain vector'
    page_text=unicodedata.normalize('NFKC',page.get_text())
    for required in ['Expression landmarks','Illustrative profile','Unresolved in time,','separated in space','Peak-time differences','Embryo 3']:
        assert required in page_text,required
    pix=page.get_pixmap(matrix=fitz.Matrix(2,2),alpha=False)
    atomic_write(ROOT/'Fig7_ver12.png',pix.tobytes('png'))
    high=page.get_pixmap(matrix=fitz.Matrix(600/72,600/72),alpha=False)
    im=Image.frombytes('RGB',(high.width,high.height),high.samples)
    buf=BytesIO();im.save(buf,format='TIFF',compression='tiff_lzw',dpi=(600,600))
    atomic_write(ROOT/'Fig7_ver12.tiff',buf.getvalue())
    with Image.open(ROOT/'Fig7_ver12.tiff') as saved:
        saved.load();assert saved.mode=='RGB';assert saved.info['dpi']==(600.,600.)
    stats={'dimensions_pt':[WIDTH,HEIGHT],'dimensions_in':[WIDTH/72,HEIGHT/72],
        'baseline_dimensions_pt':[540,856],'case_layout':'four columns in one row, C-D-E-F','case_block_translations':{label:[tx,ty] for label,start,end,tx,ty in layout},
        'original_scaling':1.0,'original_svg_elements':len(a),
        'original_paths_unchanged':len(original_paths),'all_original_object_attributes_unchanged':True,
        'only_original_text_changes_from_ver10':changed,'case_panel_letters_unchanged_from_ver11':True,'case_blocks_preserved_without_scaling':True,
        'single_layer':len(root.findall('.//s:g[@inkscape:groupmode="layer"]',{'s':SVG_NS,'inkscape':INK_NS}))==1,
        'pdf_embedded_images':len(page.get_images()),'pdf_drawing_count':len(page.get_drawings()),
        'tiff_pixels':[high.width,high.height],'tiff_dpi':600,
        'counts':verified,'pooled_comparisons':len(eligible),'candidate_rows':len(inventory),
        'schematic':'Illustrative Gaussian pulse with peak and exactly 50% rising/falling points; no biological data.',
        'purple_category_label':'Unresolved in time, separated in space',
        'source_sha256':{p.name:hashlib.sha256(p.read_bytes()).hexdigest() for p in SOURCE.iterdir() if p.is_file()},
        'output_bytes':{p.name:p.stat().st_size for p in ROOT.glob('Fig7_ver12.*')}}
    atomic_write(ROOT/'qa'/'validation.json',json.dumps(stats,indent=2).encode())
    print(json.dumps(stats,indent=2))


if __name__=='__main__':
    main()
