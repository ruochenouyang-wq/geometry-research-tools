"""Render the reader edition of the research report from its Markdown source."""
from pathlib import Path
import hashlib
import html
import json
import re

from reportlab.lib import colors
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (
    BaseDocTemplate, Frame, PageBreak, PageTemplate, Paragraph, Spacer,
    Table, TableStyle,
)
from pypdf import PdfReader
import reportlab.platypus.paragraph as paragraph_module
import reportlab.lib.textsplit as textsplit_module

# ReportLab's default CJK list primarily covers Japanese punctuation.
# Add Chinese closing punctuation to the in-process wrapping rules.
CHINESE_CLOSING = '，：；！？）》】’”'
paragraph_module.ALL_CANNOT_START += CHINESE_CLOSING
textsplit_module.ALL_CANNOT_START += CHINESE_CLOSING

ROOT = Path(__file__).resolve().parents[2]
SOURCE = Path(__file__).with_name('研究指导报告.md')
DEST = ROOT / 'output/pdf/数学研究工具_长期指导报告.pdf'
WORK = ROOT / 'tmp/pdfs/research_strategy'
WORK.mkdir(parents=True, exist_ok=True)
DEST.parent.mkdir(parents=True, exist_ok=True)
pdfmetrics.registerFont(TTFont('Heiti', '/System/Library/Fonts/STHeiti Light.ttc', subfontIndex=0))
pdfmetrics.registerFont(TTFont('HeitiBold', '/System/Library/Fonts/STHeiti Medium.ttc', subfontIndex=0))
pdfmetrics.registerFontFamily('Heiti', normal='Heiti', bold='HeitiBold', italic='Heiti', boldItalic='HeitiBold')
INK = colors.HexColor('#182C3C')
TEAL = colors.HexColor('#126B69')
MUTED = colors.HexColor('#657682')
PALE = colors.HexColor('#F1F6F5')
PW, PH = A4
WIDTH = PW - 108

STYLES = {
    'body': ParagraphStyle('body', fontName='Heiti', fontSize=10.5, leading=17.1, textColor=INK, spaceAfter=10, wordWrap='CJK', allowWidows=0, allowOrphans=0),
    'h2': ParagraphStyle('h2', fontName='HeitiBold', fontSize=19, leading=27, textColor=TEAL, spaceAfter=17, keepWithNext=True, wordWrap='CJK'),
    'h3': ParagraphStyle('h3', fontName='HeitiBold', fontSize=12, leading=18, textColor=INK, spaceBefore=9, spaceAfter=8, keepWithNext=True, wordWrap='CJK'),
    'bullet': ParagraphStyle('bullet', fontName='Heiti', fontSize=10.4, leading=16.8, textColor=INK, leftIndent=11, firstLineIndent=-10, spaceAfter=8, wordWrap='CJK'),
    'cell': ParagraphStyle('cell', fontName='Heiti', fontSize=9.5, leading=14, textColor=INK, wordWrap='CJK'),
    'head': ParagraphStyle('head', fontName='HeitiBold', fontSize=9.5, leading=14, textColor=colors.white, wordWrap='CJK'),
    'small': ParagraphStyle('small', fontName='Heiti', fontSize=9.6, leading=15.3, textColor=MUTED, spaceAfter=10, wordWrap='CJK'),
    'cover_title': ParagraphStyle('cover_title', fontName='HeitiBold', fontSize=29, leading=42, textColor=INK, spaceAfter=18, wordWrap='CJK'),
    'cover_sub': ParagraphStyle('cover_sub', fontName='HeitiBold', fontSize=16, leading=24, textColor=TEAL, spaceAfter=13, wordWrap='CJK'),
}

def normalize(s):
    s = s.replace('–', '-').replace('—', '-').replace('‑', '-')
    trans = str.maketrans('⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺', '0123456789-+')
    return re.sub(r'10([⁰¹²³⁴⁵⁶⁷⁸⁹⁻⁺]+)', lambda m: '10^(' + m[1].translate(trans) + ')', s)

def inline(s):
    s = normalize(s)
    slots = []
    def link(match):
        label, url = match.groups()
        if url.startswith('https://'):
            val = f'<link href="{html.escape(url, quote=True)}" color="#126B69">{html.escape(label)}</link>'
        else:
            val = html.escape(label) + ' <font color="#657682">(' + html.escape(url) + ')</font>'
        slots.append(val)
        return f'ZZLINK{len(slots)-1}ZZ'
    s = re.sub(r'\[([^\]]+)\]\(([^)]+)\)', link, s)
    s = html.escape(s)
    s = re.sub(r'\*\*([^*]+)\*\*', r'<b>\1</b>', s)
    s = re.sub(r'(?<!\*)\*([^*]+)\*(?!\*)', r'<i>\1</i>', s)
    s = re.sub(r'10\^\((-?\d+)\)', r'10<super>\1</super>', s)
    for i, val in enumerate(slots):
        s = s.replace(f'ZZLINK{i}ZZ', val)
    return s

class Report(BaseDocTemplate):
    def __init__(self):
        super().__init__(str(DEST), pagesize=A4, leftMargin=54, rightMargin=54,
                         topMargin=63, bottomMargin=49,
                         title='以算法放大模型的数学研究能力：面向初学者的修订稿',
                         author='Geometry Research Tools',
                         subject='研究目标、基本概念、工具框架与评估规则')
        frame = Frame(54, 49, WIDTH, PH-112, id='body', leftPadding=0,
                      rightPadding=0, topPadding=0, bottomPadding=0)
        self.addPageTemplates([PageTemplate(id='main', frames=[frame], onPage=self.decor)])
        self.sections = []

    def decor(self, canvas, doc):
        canvas.saveState()
        if doc.page == 1:
            canvas.setFillColor(colors.HexColor('#F7F9F6'))
            canvas.rect(0, 0, PW, PH, fill=1, stroke=0)
            canvas.setFillColor(TEAL)
            canvas.rect(0, 0, 12, PH, fill=1, stroke=0)
            canvas.setFillColor(MUTED)
            canvas.setFont('Heiti', 8.3)
            canvas.drawString(54, PH-40, 'GEOMETRY  /  ALGORITHMS  /  VERIFIED RESEARCH')
        else:
            canvas.setFillColor(MUTED)
            canvas.setFont('Heiti', 8)
            canvas.drawString(54, PH-32, '数学研究工具 · 长期指导报告')
            canvas.drawRightString(PW-54, PH-32, '初学者修订版 1.1')
        canvas.setFillColor(MUTED)
        canvas.setFont('Heiti', 8)
        canvas.drawString(54, 26, '2026.09.18')
        canvas.drawRightString(PW-54, 26, f'{doc.page:02d}')
        canvas.restoreState()

    def afterFlowable(self, flowable):
        if isinstance(flowable, Paragraph) and flowable.style.name == 'h2':
            title = flowable.getPlainText()
            key = f'section{len(self.sections)}'
            self.canv.bookmarkPage(key)
            self.canv.addOutlineEntry(title, key, 0, False)
            self.sections.append({'title': title, 'page': self.page})

lines = SOURCE.read_text().splitlines()
story = []
section = 'cover'
i = 0
while i < len(lines):
    line = lines[i].strip()
    if not line:
        i += 1
        continue
    if line.startswith('# '):
        story.append(Spacer(1, 26))
        story.append(Paragraph('以算法放大模型的<br/>数学研究能力', STYLES['cover_title']))
    elif line.startswith('## '):
        section = line[3:]
        story.extend([PageBreak(), Paragraph(inline(section), STYLES['h2'])])
    elif line.startswith('### '):
        story.append(Paragraph(inline(line[4:]), STYLES['h3']))
    elif line.startswith('|'):
        rows = []
        while i < len(lines) and lines[i].strip().startswith('|'):
            cells = [x.strip() for x in lines[i].strip().strip('|').split('|')]
            if not all(re.fullmatch(r':?-+:?', x) for x in cells):
                rows.append(cells)
            i += 1
        n = len(rows[0])
        weights = {2: [.34, .66], 3: [.22, .30, .48], 4: [.27, .23, .23, .27]}.get(n, [1/n]*n)
        if section.startswith('导读'):
            weights = [.25, .20, .55]
        if section.startswith('15｜'):
            weights = [.20, .38, .42]
        if section.startswith('4｜'):
            weights = [.46, .54]
        data = [[Paragraph(inline(cell), STYLES['head'] if row == 0 else STYLES['cell'])
                 for cell in cells] for row, cells in enumerate(rows)]
        table = Table(data, colWidths=[WIDTH*w for w in weights], repeatRows=1, hAlign='LEFT')
        table.setStyle(TableStyle([
            ('BACKGROUND', (0,0), (-1,0), TEAL), ('VALIGN', (0,0), (-1,-1), 'TOP'),
            ('LEFTPADDING', (0,0), (-1,-1), 7), ('RIGHTPADDING', (0,0), (-1,-1), 7),
            ('TOPPADDING', (0,0), (-1,-1), 7), ('BOTTOMPADDING', (0,0), (-1,-1), 7),
            ('ROWBACKGROUNDS', (0,1), (-1,-1), [PALE, colors.white]),
            ('LINEBELOW', (0,-1), (-1,-1), .4, colors.HexColor('#D8E4E2')),
        ]))
        story.extend([table, Spacer(1, 12)])
        continue
    elif line.startswith('- '):
        story.append(Paragraph('• ' + inline(line[2:]), STYLES['bullet']))
    elif re.match(r'^\d+\. ', line):
        story.append(Paragraph(inline(line), STYLES['bullet']))
    else:
        if section == 'cover' and line == '长期研究指导报告':
            style = STYLES['cover_sub']
        elif section == 'cover' and (line.startswith('2026 年') or line == '面向初学者的修订稿'):
            style = STYLES['small']
        elif section == 'cover' and line in ['**我们希望实现什么**', '**阅读这份报告后，应能回答三个问题**']:
            style = STYLES['h3']
        elif section.startswith('附录'):
            style = STYLES['small']
        else:
            style = STYLES['body']
        story.append(Paragraph(inline(line), style))
    i += 1

doc = Report()
doc.build(story)
reader = PdfReader(DEST)
visible = normalize(re.sub(r'\]\([^)]*\)', ']', SOURCE.read_text()))
missing = sorted({c for c in visible if ord(c)>32 and ord(c) not in pdfmetrics.getFont('Heiti').face.charToGlyph})
quality = {
    'pdf_pages': len(reader.pages), 'sections': doc.sections,
    'source_sha256': hashlib.sha256(SOURCE.read_bytes()).hexdigest(),
    'pdf_sha256': hashlib.sha256(DEST.read_bytes()).hexdigest(),
    'font_missing_characters': missing,
    'external_link_annotations': sum(1 for p in reader.pages for a in p.get('/Annots', []) if a.get_object().get('/A', {}).get('/URI')),
}
(WORK/'BUILD.json').write_text(json.dumps(quality, ensure_ascii=False, indent=2)+'\n')
print(json.dumps(quality, ensure_ascii=False, indent=2))
