"""Fill the [FIGURE PLACEHOLDER] shapes in the deck with real images.

Usage:
    python scripts/fill_deck_placeholders.py [deck.pptx]

It scans every slide for a text box containing 'FIGURE PLACEHOLDER',
reads the image path from the box, removes the box, and inserts the image
with the same geometry (aspect-ratio preserved, centred).

Expected placeholders:
    plots/Master/S1only/decay_chain.png
    plots/Master/S1only/fitted_mu_trend_master.png
"""
import os
import sys

from PIL import Image
from pptx import Presentation

DECK = sys.argv[1] if len(sys.argv) > 1 else 'Rn222_SR3_S1only_DavidTang_20260918.pptx'

prs = Presentation(DECK)
filled = 0
missing = 0
for si, s in enumerate(prs.slides, 1):
    for sh in list(s.shapes):
        if not (sh.has_text_frame and 'FIGURE PLACEHOLDER' in sh.text_frame.text):
            continue
        path = sh.text_frame.text.strip().split('\n')[-1].strip()
        if not os.path.exists(path):
            print('  slide %d: still missing -> %s' % (si, path))
            missing += 1
            continue
        L, T, W, H = sh.left, sh.top, sh.width, sh.height
        s.shapes._spTree.remove(sh._element)
        iw, ih = Image.open(path).size
        ar = iw / ih
        if (W / H) > ar:
            hh = H
            ww = int(H * ar)
        else:
            ww = W
            hh = int(W / ar)
        s.shapes.add_picture(path, L + (W - ww) // 2, T + (H - hh) // 2, ww, hh)
        print('  slide %d: inserted %s' % (si, path))
        filled += 1

prs.save(DECK)
print('done: filled %d, still missing %d' % (filled, missing))
