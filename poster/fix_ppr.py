"""Keep only the first <a:pPr> in each <a:p> (pptxgenjs writes one per run)."""
import re, sys, zipfile, shutil, os
src = sys.argv[1]
tmp = src + ".tmp"
pat_p = re.compile(r"<a:p>(.*?)</a:p>", re.S)
pat_ppr = re.compile(r"<a:pPr\b[^>]*?/>|<a:pPr\b[^>]*?>.*?</a:pPr>", re.S)
def fix_para(m):
    body = m.group(1)
    found = []
    def keep_first(mm):
        found.append(1)
        return mm.group(0) if len(found) == 1 else ""
    return "<a:p>" + pat_ppr.sub(keep_first, body) + "</a:p>"
with zipfile.ZipFile(src) as zin, zipfile.ZipFile(tmp, "w", zipfile.ZIP_DEFLATED) as zout:
    for item in zin.infolist():
        data = zin.read(item.filename)
        if re.match(r"ppt/slides/slide\d+\.xml$", item.filename):
            data = pat_p.sub(fix_para, data.decode("utf-8")).encode("utf-8")
        zout.writestr(item, data)
shutil.move(tmp, src)
print("fixed", src)
