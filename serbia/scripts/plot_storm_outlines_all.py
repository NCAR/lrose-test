#!/usr/bin/env python3
"""Plot TITAN storm outlines for 2017-08-12 16:03:32 using Pillow."""
import math
from pathlib import Path
from PIL import Image, ImageDraw, ImageFont

TARGET = (2017, 8, 12, 16, 3, 32)
FILES = ("tracks_anna", "tracks_anna_2", "tracks_mike")

def matching_records(path):
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"): continue
        f = line.split()
        if tuple(map(int, f[3:9])) == TARGET:
            yield {"source": path.name, "simple": int(f[2]),
                   "xc": float(f[13]), "yc": float(f[14]),
                   "rays": list(map(float, f[-72:])), "raw": line}

def main():
    base = Path(__file__).resolve().parent
    records, seen, duplicates = [], set(), 0
    for name in FILES:
        for record in matching_records(base / name):
            if record["raw"] in seen: duplicates += 1
            else: seen.add(record["raw"]); records.append(record)
    for r in records:
        r["points"] = [(r["xc"] + radius * math.sin(math.radians(i * 5)),
                        r["yc"] + radius * math.cos(math.radians(i * 5)))
                       for i, radius in enumerate(r["rays"])]
    xs = [x for r in records for x, _ in r["points"]]
    ys = [y for r in records for _, y in r["points"]]
    xmin, xmax, ymin, ymax = min(xs)-18, max(xs)+18, min(ys)-18, max(ys)+18
    width = height = 1200
    # Keep the plotting box square (y/x aspect ratio = 1) while retaining
    # identical pixel-per-kilometer scaling in both directions.
    left, right, top, bottom = 100, 60, 67, 93
    pw, ph = width-left-right, height-top-bottom
    scale = min(pw/(xmax-xmin), ph/(ymax-ymin))
    xe, ye = (pw-(xmax-xmin)*scale)/2, (ph-(ymax-ymin)*scale)/2
    def pixel(p):
        return (left+xe+(p[0]-xmin)*scale, top+ye+(ymax-p[1])*scale)
    image = Image.new("RGB", (width, height), "white")
    draw = ImageDraw.Draw(image, "RGBA")
    font, small = ImageFont.load_default(size=18), ImageFont.load_default(size=14)
    tick = 50
    for x in range(math.ceil(xmin/tick)*tick, int(xmax)+1, tick):
        px, _ = pixel((x, ymin)); draw.line((px,top,px,height-bottom), fill=(220,220,220,255), width=1)
        draw.text((px-10,height-bottom+9), str(x), fill=(55,55,55,255), font=small)
    for y in range(math.ceil(ymin/tick)*tick, int(ymax)+1, tick):
        _, py = pixel((xmin, y)); draw.line((left,py,width-right,py), fill=(220,220,220,255), width=1)
        draw.text((left-47,py-7), str(y), fill=(55,55,55,255), font=small)
    styles = {"anna": ((23,105,170,40),(23,105,170,255)),
              "mike": ((209,73,91,40),(209,73,91,255))}
    for r in records:
        group = "mike" if r["source"] == "tracks_mike" else "anna"
        fill, outline = styles[group]; points = [pixel(p) for p in r["points"]]
        draw.polygon(points, fill=fill); draw.line(points+[points[0]], fill=outline, width=2, joint="curve")
        cx, cy = pixel((r["xc"],r["yc"])); draw.ellipse((cx-2,cy-2,cx+2,cy+2), fill=outline)
        draw.text((cx, cy), str(r["simple"]), fill=outline, font=small,
                  anchor="mm", stroke_width=2, stroke_fill=(255,255,255,230))
    draw.rectangle((left,top,width-right,height-bottom), outline=(55,55,55,255), width=2)
    draw.text((width/2,23), "TITAN storm outlines - 2017-08-12 16:03:32", fill=(20,20,20,255), font=font, anchor="ma")
    draw.text((width/2,height-47), "TITAN X (km)", fill=(20,20,20,255), font=font, anchor="ma")
    draw.text((48,height/2), "Y (km)", fill=(20,20,20,255), font=font, anchor="mm")
    lx, ly = width-315, 83
    draw.rounded_rectangle((lx,ly,width-57,ly+77), radius=8, fill=(255,255,255,225), outline=(90,90,90,255), width=1)
    draw.line((lx+13,ly+23,lx+57,ly+23), fill=styles["anna"][1], width=3)
    draw.text((lx+67,ly+15), "tracks_anna / anna_2", fill=(20,20,20,255), font=small)
    draw.line((lx+13,ly+53,lx+57,ly+53), fill=styles["mike"][1], width=3)
    draw.text((lx+67,ly+45), "tracks_mike", fill=(20,20,20,255), font=small)
    draw.text((left+7,height-bottom-24), f"{len(records)} unique records; {duplicates} exact duplicates omitted",
              fill=(70,70,70,255), font=small)
    output = base / "storm_outlines_20170812_160332.png"
    image.save(output, dpi=(200,200))
    print(output); print(f"unique_records={len(records)} exact_duplicates_omitted={duplicates}")

if __name__ == "__main__": main()
