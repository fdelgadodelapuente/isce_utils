#!/usr/bin/env bash
set -euo pipefail

OUT=png
mkdir -p "$OUT"

OUTSIZE="20%"   # 15% casi siempre <200 KB
ZLEVEL=9
UNW_CLIP=20

CMU="$OUT/_unw.txt"
cat > "$CMU" <<'EOF'
0    0   0 130
64   0 180 255
128 255 255 255
192 255 180   0
255 130   0   0
EOF

CMW="$OUT/_wrap.txt"
cat > "$CMW" <<'EOF'
0    255   0   0
42   255 255   0
85     0 255   0
128    0 255 255
170    0   0 255
213  255   0 255
255  255   0   0
EOF

for unw in */*unw; do
  [ -f "$unw" ] || continue
  d="$(basename "$(dirname "$unw")")"

  # 1) fase chica (acelera todo)
  gdal_translate -q "$unw" "$OUT/.${d}_ph.tif" \
    -b 2 -ot Float32 -outsize "$OUTSIZE" "$OUTSIZE" \
    -a_ullr 0 1 1 0

  # 2) unwrapped -> byte -> color (pane izquierdo: x=[0,1])
  gdal_translate -q "$OUT/.${d}_ph.tif" "$OUT/.${d}_unw.tif" \
    -ot Byte -scale -$UNW_CLIP $UNW_CLIP 0 255 \
    -a_ullr 0 1 1 0
  gdaldem color-relief "$OUT/.${d}_unw.tif" "$CMW" "$OUT/.${d}_unw_c.tif"

  # 3) wrapped (desde el float chico) -> byte -> color (pane derecho: x=[1,2])
  gdal_calc.py --quiet -A "$OUT/.${d}_ph.tif" \
    --calc="(((A+3.141592653589793)%(2*3.141592653589793))-3.141592653589793)" \
    --type=Float32 --overwrite --outfile="$OUT/.${d}_wrapf.tif"

  gdal_translate -q "$OUT/.${d}_wrapf.tif" "$OUT/.${d}_wrap.tif" \
    -ot Byte -scale -3.141592653589793 3.141592653589793 0 255 \
    -a_ullr 1 1 2 0
  gdaldem color-relief "$OUT/.${d}_wrap.tif" "$CMW" "$OUT/.${d}_wrap_c.tif"

  # 4) VRT mosaico -> 1 PNG con 2 paneles (izq/der)
  gdalbuildvrt -q "$OUT/.${d}_both.vrt" "$OUT/.${d}_unw_c.tif" "$OUT/.${d}_wrap_c.tif"
  gdal_translate -q "$OUT/.${d}_both.vrt" "$OUT/${d}.png" \
    -of PNG -co ZLEVEL=$ZLEVEL

  rm -f "$OUT"/.${d}_*
done

