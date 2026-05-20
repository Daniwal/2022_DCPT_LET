#!/usr/bin/env python3
from pathlib import Path
import shutil
import re
import sys
import hashlib
from datetime import datetime

if len(sys.argv) != 2:
    print("Usage: python3 make_fluka_case_from_sh12a.py plan02_field01_geoD_mono")
    sys.exit(1)

case = sys.argv[1]
use_znarrow = True

repo_case = Path("/home/dewh/2022_DCPT_LET/data/sh12a/input") / case
if not repo_case.exists():
    raise SystemExit(f"Missing SH12A case folder: {repo_case}")

work = Path("/home/dewh/runs/let_benchmark") / f"fluka_{case}"
work.mkdir(exist_ok=True)

template = Path("/home/dewh/runs/let_benchmark/fluka_geoD_mono/geoD_sobp_step1_smoke.inp")
source_exe = Path("/home/dewh/au-fluka-tools/myfluka_sobp_pr")
source_f = Path("/home/dewh/au-fluka-tools/fluka_sobp_source/source_sampler.f")

inp_path = work / f"{case}_step1_smoke.inp"

shutil.copy2(template, inp_path)
shutil.copy2(source_exe, work / "myfluka_sobp_pr")
shutil.copy2(source_f, work / "source_sampler.f")
shutil.copy2(repo_case / "sobp.dat", work / "sobp_original.dat")

def strip_comment(line):
    return line.split("#", 1)[0].strip()

def parse_beam_dat(path):
    data = {}
    for raw in path.read_text().splitlines():
        line = strip_comment(raw)
        if not line or line.startswith("*"):
            continue
        parts = line.split()
        key = parts[0]
        vals = parts[1:]
        data[key] = vals
    return data

beam = parse_beam_dat(repo_case / "beam.dat")

def parse_mat_dat(path):
    media = []
    current = None
    for raw in path.read_text().splitlines():
        stripped = raw.strip()
        if not stripped:
            continue

        no_comment = stripped.split("#", 1)[0].strip()
        comment = stripped.split("#", 1)[1].strip() if "#" in stripped else ""

        if no_comment.startswith("MEDIUM"):
            parts = no_comment.split()
            current = {"id": parts[1], "comment": comment, "lines": []}
            media.append(current)
        elif no_comment == "END":
            current = None
        elif current is not None:
            current["lines"].append(stripped)
    return media

mat_media = parse_mat_dat(repo_case / "mat.dat")

mat_manifest_lines = []
for m in mat_media:
    mat_manifest_lines.append(f"MEDIUM {m['id']}: {m['comment']}")
    for line in m["lines"]:
        mat_manifest_lines.append(f"  {line}")
mat_manifest_text = "\n".join(mat_manifest_lines)

def parse_detect_meshes(path):
    meshes = {}
    lines = path.read_text().splitlines()
    i = 0
    while i < len(lines):
        line = lines[i].split("#", 1)[0].strip()
        if line != "Geometry Mesh":
            i += 1
            continue

        mesh = {"name": None}
        i += 1
        while i < len(lines):
            sub = lines[i].split("#", 1)[0].strip()
            if not sub:
                i += 1
                break
            if sub in ("Geometry Mesh", "Filter", "Settings", "Output"):
                break

            parts = sub.split()
            if parts and parts[0] == "Name":
                mesh["name"] = parts[1]
            elif parts and parts[0] in ("X", "Y", "Z") and len(parts) >= 4:
                mesh[parts[0]] = (float(parts[1]), float(parts[2]), int(float(parts[3])))
            i += 1

        if mesh.get("name"):
            meshes[mesh["name"]] = mesh

    return meshes

detect_meshes = parse_detect_meshes(repo_case / "detect.dat")
zn = detect_meshes.get("Z_narrow")
detect_mesh_manifest_lines = []
for name, mesh in detect_meshes.items():
    xs = mesh.get("X", "MISSING")
    ys = mesh.get("Y", "MISSING")
    zs = mesh.get("Z", "MISSING")
    detect_mesh_manifest_lines.append(f"{name}: X={xs}, Y={ys}, Z={zs}")
detect_mesh_manifest_text = "\n".join(detect_mesh_manifest_lines)




with open(repo_case / "sobp.dat") as fin, open(work / "sobpcln", "w") as fout:
    for line in fin:
        parts = line.split()
        if len(parts) in (7, 11) and not parts[0].startswith(("*", "#")):
            fout.write(line)

geo = (repo_case / "geo.dat").read_text().splitlines()
rpps = {}
for line in geo:
    parts = line.split()
    if len(parts) >= 8 and parts[0] == "RPP":
        rpps[parts[1]] = parts

rpp3 = rpps["3"]
rpp4 = rpps["4"]

phantom_zmin = float(rpp3[6])
phantom_zmax = float(rpp3[7])
slab_zmin = float(rpp4[6])
slab_zmax = float(rpp4[7])

s = inp_path.read_text()

s = re.sub(
    r"RPP PHANTOM\s+-15\.0\s+15\.0\s+-15\.0\s+15\.0\s+[-0-9.]+\s+[-0-9.]+",
    f"RPP PHANTOM   -15.0   15.0   -15.0   15.0   {phantom_zmin:7.2f} {phantom_zmax:6.2f}",
    s,
)

s = re.sub(
    r"RPP SLAB\s+-15\.0\s+15\.0\s+-15\.0\s+15\.0\s+[-0-9.]+\s+[-0-9.]+",
    f"RPP SLAB      -15.0   15.0   -15.0   15.0   {slab_zmin:7.2f} {slab_zmax:6.2f}",
    s,
)

material_block = """* --- Materials / Assignments ---
MATERIAL         1.0       0.0       1.0       0.0       0.0     1.008HYDROGEN
MATERIAL         6.0       0.0       1.0       0.0       0.0    12.011CARBON
MATERIAL         7.0       0.0       1.0       0.0       0.0    14.007NITROGEN
MATERIAL         8.0       0.0       1.0       0.0       0.0    15.999OXYGEN
MATERIAL        17.0       0.0       1.0       0.0       0.0    35.453CLMAT
MATERIAL        20.0       0.0       1.0       0.0       0.0    40.078CALCIUM

MATERIAL         0.0       0.0     1.032       0.0       0.0      0.0 SOLWATR
COMPOUND   -0.079994  HYDROGEN -0.672956    CARBON -0.023894  NITROGEN SOLWATR
COMPOUND   -0.198689    OXYGEN -0.001383     CLMAT -0.023084   CALCIUM SOLWATR

ASSIGNMA   BLCKHOLE   BLCKHOLE
ASSIGNMA   VACUUM     VAC
ASSIGNMA   AIR        AIR
ASSIGNMA   SOLWATR    SWF
ASSIGNMA   PMMA       PMMA
"""

s = re.sub(
    r"\*SETUP MATERIAL COMPOSITION\*.*?\*Biasing, closer to zero WHAT\(2\) => higher chance of collision",
    material_block + "\n*Biasing, closer to zero WHAT(2) => higher chance of collision",
    s,
    flags=re.DOTALL,
)

if use_znarrow:
    if zn is None:
        raise SystemExit("detect.dat does not define Z_narrow")
    xmin, xmax, nx = zn["X"]
    ymin, ymax, ny = zn["Y"]
    zmin_score, zmax_score, nz = zn["Z"]
else:
    xmin, xmax, nx = -15.00, 15.00, 60
    ymin, ymax, ny = -15.00, 15.00, 60
    zmin_score, zmax_score, nz = phantom_zmin, phantom_zmax, 205

s = re.sub(
    r"USRBIN\s+10\.0\s+208\.0\s+-21\.0\s+[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+EDEP_S",
    f"USRBIN          10.0     208.0     -21.0      {xmax:.2f}      {ymax:.2f}     {zmax_score:.2f}      EDEP_S",
    s,
)

s = re.sub(
    r"USRBIN\s+-15\.00\s+-15\.00\s+[-0-9.]+\s+60\.0\s+60\.0\s+205\.0&",
    f"USRBIN         {xmin:.2f}     {ymin:.2f}    {zmin_score:.2f}       {float(nx):.1f}       {float(ny):.1f}     {float(nz):.1f}       &",
    s,
    count=1,
)

s = re.sub(
    r"USRBIN\s+10\.0\s+1\.0\s+-22\.0\s+[-0-9.]+\s+[-0-9.]+\s+[-0-9.]+TLEN_P",
    f"USRBIN          10.0       1.0     -22.0      {xmax:.2f}      {ymax:.2f}     {zmax_score:.2f}      TLEN_P",
    s,
)

s = re.sub(
    r"USRBIN\s+-15\.00\s+-15\.00\s+[-0-9.]+\s+60\.0\s+60\.0\s+205\.0&",
    f"USRBIN         {xmin:.2f}     {ymin:.2f}    {zmin_score:.2f}       {float(nx):.1f}       {float(ny):.1f}     {float(nz):.1f}       &",
    s,
    count=1,
)

s = re.sub(r"^START.*", "START          10.0", s, flags=re.MULTILINE)

inp_path.write_text(s)

nstat = beam.get("NSTAT", ["100000"])[0]
prod_path = work / f"{case}_step2_prod.inp"
prod_s = re.sub(r"^START.*", f"START       {float(nstat):.1f}", s, flags=re.MULTILINE)
prod_path.write_text(prod_s)

nspots = sum(1 for _ in open(work / "sobpcln"))

def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for chunk in iter(lambda: f.read(1024 * 1024), b""):
            h.update(chunk)
    return h.hexdigest()

manifest = work / "CASE_MANIFEST.txt"
manifest.write_text(f"""# FLUKA case manifest

Generated UTC: {datetime.utcnow().isoformat(timespec='seconds')}Z
Case: {case}

Source SHIELD-HIT folder:
{repo_case}

Generated FLUKA folder:
{work}

Input files and SHA256:
beam.dat   {sha256(repo_case / 'beam.dat') if (repo_case / 'beam.dat').exists() else 'MISSING'}
geo.dat    {sha256(repo_case / 'geo.dat') if (repo_case / 'geo.dat').exists() else 'MISSING'}
mat.dat    {sha256(repo_case / 'mat.dat') if (repo_case / 'mat.dat').exists() else 'MISSING'}
detect.dat {sha256(repo_case / 'detect.dat') if (repo_case / 'detect.dat').exists() else 'MISSING'}
sobp.dat   {sha256(repo_case / 'sobp.dat') if (repo_case / 'sobp.dat').exists() else 'MISSING'}
sobpcln    {sha256(work / 'sobpcln')}

Beam.dat parsed values:
JPART0    {' '.join(beam.get('JPART0', ['MISSING']))}
TMAX0     {' '.join(beam.get('TMAX0', ['MISSING']))}
BEAMSIGMA {' '.join(beam.get('BEAMSIGMA', ['MISSING']))}
BEAMPOS   {' '.join(beam.get('BEAMPOS', ['MISSING']))}
BEAMSAD   {' '.join(beam.get('BEAMSAD', ['MISSING']))}
USECBEAM  {' '.join(beam.get('USECBEAM', ['MISSING']))}
RNDSEED   {' '.join(beam.get('RNDSEED', ['MISSING']))}
NSTAT     {' '.join(beam.get('NSTAT', ['MISSING']))}

Beam.dat physics switches recorded, not one-to-one FLUKA translated:
NEUTRLCUT {' '.join(beam.get('NEUTRLCUT', ['MISSING']))}
DELTAE    {' '.join(beam.get('DELTAE', ['MISSING']))}
DEMIN     {' '.join(beam.get('DEMIN', ['MISSING']))}
STRAGG    {' '.join(beam.get('STRAGG', ['MISSING']))}
MSCAT     {' '.join(beam.get('MSCAT', ['MISSING']))}
NUCRE     {' '.join(beam.get('NUCRE', ['MISSING']))}

Current FLUKA source interpretation:
JPART0=2 interpreted as PROTON
TMAX0 first value interpreted as nominal energy in MeV/u
BEAMPOS copied through the template/source setup
BEAMSAD handled by myfluka_sobp_pr/source_sampler logic
USECBEAM sobp.dat cleaned to sobpcln
NSTAT recorded; smoke input still uses START 10.0

Source sampler:
Executable: myfluka_sobp_pr
Source executable path: {source_exe}
Source Fortran path: {source_f}

SOBP:
Cleaned file: sobpcln
Accepted column counts: 7 or 11
Number of accepted beamlet rows: {nspots}

Geometry extracted from geo.dat:
RPP 3 phantom/SWF z range: {phantom_zmin} to {phantom_zmax} cm
RPP 4 slab/PMMA z range: {slab_zmin} to {slab_zmax} cm

SHIELD-HIT mat.dat parsed media:
{mat_manifest_text}

Material mapping used in FLUKA:
MEDIUM 1 / AIR          -> AIR
MEDIUM 2 / Solid water  -> SOLWATR, rho 1.032 g/cm3, H/C/N/O/Cl/Ca composition
MEDIUM 3 / PMMA         -> PMMA
MEDIUM 4 / Water        -> WATER reference material for LET_water context, not assigned to geometry here
MEDIUM 5 / Si           -> SILICON reference material context, not assigned to geometry here
VAC                     -> VACUUM
BLCKHOLE                -> BLCKHOLE

SHIELD-HIT detect.dat parsed meshes:
{detect_mesh_manifest_text}

Scoring:
SHIELD-HIT Z_narrow mesh from detect.dat: {zn}
FLUKA USRBIN EDEP_S: Z_narrow from detect.dat, x={xmin}..{xmax} cm, y={ymin}..{ymax} cm, z={zmin_score}..{zmax_score} cm, nz={nz}
FLUKA USRBIN TLEN_P: same mesh
Current LET approximation from postprocessing: LET ~= EDEP/TLEN

Run inputs generated:
Smoke input: {inp_path.name}
START smoke: 10.0
Production input: {prod_path.name}
START production from beam.dat NSTAT: {nstat}

Known limitations:
- This generator still patches a FLUKA template.
- beam.dat, mat.dat, and detect.dat are parsed into the manifest; USRBIN mesh values are now taken from detect.dat Z_narrow, but the overall input is still template-patched.
- Geometry currently extracts RPP 3 and RPP 4 specifically; full zone/material table parsing is still future work.
""")


print(f"Created: {work}")
print(f"Smoke input: {inp_path.name}")
print(f"Prod input:  {prod_path.name}")
print(f"Spots:   {nspots}")
print(f"Phantom z: {phantom_zmin} to {phantom_zmax} cm")
print(f"Slab z:    {slab_zmin} to {slab_zmax} cm")
print()
print("Run smoke test with:")
print(f"cd {work}")
print(f"$FLUPRO/bin/rfluka -N 0 -M 1 -e ./myfluka_sobp_pr {inp_path.name}")
print("\nRun production with:")
print(f"$FLUPRO/bin/rfluka -N 0 -M 1 -e ./myfluka_sobp_pr {prod_path.name}")
