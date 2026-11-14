#!/usr/bin/env python3
"""Write a small self-contained 3D page for whatever target the pipeline ran on.

One function, `write(outdir, data)`. No target name is baked into the HTML: the
page reads its own title, numbers and colour scale from `page.json`, so the same
file works for HER2, EGFR or anything else.

The colouring is the measurement: each residue of the antigen is coloured by how
many DISTINCT binders come within the atlas near-cutoff of it. Blue is none, red
is the most. A popular antibody solved six times cannot vote six times.
"""
import json, os

HTML = r"""<!doctype html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Binder density</title>
<style>
  :root{
    --bg:#FBFBF9; --ink:#151513; --dim:#6B6B63; --line:#DEDED6;
    --panel:#FFFFFF; --accent:#F84016;
  }
  @media (prefers-color-scheme:dark){
    :root:not([data-theme="light"]){
      --bg:#151513; --ink:#F2F2EC; --dim:#9B9B93; --line:#333330;
      --panel:#1E1E1B;
    }
  }
  :root[data-theme="dark"]{
    --bg:#151513; --ink:#F2F2EC; --dim:#9B9B93; --line:#333330; --panel:#1E1E1B;
  }
  *{box-sizing:border-box}
  [hidden]{display:none!important}
  html,body{margin:0;height:100%}
  body{background:var(--bg);color:var(--ink);
    font:15px/1.5 -apple-system,BlinkMacSystemFont,"Segoe UI",Helvetica,Arial,sans-serif;
    display:flex;flex-direction:column}
  header{padding:14px 16px 10px;border-bottom:1px solid var(--line)}
  h1{margin:0;font-size:19px;letter-spacing:-.01em}
  .sub{color:var(--dim);font-size:13px;margin-top:3px}
  main{flex:1;display:flex;min-height:0}
  #stage{flex:1;position:relative;min-width:0;min-height:340px}
  #gl{position:absolute;inset:0}
  aside{width:290px;flex:none;border-left:1px solid var(--line);
    padding:14px 16px;overflow:auto;background:var(--panel)}
  .grp{margin-bottom:16px}
  .lab{font-size:11px;letter-spacing:.08em;text-transform:uppercase;
    color:var(--dim);margin-bottom:6px}
  button{font:inherit;font-size:13px;padding:5px 10px;margin:0 4px 4px 0;
    border:1px solid var(--line);border-radius:6px;background:transparent;
    color:var(--ink);cursor:pointer}
  button[aria-pressed="true"]{background:var(--ink);color:var(--bg);
    border-color:var(--ink)}
  table{border-collapse:collapse;width:100%;font-size:13px}
  td{padding:3px 0;vertical-align:top}
  td:first-child{color:var(--dim);padding-right:10px;white-space:nowrap}
  .num{font-variant-numeric:tabular-nums}
  .bar{height:11px;border-radius:3px;margin:5px 0 3px;
    background:linear-gradient(90deg,#2B5FB0,#8EC4E8,#F2EFD8,#F0A860,#C22A14)}
  .ends{display:flex;justify-content:space-between;color:var(--dim);font-size:11px}
  ol{margin:0;padding-left:18px}
  li{margin-bottom:3px}
  .note{color:var(--dim);font-size:12px;margin-top:4px}
  #load{position:absolute;inset:0;display:flex;align-items:center;
    justify-content:center;color:var(--dim);font-size:13px}
  @media (max-width:860px){
    main{flex-direction:column}
    aside{width:auto;border-left:none;border-top:1px solid var(--line)}
    #stage{min-height:300px}
  }
</style>
</head>
<body>
<header>
  <h1 id="ttl">Binder density</h1>
  <div class="sub" id="sub"></div>
</header>
<main>
  <div id="stage"><div id="gl"></div><div id="load">loading the structure …</div></div>
  <aside>
    <div class="grp">
      <div class="lab">Antigen</div>
      <button id="bSurf" aria-pressed="true">Surface</button>
      <button id="bCart" aria-pressed="false">Backbone</button>
    </div>
    <div class="grp">
      <div class="lab">Binders</div>
      <button id="bBind" aria-pressed="false">Show all</button>
      <button id="bPatch" aria-pressed="false">Zoom to patch</button>
    </div>
    <div class="grp">
      <div class="lab">Distinct binders within reach</div>
      <div class="bar"></div>
      <div class="ends"><span>0</span><span id="dmax">—</span></div>
      <div class="note">One antibody solved six times counts once.</div>
    </div>
    <div class="grp">
      <div class="lab">What was measured</div>
      <table id="facts"></table>
    </div>
    <div class="grp">
      <div class="lab">Most-targeted patch</div>
      <div id="patch" class="num"></div>
      <ol id="who"></ol>
    </div>
  </aside>
</main>
<script src="https://cdn.jsdelivr.net/npm/3dmol@2.0.4/build/3Dmol-min.js"></script>
<script>
const $ = s => document.querySelector(s);

function ramp(t){                 // blue -> red, through a pale middle
  const stops = [[43,95,176],[142,196,232],[242,239,216],[240,168,96],[194,42,20]];
  t = Math.max(0, Math.min(1, t));
  const x = t * (stops.length - 1), i = Math.min(stops.length - 2, Math.floor(x));
  const f = x - i, a = stops[i], b = stops[i + 1];
  const c = k => Math.round(a[k] + (b[k] - a[k]) * f);
  return (c(0) << 16) | (c(1) << 8) | c(2);
}

async function main(){
  const [meta, antigen, binders] = await Promise.all([
    fetch('page.json').then(r => r.json()),
    fetch('antigen.pdb').then(r => r.text()),
    fetch('binders.pdb').then(r => r.text()).catch(() => '')
  ]);

  document.title = meta.title;
  $('#ttl').textContent = meta.title;
  $('#sub').textContent = meta.subtitle;
  const dens = meta.density || {};
  const dmax = Math.max(1, ...Object.values(dens));
  $('#dmax').textContent = dmax;
  $('#facts').innerHTML = meta.facts
    .map(([k, v]) => `<tr><td>${k}</td><td class="num">${v}</td></tr>`).join('');
  if (meta.patch){
    $('#patch').textContent = meta.patch.label;
    $('#who').innerHTML = meta.patch.binders
      .map(n => `<li>${n}</li>`).join('');
  }

  const v = $3Dmol.createViewer('gl', {backgroundColor: 'white'});
  const dark = matchMedia('(prefers-color-scheme: dark)').matches;
  v.setBackgroundColor(dark ? 0x151513 : 0xFBFBF9);

  const ag = v.addModel(antigen, 'pdb');
  let bd = null;
  if (binders.trim()) bd = v.addModel(binders, 'pdb');

  const colorByDensity = atom => ramp((dens[atom.resi] || 0) / dmax);

  function draw(mode, showBinders){
    v.removeAllSurfaces();
    v.setStyle({model: ag}, {});
    if (mode === 'surface'){
      v.setStyle({model: ag}, {cartoon: {color: 0xBBBBB4, opacity: 0.25}});
      v.addSurface($3Dmol.SurfaceType.MS,
                   {opacity: 0.95, colorfunc: colorByDensity}, {model: ag});
    } else {
      v.setStyle({model: ag}, {cartoon: {colorfunc: colorByDensity, thickness: 0.9}});
    }
    // style:"trace" is what draws a CA-only model. Without it the binders
    // load, take up memory and render nothing at all.
    if (bd) v.setStyle({model: bd}, showBinders
      ? {cartoon: {color: 0x8A8A82, style: 'trace', thickness: 0.5, opacity: 0.9}}
      : {});
    v.render();
  }

  let mode = 'surface', show = false;
  const sync = () => {
    $('#bSurf').setAttribute('aria-pressed', mode === 'surface');
    $('#bCart').setAttribute('aria-pressed', mode === 'cartoon');
    $('#bBind').setAttribute('aria-pressed', String(show));
    draw(mode, show);
  };
  $('#bSurf').onclick = () => { mode = 'surface'; sync(); };
  $('#bCart').onclick = () => { mode = 'cartoon'; sync(); };
  $('#bBind').onclick = () => { show = !show; sync(); };
  $('#bPatch').onclick = () => {
    if (!meta.patch) return;
    v.zoomTo({model: ag, resi: meta.patch.resi});
    v.zoom(0.3);           // zoomTo alone puts the camera inside the surface
    v.render();
  };

  sync();
  v.zoomTo({model: ag});
  v.render();
  $('#load').hidden = true;
}
main().catch(e => { $('#load').textContent = 'could not load: ' + e.message; });
</script>
</body>
</html>
"""


def write(outdir, data, antigen_pdb, binders_pdb=None):
    """outdir gets index.html, page.json, antigen.pdb and (if given) binders.pdb."""
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(outdir, 'index.html'), 'w') as fh:
        fh.write(HTML)
    with open(os.path.join(outdir, 'page.json'), 'w') as fh:
        json.dump(data, fh)
    with open(os.path.join(outdir, 'antigen.pdb'), 'w') as fh:
        fh.write(open(antigen_pdb).read())
    if binders_pdb and os.path.exists(binders_pdb):
        with open(os.path.join(outdir, 'binders.pdb'), 'w') as fh:
            fh.write(open(binders_pdb).read())
    return os.path.join(outdir, 'index.html')
