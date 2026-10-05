#!/usr/bin/env python3
"""Generate the browser from primary pairs, raw attributions and the frozen gene audit."""
import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]

TEMPLATE = '''<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Enzyme-resolved DDI research resource</title>
<style>
*{box-sizing:border-box}body{margin:0;color:#152235;background:#f7f9fc;font:15px system-ui,sans-serif}header{background:#173f70;color:white;padding:25px 4vw}h1{font-size:25px;margin:0 0 8px}header p{margin:0;line-height:1.5}.wrap{max-width:1300px;margin:auto;padding:20px 4vw}.note{padding:15px;background:#fff4d8;border:1px solid #e9d09a;border-radius:8px;line-height:1.6}.stats{display:flex;gap:20px;margin:20px 0}.stats b{font-size:24px}.controls{display:flex;flex-wrap:wrap;gap:10px;margin:20px 0;align-items:center}input,select,button,.button{font:inherit;padding:9px;border:1px solid #c3ceda;border-radius:6px;background:white}input[type=text]{flex:1;min-width:240px}button,.button{cursor:pointer;background:#173f70;color:white}.button{display:inline-block}.tablewrap{overflow:auto;background:white;border:1px solid #d5dce5;border-radius:8px}table{width:100%;border-collapse:collapse;font-size:14px}th,td{padding:11px;text-align:left;border-bottom:1px solid #e5e9ef}th{background:#eaf1f9;cursor:pointer;white-space:nowrap}td.ids{font:12px monospace;color:#57677a}.pill{display:inline-block;padding:3px 8px;border-radius:9px}.Major{background:#184f95;color:white}.Moderate{background:#3987e5;color:white}.Minor{background:#86b6ef}.Unknown,.NotFound{background:#d5d8dc}footer{font-size:13px;color:#59697c;margin:20px 0;line-height:1.6}.status{padding:10px 0;color:#173f70;font-size:13px}button:focus,input:focus,select:focus,.button:focus{outline:3px solid #e2a034}
</style></head><body><header><h1>Enzyme-resolved DDI research resource</h1><p>1,900 drug pairs · source protein and direction annotations · clinical grades joined locally from your own DDInter CSVs</p></header>
<main class="wrap"><div class="note"><b>Research scope.</b> Raw protein labels are preserved from the source. Accession-backed human genes are shown separately: Q4U2R8 resolves to SLC22A6/OAT1, and Q8TCC7 to SLC22A8/OAT3, despite historical SLCO labels. Resolving protein identity does not validate the drug-pair relation. The exploratory OAT1 estimate uses 12 graded pairs, all containing methotrexate (CID 4112); the single-gene sensitivity uses the same 12 pairs. It does not establish a general or drug-adjusted OAT1 effect. Two historical D3 graph candidates provide partial support; exact current-record derivation remains unestablished. DDInter grades are excluded from this public file. <b>Unknown means no recorded grade; NotFound means your loaded files contain no matching pair.</b> Neither means safe.</div>
<div class="stats"><span><b>1,900</b><br>unique pairs</span><span><b>240</b><br>drugs</span><span><b>35</b><br>source protein labels</span></div>
<div class="controls"><input id="q" type="text" aria-label="Search drug name" placeholder="Search drug name, e.g. tamoxifen"><select id="gene" aria-label="Filter accession-backed gene"><option value="">All accession-backed human genes</option>__GENE_OPTIONS__</select><select id="enz" aria-label="Filter raw source label"><option value="">All source enzymes / transporters</option>__OPTIONS__</select><select id="sev" aria-label="Filter grade" disabled><option value="">All grades</option><option>Major</option><option>Moderate</option><option>Minor</option><option>Unknown</option><option>NotFound</option></select><label for="ddf" class="button" tabindex="0">Add DDInter grades…</label><input id="ddf" type="file" accept=".csv" multiple hidden><button id="reset">Reset filters</button><button id="exp">Export filtered CSV</button></div>
<div id="loadStatus" class="status" role="status"></div><div id="count" class="status"></div><div class="tablewrap"><table><thead><tr><th data-k="drug_A">Drug A</th><th data-k="drug_B">Drug B</th><th data-k="official_genes">Accession-backed gene</th><th data-k="enzymes">Raw source protein label</th><th data-k="mechanism">Source mechanism class</th><th data-k="severity">Grade / DDInter identifiers</th></tr></thead><tbody id="tb"></tbody></table></div>
<footer>Local processing only: no external scripts, network requests or uploads. Pair-level grades are subject to DDInter's terms and are supplied by you. Identifier joins use all candidate DDInter IDs and keep the maximum matched grade under Minor &lt; Unknown &lt; Moderate &lt; Major. Consult README.md, sources/sources.csv and workflow.yaml for provenance and reuse boundaries. This research resource is not a prescribing tool.</footer></main>
<script>
const DATA=__DATA__;
const EXPECTED_FINGERPRINT='aa59130abacc7ef80249c419bc693b9d417d91d2dfc489ff25adf17fb2fc7609';
const tb=document.getElementById('tb'),cnt=document.getElementById('count');
let sortK='drug_A',sortDir=1,hasGrades=false;
const sevord={Major:0,Moderate:1,Minor:2,Unknown:3,NotFound:4,'':5};
const ORD={Minor:0,Unknown:1,Moderate:2,Major:3};
function esc(s){return String(s).replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]))}
function filt(){
 const q=document.getElementById('q').value.toLowerCase().trim(),e=document.getElementById('enz').value,g=document.getElementById('gene').value,s=document.getElementById('sev').value;
 const records=DATA.filter(r=>(!q||r.drug_A.toLowerCase().includes(q)||r.drug_B.toLowerCase().includes(q))&&(!g||r.official_genes.split(';').includes(g))&&(!e||r.enzymes.split(',').map(x=>x.trim()).includes(e))&&(!s||r.severity===s));
 records.sort((a,b)=>{let x=a[sortK],y=b[sortK];if(sortK==='severity'){x=sevord[x];y=sevord[y]}return(x>y?1:x<y?-1:0)*sortDir});
 tb.innerHTML=records.map(r=>`<tr><td>${esc(r.drug_A)}</td><td>${esc(r.drug_B)}</td><td>${esc(r.official_genes)}</td><td>${esc(r.enzymes)}</td><td>${esc(r.mechanism)}</td><td>${r.severity?`<span class="pill ${r.severity}">${esc(r.severity)}</span>`:`<span class="ids">${esc(r.ddinter_ids_A)} · ${esc(r.ddinter_ids_B)}</span>`}</td></tr>`).join('');
 cnt.textContent=records.length+' of '+DATA.length+' pairs';return records;
}
document.querySelectorAll('th').forEach(th=>th.addEventListener('click',()=>{const k=th.dataset.k;if(sortK===k)sortDir*=-1;else{sortK=k;sortDir=1}filt()}));
['q','gene','enz','sev'].forEach(id=>document.getElementById(id).addEventListener('input',filt));
document.getElementById('reset').addEventListener('click',()=>{['q','gene','enz','sev'].forEach(id=>document.getElementById(id).value='');filt()});
function exportCSV(records){const hdr=Object.keys(DATA[0]).filter(k=>k!=='severity'||hasGrades);return[hdr.join(',')].concat(records.map(r=>hdr.map(k=>'"'+String(r[k]).replace(/"/g,'""')+'"').join(','))).join('\\n')}
document.getElementById('exp').addEventListener('click',()=>{const a=document.createElement('a'),url=URL.createObjectURL(new Blob([exportCSV(filt())],{type:'text/csv;charset=utf-8'}));a.href=url;a.download='ddi_enzyme_filtered.csv';a.click();setTimeout(()=>URL.revokeObjectURL(url),0)});
function parseCSV(t){const out=[];let row=[],f='',quoted=false;for(let i=0;i<t.length;i++){const c=t[i];if(quoted){if(c==='"'){if(t[i+1]==='"'){f+='"';i++}else quoted=false}else f+=c}else if(c==='"')quoted=true;else if(c===','){row.push(f);f=''}else if(c==='\\n'||c==='\\r'){if(c==='\\r'&&t[i+1]==='\\n')i++;row.push(f);f='';if(row.length>1||row[0])out.push(row);row=[]}else f+=c}if(quoted)throw new Error('Unclosed quoted CSV field');if(f||row.length){row.push(f);out.push(row)}return out}
function joinGrades(files){const grades=new Map();for(const text of files){const rows=parseCSV(text);if(!rows.length)continue;const h=rows[0].map(x=>x.trim().toLowerCase().replace(/^\\uFEFF/,'')),ia=h.indexOf('ddinterid_a'),ib=h.indexOf('ddinterid_b'),il=h.indexOf('level');if(ia<0||ib<0||il<0)throw new Error('Expected DDInterID_A, DDInterID_B and Level columns');for(const r of rows.slice(1)){const level=(r[il]||'').trim();if(!(level in ORD))continue;const a=(r[ia]||'').trim(),b=(r[ib]||'').trim();if(!a||!b)throw new Error('Missing DDInter pair identifier');const key=[a,b].sort().join('|');if(!grades.has(key)||ORD[level]>ORD[grades.get(key)])grades.set(key,level)}}const result=DATA.map(d=>{let best='';for(const a of d.ddinter_ids_A.split(';'))for(const b of d.ddinter_ids_B.split(';')){const level=grades.get([a.trim(),b.trim()].sort().join('|'));if(level&&(!best||ORD[level]>ORD[best]))best=level}return best||'NotFound'});return result}
async function fingerprint(grades){if(!globalThis.crypto?.subtle)return null;const rows=DATA.map((d,i)=>{const cid=[+d.CID_A,+d.CID_B].sort((a,b)=>a-b);return[cid[0],cid[1],grades[i]]}).sort((a,b)=>a[0]-b[0]||a[1]-b[1]);const canonical=rows.map(r=>r.join(',')+'\\n').join('');const digest=await crypto.subtle.digest('SHA-256',new TextEncoder().encode(canonical));return Array.from(new Uint8Array(digest),x=>x.toString(16).padStart(2,'0')).join('')}
document.getElementById('ddf').addEventListener('change',async event=>{const status=document.getElementById('loadStatus');try{const texts=await Promise.all(Array.from(event.target.files,f=>f.text()));if(!texts.length)return;const grades=joinGrades(texts);const hash=await fingerprint(grades);DATA.forEach((d,i)=>d.severity=grades[i]);hasGrades=true;document.getElementById('sev').disabled=false;const counts={};grades.forEach(g=>counts[g]=(counts[g]||0)+1);status.textContent='Loaded locally: '+Object.entries(counts).map(([g,n])=>g+' '+n).join('; ')+'. '+(hash===EXPECTED_FINGERPRINT?'Exact reference fingerprint match.':hash?'Fingerprint differs from reference; download may be incomplete or a different snapshot.':'Fingerprint check unavailable in this browser context; verify with rebuild_severity.py.');filt()}catch(error){status.textContent='Import failed: '+error.message}});
filt();
</script></body></html>'''


def main():
    import html
    with (ROOT / 'data' / 'ddi_enzyme_database.csv').open(encoding='utf-8', newline='') as fh:
        data = list(csv.DictReader(fh))
    with (ROOT / 'data' / 'protein_annotation_audit.csv').open(encoding='utf-8', newline='') as fh:
        mapping = {r['source_uniprot']:r['official_primary_gene'] for r in csv.DictReader(fh)
                   if r['annotation_resolution']=='verified_unique_human_gene'}
    with (ROOT / 'data' / 'enzyme_pair_attribution.csv').open(encoding='utf-8', newline='') as fh:
        pairs = list(csv.DictReader(fh))
    resolved = {}
    key = lambda r: tuple(sorted((int(r['CID_A']), int(r['CID_B']))))
    for row in pairs:
        if row['uniprot'] in mapping:
            resolved.setdefault(key(row), set()).add(mapping[row['uniprot']])
    for record in data:
        record['official_genes'] = ';'.join(sorted(resolved.get(key(record), set())))
        record['severity'] = '' 
    enzymes = sorted({x.strip() for r in data for x in r['enzymes'].split(',') if x.strip()})
    options = ''.join(f'<option value="{html.escape(e, quote=True)}">{html.escape(e)}</option>' for e in enzymes)
    genes = sorted({g for r in data for g in r['official_genes'].split(';') if g})
    gene_options = ''.join(f'<option value="{html.escape(g, quote=True)}">{html.escape(g)}</option>' for g in genes)
    result = TEMPLATE.replace('__GENE_OPTIONS__', gene_options).replace('__DATA__', json.dumps(data, separators=(',', ':'))).replace('__OPTIONS__', options)
    (ROOT / 'docs' / 'index.html').write_text(result, encoding='utf-8')
    print(f'Generated self-contained browser with {len(data)} public records.')


if __name__ == '__main__':
    main()
