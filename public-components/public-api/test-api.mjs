import {test} from 'node:test';
import assert from 'node:assert/strict';
import worker from './worker.mjs';
const request = (path, init) => worker.fetch(new Request('https://public.example'+path, init));
test('health exposes actual release counts and immutable source provenance',async()=>{
 const r=await request('/api/v1/health');assert.equal(r.status,200);const d=await r.json();
 assert.deepEqual(d.counts,{compounds:53,foods:2,food_reports:172});
 assert.match(d.provenance.source_manifest_sha256,/^[a-f0-9]{64}$/);
 assert.equal(d.mode,'read_only_public_projection');
});
test('compound identity retains complete stereo/form key and native provenance',async()=>{
 const d=await (await request('/api/v1/compounds?limit=1')).json();assert.equal(d.total,53);assert.equal(d.items.length,1);
 const c=d.items[0];assert.match(c.full_inchikey,/^[A-Z]{14}-[A-Z]{10}-[A-Z]$/);assert.ok(c.source.field_provenance.inchi);
 const match=await (await request('/api/v1/compounds?q='+c.full_inchikey)).json();assert.equal(match.items[0].id,c.id);
});
test('food reports preserve native quantity basis and disclaim molecule occurrence',async()=>{
 const f=await (await request('/api/v1/foods')).json();assert.equal(f.total,2);
 const id=f.items[0].id;const r=await request('/api/v1/foods/'+encodeURIComponent(id)+'/reports?limit=100');assert.equal(r.status,200);
 const d=await r.json();assert.ok(d.total>0);for(const row of d.items){assert.equal(row.food_id,id);assert.equal(row.quantity.molecule_presence_asserted,false);assert.equal(row.quantity.basis.unit,'g');assert.ok(row.source);}
 assert.equal((await request('/api/v1/foods/'+encodeURIComponent(id))).status,200);
});
for(const query of ['limit=0','limit=101','limit=true','limit=1.5','offset=-1','offset=01','limit=1&limit=2','root=/private','q='+ 'x'.repeat(201)]){
 test('invalid query rejected: '+query,async()=>{const r=await request('/api/v1/compounds?'+query);assert.equal(r.status,400);assert.equal((await r.json()).error.code,'invalid_request');});
}
test('unknown identifier and route return structured 404',async()=>{
 for(const path of ['/api/v1/foods/absent','/api/v1/jobs','/api/v1/files/private']){const r=await request(path);assert.equal(r.status,404);assert.equal((await r.json()).error.code,'not_found');}
});
test('scientific worker mutations are unavailable',async()=>{
 for(const method of ['POST','PUT','PATCH','DELETE']){const r=await request('/api/v1/jobs',{method});assert.equal(r.status,405);assert.equal(r.headers.get('Allow'),'GET, HEAD');}
});
test('HEAD verifies route with no body',async()=>{const r=await request('/api/v1/health',{method:'HEAD'});assert.equal(r.status,200);assert.equal(await r.text(),'');});
test('raw exports preserve approved counts and source snapshot hash',async()=>{
 const c=await(await request('/api/v1/exports/compounds')).json();assert.equal(c.compounds.length,53);assert.equal(c.foods.length,0);
 const f=await(await request('/api/v1/exports/foods')).json();assert.equal(f.foods.length,2);assert.equal(f.reports.length,172);assert.match(f.source_dataset_sha256,/^[a-f0-9]{64}$/);
});
test('component API distinguishes code from unverified model inference',async()=>{
 const d=await(await request('/api/v1/components')).json();assert.equal(d.components.length,7);
 for(const c of d.components){assert.equal(c.public_runtime,'unverified');assert.ok(c.component);}
});
test('native nutrient and measure references are defined and separately typed',async()=>{
 const d=await(await request('/api/v1/exports/foods')).json();assert.equal(d.measures.length,2);
 const ids=new Map([...d.nutrients,...d.measures].map(x=>[x.id,x.entity_type]));
 for(const r of d.reports){assert.equal(ids.get(r.target_id),r.target_entity_type);}
});
test('OpenAPI contract describes working public routes without an invented server',async()=>{
 const r=await request('/api/v1/openapi.json');assert.equal(r.status,200);const d=await r.json();
 assert.equal(d.openapi,'3.1.0');assert.deepEqual(d.servers,[]);
 for(const p of ['/api/v1/health','/api/v1/foods/{id}/reports','/api/v1/exports/foods'])assert.ok(d.paths[p].get.responses['200']);
});
