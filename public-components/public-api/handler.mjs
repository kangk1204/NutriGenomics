// Public projection only. No research queues, model downloads, or file-system access.
export function createHandler({compounds, foods, components, provenance, openapi}) {
  function reply(request, status, payload, extra = {}) {
    const headers = {'Content-Type':'application/json; charset=utf-8',
      'Cache-Control':'public, max-age=60', 'X-Content-Type-Options':'nosniff', ...extra};
    return new Response(request.method === 'HEAD' ? null : JSON.stringify(payload), {status, headers});
  }
  function error(request, status, code, message, extra) {
    return reply(request, status, {schema_version:'nutri-public-api-v1', error:{code,message}}, extra);
  }
  function page(url, allowed = ['limit','offset','q']) {
    for (const key of url.searchParams.keys()) {
      if (!allowed.includes(key) || url.searchParams.getAll(key).length !== 1) throw new Error('Unsupported or repeated query parameter');
    }
    const integer = (name, fallback, low, high) => {
      const value = url.searchParams.get(name);
      if (value === null) return fallback;
      if (!/^(0|[1-9][0-9]*)$/.test(value)) throw new Error('Invalid '+name);
      const n = Number(value);
      if (!Number.isSafeInteger(n) || n < low || n > high) throw new Error('Invalid '+name);
      return n;
    };
    const q=(url.searchParams.get('q') || '').trim();
    if (q.length > 200 || /[\u0000-\u001f]/.test(q)) throw new Error('Invalid search text');
    return {limit:integer('limit',20,1,100),offset:integer('offset',0,0,1000000),q:q.toLowerCase()};
  }
  function listing(items, paging, release) {
    return {schema_version:'nutri-public-api-v1',release_id:release,total:items.length,
      offset:paging.offset,limit:paging.limit,items:items.slice(paging.offset,paging.offset+paging.limit)};
  }
  return async request => {
    if (!['GET','HEAD'].includes(request.method)) return error(request,405,'method_not_allowed','Only GET and HEAD are available',{Allow:'GET, HEAD'});
    let url;
    try { url=new URL(request.url); } catch { return error(request,400,'invalid_request','Invalid URL'); }
    const path=url.pathname.replace(/\/$/,'') || '/';
    try {
      if (path==='/api/v1/health') {page(url,[]);return reply(request,200,{schema_version:'nutri-public-api-v1',status:'ok',mode:'read_only_public_projection',counts:{compounds:compounds.compounds.length,foods:foods.foods.length,food_reports:foods.reports.length},provenance});}
      if (path==='/api/v1/components') {page(url,[]);return reply(request,200,{schema_version:'nutri-public-api-v1',components});}
      if (path==='/api/v1/provenance') {page(url,[]);return reply(request,200,provenance);}
      if (path==='/api/v1/openapi.json') {page(url,[]);return reply(request,200,openapi);}
      if (path==='/api/v1/compounds') {
        const paging=page(url); const items=compounds.compounds.filter(c=>!paging.q || [c.id,c.full_inchikey,c.iupac_name].some(x=>x.toLowerCase().includes(paging.q)));
        return reply(request,200,listing(items,paging,compounds.release_id));
      }
      if (path==='/api/v1/foods') {
        const paging=page(url); const items=foods.foods.filter(c=>!paging.q || [c.id,c.name_en].some(x=>x.toLowerCase().includes(paging.q)));
        return reply(request,200,listing(items,paging,foods.release_id));
      }
      const match=path.match(/^\/api\/v1\/foods\/([^/]+)(\/reports)?$/);
      if (match) {
        let id;try{id=decodeURIComponent(match[1]);}catch{return error(request,400,'invalid_request','Invalid food identifier');}
        const item=foods.foods.find(x=>x.id===id);
        if (!item) return error(request,404,'not_found','Food identifier is absent from this release');
        if (match[2]) {const paging=page(url,['limit','offset']);return reply(request,200,listing(foods.reports.filter(x=>x.food_id===id),paging,foods.release_id));}
        page(url,[]);return reply(request,200,{schema_version:'nutri-public-api-v1',release_id:foods.release_id,item});
      }
      if (path==='/api/v1/exports/compounds') {page(url,[]);return reply(request,200,compounds,{'Content-Disposition':'attachment; filename="compounds.json"'});}
      if (path==='/api/v1/exports/foods') {page(url,[]);return reply(request,200,foods,{'Content-Disposition':'attachment; filename="foods-pilot.json"'});}
      return error(request,404,'not_found','Route is unavailable');
    } catch { return error(request,400,'invalid_request','Invalid query parameters'); }
  };
}
