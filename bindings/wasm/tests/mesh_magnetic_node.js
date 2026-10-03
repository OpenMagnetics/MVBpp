// node tests/mesh_magnetic_node.js <mas.json> <reference.msh> <out.msh>
// Meshes the MAS file's magnetic in the WASM module with the recipe the reference run recorded
// (<reference.msh>.recipe.json minus provenance, units and defaults) and writes the msh2 and its
// sidecars, for mvbpp_mesh_diff against a native run. Not part of the quick suite: a real mesh
// takes minutes.
const fs = require('fs');
const [mas, ref, out] = process.argv.slice(2);
const JS = process.env.MVBPP_WASM_JS || require('path').join(__dirname, '..', 'build-wasm', 'mvbpp.js');
const eff = JSON.parse(fs.readFileSync(ref + '.recipe.json', 'utf8'));
const recipe = {};
for (const [sec, v] of Object.entries(eff)) {
    if (sec === 'provenance' || sec === 'recipe_version') continue;
    recipe[sec] = {};
    for (const [k, x] of Object.entries(v))
        if (!k.startsWith('_unit_') && !k.endsWith('_default')) recipe[sec][k] = x;
}
require(JS)().then(m => {
    const t0 = Date.now();
    const r = JSON.parse(m.meshMagnetic(JSON.stringify(JSON.parse(fs.readFileSync(mas, 'utf8')).magnetic),
                                        JSON.stringify(recipe), 'msh2', 'm'));
    fs.writeFileSync(out, r.mesh);
    for (const [suffix, text] of Object.entries(r.sidecars)) fs.writeFileSync(out + suffix, text);
    console.log(`meshed in ${((Date.now() - t0) / 1000).toFixed(1)} s -> ${out}`);
}).catch(e => { console.error('FAILED', e.message || e); process.exit(1); });
