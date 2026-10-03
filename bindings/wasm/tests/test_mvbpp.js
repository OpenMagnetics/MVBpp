/**
 * mvbpp WASM binding tests — unified 7-function API.
 *
 * Run with:  node bindings/wasm/tests/test_mvbpp.js
 * Requires:  bindings/wasm/dist/mvbpp.js  (built via `npm run build`)
 *
 * ABT #631: dist/ is the SHIPPED module, not necessarily the one you just built —
 * a green run here says nothing about build-wasm/ until it has been copied over.
 * Point the suite at a fresh build before copying:
 *     MVBPP_WASM_JS=build-wasm/mvbpp.js node tests/test_mvbpp.js
 */

'use strict';

const fs   = require('fs');
const path = require('path');

const DIST_JS   = process.env.MVBPP_WASM_JS
    ? path.resolve(process.env.MVBPP_WASM_JS)
    : path.resolve(__dirname, '../dist/mvbpp.js');
const TESTS_DIR = __dirname;

// ── tiny test harness ─────────────────────────────────────────────────────────

let passed = 0, failed = 0;

function test(name, fn) {
    try {
        fn();
        console.log(`  ✓ ${name}`);
        passed++;
    } catch (e) {
        console.error(`  ✗ ${name}`);
        console.error(`    ${e.message}`);
        failed++;
    }
}

function assert(cond, msg) {
    if (!cond) throw new Error(msg || 'Assertion failed');
}

function assertEqual(a, b, msg) {
    if (a !== b) throw new Error(msg || `Expected ${b}, got ${a}`);
}

// ── helpers ───────────────────────────────────────────────────────────────────

function loadMagneticJson(filename) {
    const raw  = fs.readFileSync(path.join(TESTS_DIR, filename), 'utf8');
    const data = JSON.parse(raw);
    return JSON.stringify(data.magnetic);
}

function isStepBytes(buf) {
    const header = Buffer.from(buf.slice(0, 12));
    return header.toString('ascii') === 'ISO-10303-21';
}

function bufferToString(buf) {
    return Buffer.from(buf).toString('latin1');
}

// Default parameter pack for drawXxx (all 9 trailing args)
const DEF_3D = ['3D', 'XY', 0.0, 'step', 1.0, 32, 'none', '+X+Y+Z'];
const DEF_2D = ['2D', 'XY', 0.0, 'step', 1.0, 32, 'none', ''];
const NO_SIDE = ['3D', 'XY', 0.0, 'step', 1.0, 32, 'none', ''];

// ── main ──────────────────────────────────────────────────────────────────────

async function main() {
    if (!fs.existsSync(DIST_JS)) {
        console.error(`ERROR: ${DIST_JS} not found. Run 'npm run build' in bindings/wasm/ first.`);
        process.exit(1);
    }

    const createMvbpp = require(DIST_JS);
    const mvbpp = await createMvbpp();

    // drawMagnetic/drawTurns/drawWinding (+ *ToPath) gained a trailing
    // `paintCoating` argument: true → OUTER (insulation) footprint, false →
    // CONDUCTING (copper) footprint for FEM winding-loss meshing. embind
    // enforces exact arity, so every call must now pass it. The frontend
    // default is coating (true); wrap these so the existing positional calls
    // below keep exercising that default. rawDraw* keep the unwrapped handles
    // for the explicit copper-vs-coating test.
    const rawDrawMagnetic = mvbpp.drawMagnetic.bind(mvbpp);
    // Bind BEFORE the wrappers below, or "raw" would still carry the appended paintCoating.
    const rawDrawTurns = mvbpp.drawTurns.bind(mvbpp);
    for (const fn of ['drawMagnetic', 'drawMagneticToPath',
                      'drawTurns', 'drawTurnsToPath',
                      'drawWinding', 'drawWindingToPath']) {
        const orig = mvbpp[fn].bind(mvbpp);
        mvbpp[fn] = (...a) => orig(...a, /*paintCoating=*/true);
    }
    // drawMagnetic / drawTurns (+ *ToPath) ALSO take two trailing real-winding arguments
    // (useRealWindingGeometry, femReady). They are documented as optional, but embind enforces
    // exact arity, so they must be PASSED -- undefined selects the documented default (per-turn
    // closed loops, fast drawing compound). Without this every such test fails with
    // "called with 10 arguments, expected 12", which looks like a broken binding rather than a
    // stale test. rawDrawTurns (bound above) keeps the unwrapped handle for the
    // real-winding test.
    for (const fn of ['drawMagnetic', 'drawMagneticToPath',
                      'drawTurns', 'drawTurnsToPath']) {
        const orig = mvbpp[fn].bind(mvbpp);
        mvbpp[fn] = (...a) => orig(...a, /*useRealWindingGeometry=*/undefined,
                                        /*femReady=*/undefined);
    }

    // concentric_basic.json is the BY-NAME fixture (bobbin "basic", shape "E 19/8/5"),
    // the shape of design the frontend actually sends. It used to die in the WASM build
    // with CORE_SHAPE_NOT_FOUND: EI 101/50 (ABT #631) and was swapped out for the
    // fully-inline fixture, which hid the whole by-name path from these tests. Keep BOTH:
    // the by-name one is the regression guard.
    const basicJson    = loadMagneticJson('concentric_rectangular_column_one_turn.json');
    const byNameJson   = loadMagneticJson('concentric_basic.json');
    const etdJson      = loadMagneticJson('ETD49_N87_10uH_5T.json');

    // ── drawMagnetic ─────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawMagnetic');
    console.log('─'.repeat(50));

    test('drawMagnetic returns Uint8Array for basic fixture', () => {
        const result = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        assert(result instanceof Uint8Array, 'Expected Uint8Array');
    });

    test('drawMagnetic output is non-empty', () => {
        const result = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        assert(result.length > 1000, `Too small: ${result.length} bytes`);
    });

    test('drawMagnetic output has valid STEP header', () => {
        const result = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        assert(isStepBytes(result), 'Does not start with ISO-10303-21');
    });

    test('drawMagnetic STEP contains solid geometry markers', () => {
        const text = bufferToString(mvbpp.drawMagnetic(basicJson, ...NO_SIDE));
        assert(
            text.includes('CLOSED_SHELL') || text.includes('ADVANCED_BREP'),
            'No solid geometry markers found in STEP output'
        );
    });

    // ABT #631: a design whose bobbin AND core shape are given BY NAME must build.
    // This is what the frontend sends, and what the whole browser build died on.
    test('drawMagnetic builds a by-name design (bobbin "basic", shape by name) — ABT #631', () => {
        const result = mvbpp.drawMagnetic(byNameJson, ...NO_SIDE);
        assert(result instanceof Uint8Array, 'Expected Uint8Array');
        assert(result.length > 1000, `Too small: ${result.length} bytes`);
        assert(isStepBytes(result), 'Does not start with ISO-10303-21');
    });

    test('ETD49 produces valid STEP', () => {
        const result = mvbpp.drawMagnetic(etdJson, ...NO_SIDE);
        assert(isStepBytes(result), 'ETD49 STEP header invalid');
        assert(result.length > 10_000, `ETD49 output too small: ${result.length}`);
    });

    test('ETD49 produces more geometry than basic', () => {
        const basic = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        const etd   = mvbpp.drawMagnetic(etdJson, ...NO_SIDE);
        assert(etd.length > basic.length,
            `ETD49 (${etd.length}) should be larger than basic (${basic.length})`);
    });

    // ── paintCoating (ABT #7) ────────────────────────────────────────────────

    console.log('\nmvbpp WASM — paintCoating (conductor vs insulation)');
    console.log('─'.repeat(50));

    const STL_3D = ['3D', 'XY', 0.0, 'stl', 1.0, 16, 'none', ''];

    test('paintCoating=false draws copper, differs from coating', () => {
        // rawDrawMagnetic is the UNWRAPPED handle, so it needs the full 12-argument signature:
        // the two trailing real-winding args (useRealWindingGeometry, femReady) must be passed
        // even though undefined selects their defaults -- embind enforces exact arity.
        const coating = rawDrawMagnetic(etdJson, ...STL_3D, true,  undefined, undefined);
        const copper  = rawDrawMagnetic(etdJson, ...STL_3D, false, undefined, undefined);
        assert(coating.length > 0 && copper.length > 0, 'both outputs non-empty');
        assert(coating.length !== copper.length,
            `copper geometry (${copper.length}) must differ from coating (${coating.length})`);
    });

    test('paintCoating=undefined defaults to coating (frontend default)', () => {
        const coating = rawDrawMagnetic(etdJson, ...STL_3D, true,      undefined, undefined);
        const omitted = rawDrawMagnetic(etdJson, ...STL_3D, undefined, undefined, undefined);
        assert(coating.length === omitted.length,
            'omitted paintCoating must equal explicit coating (true)');
    });

    // ── drawMagneticToPath ───────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawMagneticToPath');
    console.log('─'.repeat(50));

    test('writes file to emscripten virtual FS', () => {
        const out  = '/tmp/test_basic.step';
        mvbpp.drawMagneticToPath(basicJson, out, ...NO_SIDE);
        const data = mvbpp.FS.readFile(out);
        assert(data.length > 1000, `Virtual FS file too small: ${data.length}`);
        mvbpp.FS.unlink(out);
    });

    test('virtual FS file has valid STEP header', () => {
        const out = '/tmp/test_header.step';
        mvbpp.drawMagneticToPath(basicJson, out, ...NO_SIDE);
        const data = mvbpp.FS.readFile(out);
        assert(isStepBytes(data), 'Virtual FS file does not start with ISO-10303-21');
        mvbpp.FS.unlink(out);
    });

    test('virtual FS file matches buffer output', () => {
        const out = '/tmp/test_match.step';
        mvbpp.drawMagneticToPath(basicJson, out, ...NO_SIDE);
        const fileData   = mvbpp.FS.readFile(out);
        const bufferData = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        mvbpp.FS.unlink(out);
        const ratio = Math.abs(fileData.length - bufferData.length) /
                      Math.max(fileData.length, bufferData.length);
        assert(ratio < 0.01, `File vs buffer size differs by ${(ratio * 100).toFixed(1)}%`);
    });

    // ── drawCore ─────────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawCore');
    console.log('─'.repeat(50));

    test('drawCore returns valid STEP for enriched core', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const coreJson = JSON.stringify(j.core);
        const result = mvbpp.drawCore(coreJson, ...NO_SIDE);
        assert(isStepBytes(result), 'drawCore STEP header invalid');
        assert(result.length > 500, `drawCore output too small: ${result.length}`);
    });

    test('drawCore with mode=2D returns planar faces', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const coreJson = JSON.stringify(j.core);
        const result = mvbpp.drawCore(coreJson, ...DEF_2D);
        assert(isStepBytes(result), 'drawCore 2D STEP header invalid');
    });

    // ── drawCorePiece ────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawCorePiece');
    console.log('─'.repeat(50));

    test('drawCorePiece returns valid STEP', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const shape = j.core.functionalDescription.shape;
        const shapeJson = JSON.stringify(shape);
        // side="" because a single piece may be centered
        const result = mvbpp.drawCorePiece(shapeJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'none', '');
        assert(isStepBytes(result), 'drawCorePiece STEP header invalid');
    });

    // ── drawBobbin ───────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawBobbin');
    console.log('─'.repeat(50));

    test('drawBobbin returns valid STEP', () => {
        // A core WITH a real bobbin: non-zero wall + column thickness, so
        // buildBobbin renders a hollow body and flanges. (basicJson's enriched
        // bobbin is a zero-thickness column footprint — i.e. no real bobbin — so
        // it is intentionally NOT used here; see the no-bobbin test below.)
        const realBobbin = JSON.stringify({
            processedDescription: {
                columnShape: 'rectangular', columnWidth: 0.005, columnDepth: 0.005,
                columnThickness: 0.001, wallThickness: 0.001, coordinates: [0, 0, 0],
                windingWindows: [{
                    coordinates: [0.0075, 0, 0], height: 0.01, width: 0.005,
                    shape: 'rectangular'
                }]
            }
        });
        const result = mvbpp.drawBobbin(realBobbin, ...NO_SIDE);
        assert(isStepBytes(result), 'drawBobbin STEP header invalid');
    });

    test('drawBobbin on a no-bobbin (planar) core returns empty, not an error', () => {
        // Planar / PCB windings have no physical bobbin: MKF emits a column
        // footprint with columnThickness = wallThickness = 0, so buildBobbin
        // yields an empty shape. drawBobbin must return an empty buffer rather
        // than throw "[mvbpp] unknown C++ exception".
        const planarBobbin = JSON.stringify({
            processedDescription: {
                columnDepth: 0.004505, columnShape: 'oblong', columnThickness: 0,
                columnWidth: 0.00196, coordinates: [0, 0, 0], wallThickness: 0,
                windingWindows: [{
                    area: 2.7e-05, coordinates: [0.00421, 0, 0], height: 0.006,
                    sectionsAlignment: 'centered', sectionsOrientation: 'contiguous',
                    shape: 'rectangular', width: 0.0045
                }]
            }
        });
        let result, threw = false;
        try { result = mvbpp.drawBobbin(planarBobbin, ...NO_SIDE); }
        catch (e) { threw = true; }
        assert(!threw, 'drawBobbin must not throw on a no-bobbin (planar) core');
        assert(result.length === 0, 'no-bobbin core should yield empty bobbin output');
    });

    // ── drawTurns ────────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawTurns');
    console.log('─'.repeat(50));

    test('drawTurns returns valid STEP', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const turnsJson = JSON.stringify(j.coil.turnsDescription);
        const result = mvbpp.drawTurns(turnsJson, ...NO_SIDE);
        assert(isStepBytes(result), 'drawTurns STEP header invalid');
    });

    // drawTurns gained the real-winding form so a viewer that draws core, bobbin and
    // turns as SEPARATE meshes can show the real copper — the single-call drawMagnetic
    // assembly is not an option there, since each component needs its own colour and
    // visibility toggle.
    test('drawTurns builds real-winding conductors, and they differ from the loops', () => {
        const ideal = rawDrawTurns(etdJson, ...NO_SIDE, true, false, false);
        const real  = rawDrawTurns(etdJson, ...NO_SIDE, true, true,  false);
        assert(isStepBytes(ideal), 'ideal drawTurns STEP header invalid');
        assert(isStepBytes(real), 'real-winding drawTurns STEP header invalid');
        assert(real.length > 1000, `real-winding turns too small: ${real.length}`);
        // One continuous body per (winding, parallel) vs one closed loop per turn is a
        // different topology; identical bytes would mean the flag did nothing.
        assert(real.length !== ideal.length,
               'real-winding turns are byte-identical to the idealised ones');
    });

    test('drawTurns real winding refuses a bare turns array instead of silently idealising', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const turnsJson = JSON.stringify(JSON.parse(enriched).coil.turnsDescription);
        let threw = null;
        try { rawDrawTurns(turnsJson, ...NO_SIDE, true, true, false); }
        catch (e) { threw = String(e.message || e); }
        assert(threw !== null, 'expected a refusal for a bare turnsDescription array');
        assert(/useRealWindingGeometry/.test(threw),
               `error should name the flag, got: ${threw}`);
    });

    // ── drawWinding ──────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawWinding');
    console.log('─'.repeat(50));

    test('drawWinding returns valid STEP', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const coilJson = JSON.stringify(j.coil);
        const result = mvbpp.drawWinding(coilJson, 'primary', ...NO_SIDE);
        assert(isStepBytes(result), 'drawWinding STEP header invalid');
    });

    test('drawWinding throws for unknown winding name', () => {
        const enriched = mvbpp._enrichMagnetic(basicJson);
        const j = JSON.parse(enriched);
        const coilJson = JSON.stringify(j.coil);
        let threw = false;
        try {
            mvbpp.drawWinding(coilJson, 'nonexistent', ...NO_SIDE);
        } catch (e) {
            threw = true;
        }
        assert(threw, 'Expected exception for unknown winding name');
    });

    // ── drawView ─────────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — drawView');
    console.log('─'.repeat(50));

    test('drawView returns SVG string', () => {
        const result = mvbpp.drawView(basicJson, true, 'XY', 0.0, 800, 'svg');
        assert(typeof result === 'string', 'Expected string');
        assert(result.startsWith('<?xml'), 'Expected XML/SVG header');
        assert(result.includes('<svg'), 'Expected <svg tag');
    });

    test('drawView without dimensions returns plain SVG', () => {
        const result = mvbpp.drawView(basicJson, false, 'XY', 0.0, 800, 'svg');
        assert(typeof result === 'string', 'Expected string');
        assert(result.includes('<svg'), 'Expected <svg tag');
    });

    test('drawViewToPath writes SVG to virtual FS', () => {
        const out = '/tmp/test_view.svg';
        mvbpp.drawViewToPath(basicJson, out, true, 'XY', 0.0, 800, 'svg');
        const data = mvbpp.FS.readFile(out, { encoding: 'utf8' });
        assert(data.includes('<svg'), 'Expected <svg tag in file');
        mvbpp.FS.unlink(out);
    });

    // ── symmetry / side ──────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — symmetry & side');
    console.log('─'.repeat(50));

    test('symmetry=half reduces output size', () => {
        const full = mvbpp.drawMagnetic(basicJson, ...NO_SIDE);
        const half = mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'half', '+X+Y+Z');
        assert(half.length < full.length,
            `Half symmetry (${half.length}) should be smaller than full (${full.length})`);
    });

    test('symmetry=quarter reduces more than half', () => {
        const half    = mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'half', '+X+Y+Z');
        const quarter = mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'quarter', '+X+Y+Z');
        assert(quarter.length < half.length,
            `Quarter symmetry (${quarter.length}) should be smaller than half (${half.length})`);
    });

    test('side=+X filters geometry', () => {
        const unfiltered = mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'none', '');
        const filtered   = mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'none', '+X');
        assert(filtered.length <= unfiltered.length,
            `Filtered (${filtered.length}) should not be larger than unfiltered (${unfiltered.length})`);
    });

    test('invalid symmetry throws', () => {
        let threw = false;
        try {
            mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'sextant', '+X+Y+Z');
        } catch (e) {
            threw = true;
        }
        assert(threw, 'Expected exception for invalid symmetry');
    });

    test('invalid side throws', () => {
        let threw = false;
        try {
            mvbpp.drawMagnetic(basicJson, '3D', 'XY', 0.0, 'step', 1.0, 32, 'none', '+Q');
        } catch (e) {
            threw = true;
        }
        assert(threw, 'Expected exception for invalid side');
    });

    // ── mesher ────────────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — mesher');
    console.log('─'.repeat(50));

    test('meshBuildRevision reports the linked gmsh fork and MMG', () => {
        const rev = JSON.parse(mvbpp.meshBuildRevision());
        // gmsh answers this itself at run time: the fork's om-4.15.2 build, not a checkout.
        // The commit is abbreviated as git chooses (7 in a fresh clone, more in a big repo).
        assert(/^4\.15\.2-git-[0-9a-f]{7,40}$/.test(rev.gmsh), `gmsh version: ${rev.gmsh}`);
        assert(!/\bMetis\b/.test(rev.gmshBuildOptions), `gmsh built with Metis: ${rev.gmshBuildOptions}`);
        assert(rev.mmg === '5.8.0', `mmg version: ${rev.mmg}`);
        for (const k of ['mvbpp', 'mkf', 'mas'])
            assert(/^[0-9a-f]{40}(-dirty)?$/.test(rev[k]), `${k} revision: ${rev[k]}`);
    });

    // ── error handling ────────────────────────────────────────────────────────

    console.log('\nmvbpp WASM — error handling');
    console.log('─'.repeat(50));

    test('invalid JSON throws', () => {
        let threw = false;
        try { mvbpp.drawMagnetic('not json', ...NO_SIDE); } catch { threw = true; }
        assert(threw, 'Expected exception for invalid JSON');
    });

    test('empty object throws', () => {
        let threw = false;
        try { mvbpp.drawMagnetic('{}', ...NO_SIDE); } catch { threw = true; }
        assert(threw, 'Expected exception for empty object');
    });

    test('2D mode with STL format throws', () => {
        let threw = false;
        try {
            mvbpp.drawMagnetic(basicJson, '2D', 'XY', 0.0, 'stl', 1.0, 32, 'none', '');
        } catch (e) {
            threw = true;
        }
        assert(threw, 'Expected exception for 2D+STL');
    });

    // ── summary ───────────────────────────────────────────────────────────────

    console.log('\n' + '─'.repeat(50));
    console.log(`${passed + failed} tests: ${passed} passed, ${failed} failed`);

    if (failed > 0) process.exit(1);
}

main().catch(e => { console.error(e); process.exit(1); });
