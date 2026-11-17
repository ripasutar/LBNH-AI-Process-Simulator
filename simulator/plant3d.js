/* LBNH plant in 3D (three.js r147, works offline).
   Plant3D.init(el)        build the scene inside element el
   Plant3D.update(r)       new solver result {base,new,rev_new,...}
   Plant3D.solving(bool)   pulse while the solver runs
   Plant3D.tour()          guided camera tour with captions
   Plant3D.view(name)      fly the camera to one unit                    */
window.Plant3D = (() => {
  const C = {gas: 0x7895A8, lng: 0x1E9BD7, bog: 0x7B5CD6, n2: 0x2E9E48, he: 0xDB8500, mr: 0xD6453A, c3: 0xF08A5D};
  let el, renderer, scene, camera, controls, clock, labels = [], streams = {}, spinners = [], glows = [], frost = [];
  let R = null, solvingOn = false, fly = null, tourI = -1, tourTimer = null, autoRot = true, ready = false;

  /* ---------- small builders ---------- */
  const M = (c, o = {}) => new THREE.MeshStandardMaterial(Object.assign({color: c, roughness: .6, metalness: .08}, o));
  function add(mesh, x, y, z, shadow = true){ mesh.position.set(x, y, z); mesh.castShadow = shadow; mesh.receiveShadow = true; scene.add(mesh); return mesh; }
  const box = (w, h, d, c, x, y, z, o) => add(new THREE.Mesh(new THREE.BoxGeometry(w, h, d), M(c, o)), x, y + h / 2, z);
  const cyl = (r, h, c, x, y, z, o, seg = 40) => add(new THREE.Mesh(new THREE.CylinderGeometry(r, r, h, seg), M(c, o)), x, y + h / 2, z);
  function dome(r, c, x, y, z, flat = 1){ const m = add(new THREE.Mesh(new THREE.SphereGeometry(r, 40, 16, 0, Math.PI * 2, 0, Math.PI / 2), M(c)), x, y, z); m.scale.y = flat; return m; }
  function hcyl(r, l, c, x, y, z, axis = 'x'){ const m = add(new THREE.Mesh(new THREE.CylinderGeometry(r, r, l, 28), M(c, {metalness: .45})), x, y, z); if (axis === 'x') { m.rotation.order = 'ZYX'; m.rotation.z = Math.PI / 2; } else m.rotation.x = Math.PI / 2; return m; }

  function compressor(x, z, accent, axis = 'x', scale = 1){
    const s = scale;
    box(9 * s, .5, 4.5 * s, 0xBFC6CC, x, 0, z);                               // skid
    hcyl(1.3 * s, 4.6 * s, accent, axis === 'x' ? x - 1.8 * s : x, 1.9 * s, axis === 'x' ? z : z - 1.8 * s, axis);   // casing
    box(axis === 'x' ? 2.8 * s : 2.4 * s, 2.4 * s, axis === 'x' ? 2.4 * s : 2.8 * s, 0x34495E, axis === 'x' ? x + 2.6 * s : x, .5, axis === 'x' ? z : z + 2.6 * s);   // motor
    const disk = hcyl(.9 * s, .25, 0xF2C14E, axis === 'x' ? x + 0.9 * s : x, 1.9 * s, axis === 'x' ? z : z + 0.9 * s, axis);   // coupling
    const stripe = new THREE.Mesh(new THREE.BoxGeometry(.25, 1.9 * s, .3), M(0x222222)); disk.add(stripe);
    spinners.push({m: disk, axis: axis === 'x' ? 'y' : 'y', speed: 6, sec: null});
    return disk;
  }

  /* ---------- pipes with moving particles ---------- */
  function stream(key, pts, color, r = .45, n = 26, size = .55){
    const path = new THREE.CurvePath();
    for (let i = 0; i < pts.length - 1; i++) path.add(new THREE.LineCurve3(new THREE.Vector3(...pts[i]), new THREE.Vector3(...pts[i + 1])));
    const segs = pts.length * 16;
    const tube = new THREE.Mesh(new THREE.TubeGeometry(path, segs, r, 10, false), M(color, {transparent: true, opacity: .55, roughness: .3}));
    scene.add(tube);
    pts.slice(1, -1).forEach(p => { const j = new THREE.Mesh(new THREE.SphereGeometry(r * 1.05, 12, 8), M(color, {transparent: true, opacity: .55})); j.position.set(...p); scene.add(j); });
    const geo = new THREE.SphereGeometry(size * 1.2, 10, 8);
    const mat = new THREE.MeshStandardMaterial({color, emissive: color, emissiveIntensity: 1.4});
    const inst = new THREE.InstancedMesh(geo, mat, n); scene.add(inst);
    const len = path.getLength();
    streams[key] = {path, inst, n, off: Math.random(), speed: 7 / len, rate: 1, tube};
  }

  /* ---------- labels (HTML on top of the canvas) ---------- */
  function label(key, x, y, z, name, cls){
    const d = document.createElement('div'); d.className = 'lb k-' + key + ' ' + (cls || ''); d.dataset.key = key;
    d.innerHTML = `<div class="ln">${name}</div><div class="lv" id="lv-${key}"></div>`;
    d.onclick = () => view(key);
    el.querySelector('.p3-labels').appendChild(d);
    labels.push({d, p: new THREE.Vector3(x, y, z)});
  }

  /* ---------- camera views ---------- */
  const VIEWS = {
    overview: {p: [-52, 64, 100], t: [2, 2, -8]},
    feed: {p: [-80, 16, 26], t: [-56, 4, 0]},
    liq: {p: [-70, 34, 52], t: [-38, 12, 2]},
    tank: {p: [-28, 28, 46], t: [-5, 8, 0]},
    bog: {p: [4, 14, 30], t: [14, 3, 6]},
    nru: {p: [4, 30, 44], t: [28, 12, -2]},
    hls: {p: [34, 20, 36], t: [52, 6, -2]},
    ship: {p: [-30, 30, 2], t: [6, 4, -44]},
  };
  function view(name, dur = 1.6){
    const v = VIEWS[name] || VIEWS.overview; autoRot = false;
    fly = {p0: camera.position.clone(), t0: controls.target.clone(), p1: new THREE.Vector3(...v.p), t1: new THREE.Vector3(...v.t), k: 0, dur, start: performance.now()};
  }

  /* ---------- tour ---------- */
  const fm = (x, d = 0) => Number(x).toLocaleString('en-US', {minimumFractionDigits: d, maximumFractionDigits: d});
  const TOUR = [
    ['overview', r => `This is my LNG and helium plant, in 3D. Natural gas comes in on the left; LNG, liquid helium and nitrogen leave on the right.`],
    ['feed', r => `<b>Feed:</b> ${fm(r.new.feed)} t/h of natural gas arrives. Only <b>${fm(r.he_wt, 2)} %</b> of it is helium.`],
    ['liq', r => `<b>1 · LNG liquefaction:</b> propane and a mixed refrigerant cool the gas to about −162 °C in the tall main exchanger. The compressors use <b>${fm(r.new.LNGP, 1)} MW</b>, most of the plant's power.`],
    ['tank', r => `<b>2 · LNG tank:</b> ${fm(r.new.lng, 1)} t/h of LNG is stored. A little boils off: <b>${fm(r.new.bog, 1)} t/h</b> of boil-off gas, and the helium collects there (${fm(r.new.bog_he, 2)} % helium).`],
    ['bog', r => `<b>Boil-off gas compressor:</b> sends the boil-off gas to helium recovery. <b>${fm(r.new.BOGC, 2)} MW</b>.`],
    ['nru', r => `<b>3 · Nitrogen rejection:</b> a cold column separates nitrogen from helium. <b>${fm(r.new.NRHR, 1)} MW</b>, and ${fm(r.new.n2, 1)} t/h of nitrogen product.`],
    ['hls', r => `<b>4 · Helium liquefier:</b> the helium is cooled to −269 °C and stored as liquid: <b>${fm(r.new.he)} kg/h</b>, about ${fm(r.new.he * 24 / 1000, 1)} tonnes per day.`],
    ['ship', r => `<b>Products:</b> LNG is shipped by sea. Helium is only ${fm(r.he_wt, 2)} % of the feed, but it brings <b>${fm(r.rev_new.he / r.rev_new.total * 100)} %</b> of the revenue: $${fm(r.rev_new.he)} M per year.`],
  ];
  function tourStep(i){
    tourI = Math.max(0, Math.min(TOUR.length - 1, i));
    const [v, t] = TOUR[tourI]; view(v);
    const cap = el.querySelector('.p3-cap'); cap.classList.add('show');
    cap.querySelector('.cn').textContent = (tourI + 1) + '/' + TOUR.length; cap.querySelector('.ct').innerHTML = t(R);
    labels.forEach(l => l.d.classList.toggle('dim', v !== 'overview' && l.d.dataset.key !== v && !(v === 'ship' && l.d.dataset.key === 'heprod')));
  }
  function tour(){ clearInterval(tourTimer); tourStep(0); tourTimer = setInterval(() => { if (tourI >= TOUR.length - 1) clearInterval(tourTimer); else tourStep(tourI + 1); }, 6500); }
  function stopTour(){ clearInterval(tourTimer); tourI = -1; el.querySelector('.p3-cap').classList.remove('show'); labels.forEach(l => l.d.classList.remove('dim')); }

  /* ---------- build the plant ---------- */
  function build(){
    const sk = document.createElement('canvas'); sk.width = 4; sk.height = 256; const sg = sk.getContext('2d'), gr = sg.createLinearGradient(0, 0, 0, 256);
    gr.addColorStop(0, '#3D8FD6'); gr.addColorStop(.55, '#8CC8F0'); gr.addColorStop(1, '#E6F5FD'); sg.fillStyle = gr; sg.fillRect(0, 0, 4, 256);
    const skyTex = new THREE.CanvasTexture(sk); scene.background = skyTex;
    scene.fog = new THREE.Fog(0xCFE9F8, 190, 380);
    scene.add(new THREE.HemisphereLight(0xFFFFFF, 0x8A9A80, .7));
    const sun = new THREE.DirectionalLight(0xFFF3E0, .75); sun.position.set(-60, 110, 70); sun.castShadow = true;
    sun.shadow.mapSize.set(2048, 2048); Object.assign(sun.shadow.camera, {left: -110, right: 110, top: 90, bottom: -90, far: 320}); scene.add(sun);
    // ground, roads, water
    const g = new THREE.Mesh(new THREE.PlaneGeometry(420, 300), M(0x5DAA4E, {roughness: 1, metalness: 0})); g.rotation.x = -Math.PI / 2; g.receiveShadow = true; g.position.z = 20; scene.add(g);
    const yard = new THREE.Mesh(new THREE.PlaneGeometry(166, 72), M(0xC3CBD1, {roughness: .95, metalness: 0})); yard.rotation.x = -Math.PI / 2; yard.position.set(-4, .03, 0); yard.receiveShadow = true; scene.add(yard);
    const grid = new THREE.GridHelper(166, 33, 0xAEB7BE, 0xB4BDC3); grid.scale.z = 72 / 166; grid.position.set(-4, .05, 0); scene.add(grid);
    const water = new THREE.Mesh(new THREE.PlaneGeometry(300, 90), M(0x2E8BC9, {roughness: .12, metalness: .25})); water.rotation.x = -Math.PI / 2; water.position.set(0, .06, -82); scene.add(water);
    box(300, .4, 4, 0x8E989F, 0, 0, -35);                                      // quay
    box(170, .15, 6, 0x4F5961, -4, 0, 39);                                     // front road
    for (let x = -84; x < 80; x += 8) box(4, .05, .4, 0xF2D04B, x, .16, 39);    // road marking
    for (let i = 0; i < 22; i++){ const x = -84 + i * 8 + (i % 2) * 3, z = 47 + (i % 3) * 3;                // trees
      cyl(.35, 2.2, 0x7A5A3A, x, 0, z); const c = add(new THREE.Mesh(new THREE.ConeGeometry(2.1, 5, 8), M([0x3E9A4E, 0x4FAE5A, 0x2F8A45][i % 3], {roughness: .9, metalness: 0})), x, 4.6, z); }
    // section pads
    [[-58, 0, 20, 26, 0xDDE3E8], [-38, 6, 34, 56, 0xA9D8F2], [-5, 0, 30, 32, 0xC7E6F8], [12, 7, 12, 14, 0xD6CBF5], [26, 0, 22, 40, 0xB5E3C2], [51, 2, 20, 34, 0xF8D99B]]
      .forEach(([x, z, w, d, c]) => box(w, .2, d, c, x, 0, z, {roughness: .9}));

    // --- feed: pipeline, inlet skid, scrubber
    box(6, 3, 4, 0x3C6E91, -64, 0, 0);
    cyl(1.3, 11, 0x8FC3E6, -55, 0, -6); dome(1.3, 0x8FC3E6, -55, 11, -6);
    hcyl(1.1, 5, 0x8FC3E6, -55, 1.6, 5, 'x');
    // --- 1 liquefaction: main exchanger (two bundles), propane kettles, compressors, air coolers
    cyl(3.3, 12, 0xE4E9ED, -38, 0, 0, {metalness: .35}); cyl(2.7, 9, 0xE4E9ED, -38, 12.3, 0, {metalness: .35}); dome(2.7, 0xE4E9ED, -38, 21.3, 0, .55);
    for (let i = 0; i < 3; i++){ const fr = cyl(3.36, .55, 0x1E9BD7, -38, 3 + i * 3.5, 0, {emissive: 0x1E9BD7, emissiveIntensity: .35}); frost.push(fr); }
    cyl(.4, 23, 0x9AA5AD, -34.6, 0, -2.5); cyl(.4, 23, 0x9AA5AD, -41.4, 0, -2.5);      // support legs
    for (let i = 0; i < 3; i++) { hcyl(1.2, 6, 0xF08A5D, -48 + i * 8, 1.6, -11, 'x'); box(.4, 1, .4, 0x8A959C, -50 + i * 8, 0, -11); box(.4, 1, .4, 0x8A959C, -46 + i * 8, 0, -11); }
    compressor(-48, 13, 0xF08A5D); compressor(-33, 13, 0xD6453A);
    // air-cooler rack with fans
    const rack = box(30, .7, 7, 0x4E7FA8, -40, 7, 24); for (const [x, z] of [[-54, 21], [-26, 21], [-54, 27], [-26, 27], [-40, 21], [-40, 27]]) box(.5, 7, .5, 0x4E7FA8, x, 0, z);
    for (let i = 0; i < 6; i++){ const f = hcyl(1.5, .2, 0x52616B, -51 + i * 4.4, 8.2, 24, 'y'); f.rotation.set(0, 0, 0);
      const blade = new THREE.Mesh(new THREE.BoxGeometry(2.8, .1, .45), M(0x2B3940)); f.add(blade); const b2 = blade.clone(); b2.rotation.y = Math.PI / 2; f.add(b2);
      spinners.push({m: f, axis: 'y', speed: 4, sec: 'LNGP'}); }
    // --- 2 LNG tank + jetty
    cyl(11, 13, 0xF6F8FA, -5, 0, 0, {metalness: .1, roughness: .7}); dome(11, 0xF6F8FA, -5, 13, 0, .26);
    cyl(11.15, .7, 0x1E9BD7, -5, 12.2, 0, {emissive: 0x1E9BD7, emissiveIntensity: .35}); cyl(11.15, .5, 0x1E9BD7, -5, 1, 0);
    for (let k = 0; k < 9; k++) box(.6, .6, .6, 0xF2D04B, -5 + 11.3 * Math.cos(-.4 + k * .09), 1.4 + k * 1.3, 11.3 * Math.sin(-.4 + k * .09));   // stairs
    for (let i = 0; i < 5; i++) box(.3, .3, .3, 0x8A959C, -5 + 11 * Math.cos(i), 13, 11 * Math.sin(i));
    box(14, 1, 6, 0x9AA5AD, 6, 0, -38);                                         // jetty head
    // --- LNG ship (moss type)
    const hull = box(64, 4.5, 11, 0xC0392B, 8, 0, -48, {roughness: .5}); box(64, 2.4, 11, 0x1F3A5F, 8, 4.5, -48);
    box(8, 7, 9, 0xF2F4F5, 37, 6.9, -48); box(5, 2, 6, 0x34495E, 37, 13.9, -48);
    for (let i = 0; i < 4; i++){ const sx = -18 + i * 13; dome(5, 0xF3E9D2, sx, 6.9, -48, 1); cyl(5, 1.2, 0xE8DCC2, sx, 5.9, -48); }
    // --- BOG compressor
    compressor(13, 7, 0x7B5CD6, 'x', .8);
    // --- 3 nitrogen rejection: cold box, column, N2 cycle compressor, N2 product tank
    box(8, 22, 8, 0x8FD1A3, 24, 0, -2, {metalness: .3}); box(8.2, .8, 8.2, 0x2E9E48, 24, 21.4, -2); box(8.2, .6, 8.2, 0x2E9E48, 24, 0, -2);
    cyl(1.5, 27, 0xEAF6EE, 31, 0, -2, {metalness: .3}); dome(1.5, 0x2E9E48, 31, 27, -2); for (let k = 0; k < 4; k++) cyl(1.56, .4, 0x2E9E48, 31, 5 + k * 6, -2);
    compressor(24, -16, 0x2E9E48, 'x', .9);
    cyl(3, 8, 0xEAF6EE, 30, 0, 14); dome(3, 0x2E9E48, 30, 8, 14, .4); cyl(3.05, .7, 0x2E9E48, 30, 6.3, 14, {emissive: 0x2E9E48, emissiveIntensity: .35});
    // --- 4 helium liquefier: cold box, compressors, dewar, ISO container
    box(6, 14, 6, 0xF7C873, 48, 0, -4, {metalness: .3}); box(6.2, .8, 6.2, 0xDB8500, 48, 13.4, -4); box(6.2, .6, 6.2, 0xDB8500, 48, 0, -4);
    compressor(46, -17, 0xDB8500, 'x', .7); compressor(58, -17, 0xDB8500, 'x', .7);
    cyl(2.2, 7, 0xEDEFF1, 56, 0, -4, {metalness: .5}); dome(2.2, 0xDB8500, 56, 7, -4, .5); cyl(2.26, .5, 0xDB8500, 56, 3.5, -4);
    const iso = box(12, 2.8, 2.8, 0xF29A1E, 55, 1.2, 11); box(12.6, .5, 3.1, 0x6B4A00, 55, .7, 11);
    for (const x of [49, 61]) box(.4, 3.4, 3.2, 0x6B4A00, x, .7, 11);
    hcyl(1.2, 10.5, 0xF4F6F8, 55, 2.6, 11, 'x');

    // --- streams (pipes with moving particles)
    stream('feed', [[-82, 1.6, 0], [-67, 1.6, 0], [-61, 1.6, 0], [-55, 1.6, -3], [-55, 9, -4.7], [-50, 9, -2], [-41.3, 9, 0]], C.gas, .55, 34, .6);
    stream('lng', [[-38, 22.5, 0], [-38, 24, 0], [-20, 24, 0], [-14, 16, 0], [-10, 14.5, 0]], C.lng, .5, 30);
    stream('export', [[-5, 1.4, -11], [-5, 1.4, -24], [2, 3, -32], [6, 3, -37], [6, 7.5, -42.5]], C.lng, .5, 22);
    stream('bog', [[-5, 16.2, -2], [-5, 19, -2], [7, 19, -2], [9, 5, 7], [9.5, 1.9, 7]], C.bog, .42, 22, .5);
    stream('bog2', [[16.6, 1.9, 7], [19, 5, 7], [19, 5, -2], [20, 5, -2]], C.bog, .42, 12, .5);
    stream('n2', [[28, 3, 2], [30, 3, 8], [30, 3, 11]], C.n2, .4, 10, .5);
    stream('hecrude', [[24, 22, -2], [24, 25, -2], [48, 25, -4], [48, 14, -4]], C.he, .35, 20, .45);
    stream('lhe', [[51, 3, -4], [53.8, 3, -4], [56, 3, -1.8], [56, 3, 6], [55, 3.2, 9.6]], C.he, .35, 14, .45);
    stream('mr', [[-35.5, 2, 2.5], [-35.5, 2, 10.5], [-35, 2, 13], [-31, 4.5, 13], [-31, 10, 6], [-35, 10, 1.6]], C.mr, .45, 22, .5);
    stream('c3', [[-45.5, 2, 13], [-45.5, 1.5, 7], [-40, 1.5, -6], [-40, 1.6, -11], [-50, 1.6, -11]], C.c3, .42, 18, .5);

    // --- labels
    label('feed', -64, 8, 0, 'Natural gas in');
    label('liq', -38, 28, 0, '1 · LNG liquefaction', 'big');
    label('tank', -5, 21, 0, '2 · LNG tank', 'big');
    label('bog', 11, 3.5, 10, 'Boil-off compressor');
    label('nru', 27, 31, -2, '3 · N₂ rejection', 'big');
    label('hls', 50, 18, -4, '4 · Helium liquefier', 'big');
    label('ship', 6, 17, -48, 'LNG to market (ship)');
    label('heprod', 70, 2.5, 15, 'Liquid helium to market', 'he');
    label('n2prod', 31, 15, 15, 'Nitrogen product');
  }

  /* ---------- data in ---------- */
  function update(r){
    R = r; if (!ready) return;
    const B = r.base, N = r.new, RN = r.rev_new;
    const ch = (a, b) => { const p = (b - a) / a * 100; return Math.abs(p) < .05 ? '' : ` <span class="${p > 0 ? 'up' : 'dn'}">${p > 0 ? '+' : ''}${fm(p, 1)} %</span>`; };
    const set = (k, h) => { const e = document.getElementById('lv-' + k); if (e) e.innerHTML = h; };
    set('feed', `${fm(N.feed)} t/h · He ${fm(r.he_wt, 2)} %${ch(B.feed, N.feed)}`);
    set('liq', `<b>${fm(N.LNGP, 1)} MW</b>${ch(B.LNGP, N.LNGP)}`);
    set('tank', `LNG <b>${fm(N.lng, 1)} t/h</b>${ch(B.lng, N.lng)}<br>boil-off ${fm(N.bog, 1)} t/h`);
    set('bog', `${fm(N.BOGC, 2)} MW${ch(B.BOGC, N.BOGC)}`);
    set('nru', `<b>${fm(N.NRHR, 1)} MW</b>${ch(B.NRHR, N.NRHR)}`);
    set('hls', `<b>${fm(N.HLS, 2)} MW</b>${ch(B.HLS, N.HLS)}`);
    set('ship', `${fm(N.lng, 1)} t/h · <b>$${fm(RN.lng)} M/yr</b>`);
    set('heprod', `${fm(N.he)} kg/h · <b>$${fm(RN.he)} M/yr</b>${ch(B.he, N.he)}`);
    set('n2prod', `${fm(N.n2, 1)} t/h · $${fm(RN.n2)} M/yr`);
    const rate = (k, a, b) => { if (streams[k]) streams[k].rate = Math.max(.25, Math.min(3, b / a)); };
    rate('feed', B.feed, N.feed); rate('lng', B.lng, N.lng); rate('export', B.lng, N.lng); rate('bog', B.bog, N.bog); rate('bog2', B.bog, N.bog);
    rate('n2', B.n2, N.n2); rate('hecrude', B.he, N.he); rate('lhe', B.he, N.he); rate('mr', B.LNGP, N.LNGP); rate('c3', B.LNGP, N.LNGP);
    spinners.forEach(s => s.rate = s.sec ? N[s.sec] / B[s.sec] : N.total / B.total);
    if (tourI >= 0) tourStep(tourI);
  }
  function solving(on){ solvingOn = on; el.classList.toggle('solving', on); }

  /* ---------- loop ---------- */
  const tmp = new THREE.Object3D(), v3 = new THREE.Vector3();
  let onScreen = true;
  function frame(){
    requestAnimationFrame(frame);
    if (!el.offsetWidth || !onScreen) { clock.getDelta(); return; }
    const dt = Math.min(clock.getDelta(), .05), t = clock.elapsedTime;
    for (const k in streams){
      const s = streams[k]; s.off = (s.off + dt * s.speed * s.rate * (solvingOn ? 2.5 : 1)) % 1;
      for (let i = 0; i < s.n; i++){ s.path.getPointAt((s.off + i / s.n) % 1, v3); tmp.position.copy(v3); tmp.updateMatrix(); s.inst.setMatrixAt(i, tmp.matrix); }
      s.inst.instanceMatrix.needsUpdate = true;
    }
    spinners.forEach(s => s.m.rotation[s.axis] += dt * s.speed * (s.rate || 1) * (s.axis === 'y' && s.sec ? 1 : 1));
    frost.forEach((f, i) => f.material.emissiveIntensity = .2 + .15 * Math.sin(t * 2 + i));
    if (fly){
      fly.k = Math.min(1, (performance.now() - fly.start) / 1000 / fly.dur); const e = fly.k < .5 ? 4 * fly.k ** 3 : 1 - Math.pow(-2 * fly.k + 2, 3) / 2;
      camera.position.lerpVectors(fly.p0, fly.p1, e); controls.target.lerpVectors(fly.t0, fly.t1, e); if (fly.k >= 1) fly = null;
    }
    controls.autoRotate = autoRot && !fly && tourI < 0; controls.update();
    renderer.render(scene, camera);
    const W = el.offsetWidth, H = el.offsetHeight;
    labels.forEach(l => { v3.copy(l.p).project(camera); const vis = v3.z < 1;
      l.d.style.display = vis ? '' : 'none'; l.d.style.transform = `translate(-50%,-100%) translate(${(v3.x * .5 + .5) * W}px,${(-v3.y * .5 + .5) * H}px)`; });
  }

  function init(container){
    el = container;
    try {
      renderer = new THREE.WebGLRenderer({antialias: true});
    } catch (e) { el.querySelector('.p3-fallback').style.display = 'block'; return false; }
    renderer.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2)); renderer.shadowMap.enabled = true; renderer.shadowMap.type = THREE.PCFSoftShadowMap;
    renderer.toneMapping = THREE.NoToneMapping;
    el.querySelector('.p3-canvas').appendChild(renderer.domElement);
    scene = new THREE.Scene(); camera = new THREE.PerspectiveCamera(40, 2, .5, 600);
    camera.position.set(...VIEWS.overview.p);
    controls = new THREE.OrbitControls(camera, renderer.domElement); controls.target.set(...VIEWS.overview.t);
    controls.enableDamping = true; controls.dampingFactor = .08; controls.maxPolarAngle = Math.PI * .48; controls.minDistance = 15; controls.maxDistance = 220; controls.autoRotateSpeed = .35;
    controls.addEventListener('start', () => { autoRot = false; fly = null; });
    clock = new THREE.Clock();
    build();
    const resize = () => { const w = el.offsetWidth, h = el.offsetHeight; if (!w) return; renderer.setSize(w, h); camera.aspect = w / h; camera.updateProjectionMatrix(); };
    new ResizeObserver(resize).observe(el); resize();
    el.querySelector('[data-a=overview]').onclick = () => { stopTour(); view('overview'); setTimeout(() => autoRot = true, 1700); };
    el.querySelector('[data-a=tour]').onclick = tour;
    el.querySelector('[data-a=rotate]').onclick = () => { autoRot = !autoRot; fly = null; };
    el.querySelector('.p3-cap .nx').onclick = () => { clearInterval(tourTimer); tourStep(tourI + 1); };
    el.querySelector('.p3-cap .pv').onclick = () => { clearInterval(tourTimer); tourStep(tourI - 1); };
    el.querySelector('.p3-cap .cx').onclick = stopTour;
    try { new IntersectionObserver(es => { onScreen = es[0].isIntersecting; }, {threshold: 0.02}).observe(el); } catch (e) {}
    ready = true; frame(); if (R) update(R);
    return true;
  }
  return {init, update, solving, tour, view, stopTour};
})();
