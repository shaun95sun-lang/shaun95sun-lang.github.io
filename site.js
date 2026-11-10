/* Shared behaviour. Every block checks that its elements exist, so each page only runs what it has. */
(function () {
  var reduce = window.matchMedia('(prefers-reduced-motion: reduce)').matches;
  var params = new URLSearchParams(location.search);

  /* ---------- name: split into letters so they can type in one by one ---------- */
  var ci = 0;
  document.querySelectorAll('#name .ln').forEach(function (ln) {
    ln.innerHTML = ln.textContent.split('').map(function (c) {
      return '<span class="ch" aria-hidden="true" style="--i:' + (ci++) + '">' + c + '</span>';
    }).join('');
  });

  /* ---------- wind tunnel: NACA 4412, potential flow from a Hess-Smith panel method ---------- */
  (function () {
    var cv = document.getElementById('flow');
    if (!cv) return;
    var ctx = cv.getContext('2d'), view = cv.parentNode;
    var input = document.getElementById('aoa-in'), out = document.getElementById('aoa-out'), stateEl = document.getElementById('state');
    var W = 0, H = 0, aoa = 4, target = 4, colors = {};
    var AMIN = -4, AMAX = 16, STALL = 12, RULER = 46;   // STALL: angle where the animation starts showing separation (a display setting, not a measured value)

    // surface nodes, clockwise: trailing edge -> lower surface -> leading edge -> upper surface -> trailing edge
    var M = 50, NX = [], NY = [];
    function naca(x) {
      var m = 0.04, p = 0.4, t = 0.12;
      var yt = 5 * t * (0.2969 * Math.sqrt(x) - 0.126 * x - 0.3516 * x * x + 0.2843 * x * x * x - 0.1036 * x * x * x * x);
      var yc = x < p ? m / (p * p) * (2 * p * x - x * x) : m / ((1 - p) * (1 - p)) * (1 - 2 * p + 2 * p * x - x * x);
      return [yc + yt, yc - yt];
    }
    for (var i = 0; i <= M; i++) { var x = (1 + Math.cos(Math.PI * i / M)) / 2; NX.push(x); NY.push(naca(x)[1]); }
    for (i = M - 1; i >= 0; i--) { x = (1 + Math.cos(Math.PI * i / M)) / 2; NX.push(x); NY.push(naca(x)[0]); }
    var NP = NX.length - 1, PX = [], PY = [], TX = [], TY = [], PL = [], CX = [], CY = [];
    for (i = 0; i < NP; i++) {
      var dx = NX[i + 1] - NX[i], dy = NY[i + 1] - NY[i], L = Math.hypot(dx, dy);
      PX.push(NX[i]); PY.push(NY[i]); TX.push(dx / L); TY.push(dy / L); PL.push(L);
      CX.push(NX[i] + dx / 2); CY.push(NY[i] + dy / 2);
    }
    // velocity at (px, py) induced by panel j with unit source (s) and unit counter-clockwise vortex (v) strength
    var INF = [0, 0, 0, 0];
    function influence(j, px, py, self) {
      var ex = px - PX[j], ey = py - PY[j], tx = TX[j], ty = TY[j], nx = -ty, ny = tx, L = PL[j];
      var xl = ex * tx + ey * ty, yl = ex * nx + ey * ny;
      var lr = 0.5 * Math.log((xl * xl + yl * yl) / ((xl - L) * (xl - L) + yl * yl)) / (2 * Math.PI);
      var be = (self ? Math.PI : Math.atan2(yl, xl - L) - Math.atan2(yl, xl)) / (2 * Math.PI);
      INF[0] = lr * tx + be * nx; INF[1] = lr * ty + be * ny;      // source
      INF[2] = -be * tx + lr * nx; INF[3] = -be * ty + lr * ny;    // vortex = source rotated 90 deg
      return INF;
    }
    function solve(A, b) {   // Gaussian elimination with partial pivoting
      var n = b.length;
      for (var c = 0; c < n; c++) {
        var p = c;
        for (var r = c + 1; r < n; r++) if (Math.abs(A[r][c]) > Math.abs(A[p][c])) p = r;
        var t = A[c]; A[c] = A[p]; A[p] = t; t = b[c]; b[c] = b[p]; b[p] = t;
        for (r = c + 1; r < n; r++) {
          var f = A[r][c] / A[c][c];
          for (var k = c; k < n; k++) A[r][k] -= f * A[c][k];
          b[r] -= f * b[c];
        }
      }
      var x = new Array(n);
      for (r = n - 1; r >= 0; r--) { var s = b[r]; for (k = r + 1; k < n; k++) s -= A[r][k] * x[k]; x[r] = s / A[r][r]; }
      return x;
    }
    // unknowns: NP source strengths + one vortex strength; rows: zero normal flow at each control point + Kutta condition
    function strengths(ux, uy) {
      var A = [], b = [], kt = new Array(NP + 1).fill(0);
      for (var i = 0; i < NP; i++) {
        var row = new Array(NP + 1).fill(0), nx = -TY[i], ny = TX[i];
        for (var j = 0; j < NP; j++) {
          var v = influence(j, CX[i], CY[i], i === j);
          row[j] = v[0] * nx + v[1] * ny;
          row[NP] += v[2] * nx + v[3] * ny;
          if (i === 0 || i === NP - 1) { kt[j] += v[0] * TX[i] + v[1] * TY[i]; kt[NP] += v[2] * TX[i] + v[3] * TY[i]; }
        }
        A.push(row); b.push(-(ux * nx + uy * ny));
      }
      A.push(kt); b.push(-(ux * (TX[0] + TX[NP - 1]) + uy * (TY[0] + TY[NP - 1])));
      return solve(A, b);
    }
    function inside(px, py) {
      var c = false;
      for (var i = 0, j = NP - 1; i < NP; j = i++) {
        if ((NY[i] > py) !== (NY[j] > py) && px < (NX[j] - NX[i]) * (py - NY[i]) / (NY[j] - NY[i]) + NX[i]) c = !c;
      }
      return c;
    }
    // body-frame velocity grids for free streams (1,0) and (0,1); any angle of attack is cos a * F0 + sin a * F90
    var G = { x0: -1.1, y0: -1.05, h: 0.0125, nx: 0, ny: 0 };
    G.nx = Math.round(3.8 / G.h) + 1; G.ny = Math.round(2.1 / G.h) + 1;
    var F0 = new Float32Array(G.nx * G.ny * 2), F90 = new Float32Array(G.nx * G.ny * 2), SOLID = new Uint8Array(G.nx * G.ny);
    (function () {
      var s0 = strengths(1, 0), s9 = strengths(0, 1);
      for (var gj = 0; gj < G.ny; gj++) for (var gi = 0; gi < G.nx; gi++) {
        var px = G.x0 + gi * G.h, py = G.y0 + gj * G.h, o = gj * G.nx + gi;
        if (px > -0.01 && px < 1.01 && py > -0.06 && py < 0.11 && inside(px, py)) { SOLID[o] = 1; continue; }
        var a = 1, b = 0, c = 0, d = 1;   // u0, v0, u90, v90 start from the free stream
        for (var j = 0; j < NP; j++) {
          var v = influence(j, px, py, false);
          a += s0[j] * v[0] + s0[NP] * v[2]; b += s0[j] * v[1] + s0[NP] * v[3];
          c += s9[j] * v[0] + s9[NP] * v[2]; d += s9[j] * v[1] + s9[NP] * v[3];
        }
        F0[o * 2] = a; F0[o * 2 + 1] = b; F90[o * 2] = c; F90[o * 2 + 1] = d;
      }
    })();
    var VEL = [0, 0];
    // returns false inside the foil; skips solid corners so the surface stays crisp
    function velBody(bx, by, ca, sa) {
      var fx = (bx - G.x0) / G.h, fy = (by - G.y0) / G.h;
      if (fx < 0 || fy < 0 || fx >= G.nx - 1 || fy >= G.ny - 1) { VEL[0] = ca; VEL[1] = sa; return true; }
      var i0 = Math.floor(fx), j0 = Math.floor(fy), tx = fx - i0, ty = fy - j0, u = 0, v = 0, ws = 0;
      for (var c = 0; c < 4; c++) {
        var di = c & 1, dj = c >> 1, o = (j0 + dj) * G.nx + i0 + di;
        if (SOLID[o]) continue;
        var w = (di ? tx : 1 - tx) * (dj ? ty : 1 - ty);
        u += w * (ca * F0[o * 2] + sa * F90[o * 2]); v += w * (ca * F0[o * 2 + 1] + sa * F90[o * 2 + 1]); ws += w;
      }
      if (ws < 0.35) return false;
      VEL[0] = u / ws; VEL[1] = v / ws;
      return true;
    }
    // wind frame: free stream along +x, foil pitched nose-up about the quarter chord
    var CA = 1, SA = 0, B = [0, 0];
    function toBody(wx, wy) { var x = wx - 0.25; B[0] = x * CA - wy * SA + 0.25; B[1] = x * SA + wy * CA; return B; }
    function velWind(wx, wy) {
      toBody(wx, wy);
      if (!velBody(B[0], B[1], CA, SA)) return false;
      var u = VEL[0], v = VEL[1];
      VEL[0] = u * CA + v * SA; VEL[1] = -u * SA + v * CA;
      return true;
    }

    function readColors() {
      var s = getComputedStyle(cv);
      function v(n) { return s.getPropertyValue(n).trim(); }
      colors = { bg: v('--c-bg'), grid: v('--c-grid'), line: v('--c-line'), hot: v('--c-hot'), foil: v('--c-foil'), alert: v('--c-alert'), tick: v('--c-tick') };
    }
    function size() {
      var dpr = Math.min(window.devicePixelRatio || 1, 2), r = cv.getBoundingClientRect();
      W = r.width; H = r.height;
      cv.width = Math.round(W * dpr); cv.height = Math.round(H * dpr);
      ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
      readColors(); ctx.fillStyle = colors.bg; ctx.fillRect(0, 0, W, H);
    }
    function chord() { return Math.min((W - RULER) * 0.44, H * 1.5); }
    function x0() { return RULER + (W - RULER) * 0.28; }
    function sx(wx) { return x0() + wx * chord(); }
    function sy(wy) { return H * 0.5 - wy * chord(); }
    function aoaToY(a) { return H - 12 - (a - AMIN) / (AMAX - AMIN) * (H - 24); }
    function xLeft() { return (RULER - x0()) / chord(); }
    function xRight() { return (W - x0()) / chord(); }
    function yHalf() { return H * 0.5 / chord(); }

    // grid + angle-of-attack ruler, redrawn every frame so the trail fade never builds them up
    function drawGrid() {
      ctx.strokeStyle = colors.grid; ctx.lineWidth = 1;
      ctx.beginPath();
      for (var x = RULER + 0.5; x < W; x += 40) { ctx.moveTo(x, 0); ctx.lineTo(x, H); }
      for (var y = (H / 2 % 40) + 0.5; y < H; y += 40) { ctx.moveTo(RULER, y); ctx.lineTo(W, y); }
      ctx.stroke();
      ctx.fillStyle = colors.bg; ctx.fillRect(0, 0, RULER, H);
      ctx.font = '600 12px "Saira Extra Condensed", "Arial Narrow", sans-serif'; ctx.textBaseline = 'middle';
      for (var a = AMIN; a <= AMAX; a++) {
        var ty = Math.round(aoaToY(a)) + 0.5, major = a % 4 === 0, c = a > STALL ? colors.alert : colors.tick;
        ctx.strokeStyle = c;
        ctx.beginPath(); ctx.moveTo(40, ty); ctx.lineTo(major ? 28 : 34, ty); ctx.stroke();
        if (major) { ctx.fillStyle = c; ctx.fillText((a > 0 ? '+' : '') + a, 6, ty); }
      }
      ctx.strokeStyle = colors.tick; ctx.beginPath(); ctx.moveTo(40.5, 12); ctx.lineTo(40.5, H - 12); ctx.stroke();
      var my = aoaToY(target);
      ctx.fillStyle = colors.hot; ctx.beginPath(); ctx.moveTo(42, my); ctx.lineTo(52, my - 6); ctx.lineTo(52, my + 6); ctx.closePath(); ctx.fill();
    }
    function drawFoil() {
      ctx.beginPath();
      for (var i = 0; i <= NP; i++) {
        var x = NX[i] - 0.25, y = NY[i], wx = x * CA + y * SA + 0.25, wy = -x * SA + y * CA;
        i ? ctx.lineTo(sx(wx), sy(wy)) : ctx.moveTo(sx(wx), sy(wy));
      }
      ctx.closePath();
      ctx.fillStyle = colors.bg; ctx.fill();
      ctx.lineWidth = 1.6; ctx.strokeStyle = colors.foil; ctx.stroke();
    }
    function setAngle(a) { var r = a * Math.PI / 180; CA = Math.cos(r); SA = Math.sin(r); }

    // particles live in the wind frame; above the stall angle the ones over the aft upper surface get a separated wake
    var P = [];
    function spawn(p, anywhere) {
      p.x = anywhere ? xLeft() + Math.random() * (xRight() - xLeft()) : xLeft() - Math.random() * 0.2;
      p.y = (Math.random() * 2 - 1) * yHalf();
      p.j = 0;
      return p;
    }
    function upperGap(wx, wy) {   // height above the upper surface in chord units, or -1 when not over the aft upper surface
      toBody(wx, wy);
      if (B[0] < 0.3 || B[0] > 1.8) return -1;
      var top = B[0] <= 1 ? naca(B[0])[0] : 0;
      var g = B[1] - top;
      return g >= 0 ? g : -1;
    }
    function frame(dt) {
      aoa += (target - aoa) * Math.min(1, dt * 6);
      setAngle(aoa);
      var sep = Math.max(0, aoa - STALL) / 4, k = 0.12 * (W / chord()) * dt;
      ctx.globalAlpha = 0.16; ctx.fillStyle = colors.bg; ctx.fillRect(0, 0, W, H); ctx.globalAlpha = 1;
      drawGrid();
      ctx.lineWidth = 1.3; ctx.lineCap = 'round';
      for (var i = 0; i < P.length; i++) {
        var p = P[i], ox = p.x, oy = p.y;
        if (!velWind(p.x, p.y)) { spawn(p, false); continue; }
        var u = VEL[0], v = VEL[1], sp = Math.hypot(u, v), wake = false;
        // midpoint step keeps particles hugging the leading edge instead of cutting through it
        if (!velWind(p.x + u * k * 0.5, p.y + v * k * 0.5)) { spawn(p, false); continue; }
        p.x += VEL[0] * k; p.y += VEL[1] * k;
        if (sep > 0) {
          var g = upperGap(p.x, p.y);
          if (g >= 0 && g < 0.3) {
            p.j += (Math.random() - 0.5) * 0.02 * sep * dt * 60; p.j *= 0.97;
            p.x -= 0.55 * Math.min(1, sep) * VEL[0] * k * Math.exp(-g / 0.12);   // slow, recirculating air
            p.y += p.j * 0.3 + 0.02 * sep * k;
            wake = true;
          }
        }
        if (ox < xLeft()) continue;
        ctx.strokeStyle = wake ? colors.alert : sp > 1.25 ? colors.hot : colors.line;
        ctx.beginPath(); ctx.moveTo(sx(ox), sy(oy)); ctx.lineTo(sx(p.x), sy(p.y)); ctx.stroke();
        if (p.x > xRight() + 0.05 || Math.abs(p.y) > yHalf() + 0.1) spawn(p, false);
      }
      drawFoil();
    }
    function drawStatic() {   // reduced motion: integrate whole streamlines once
      setAngle(aoa);
      ctx.fillStyle = colors.bg; ctx.fillRect(0, 0, W, H);
      drawGrid();
      ctx.lineWidth = 1;
      for (var y = -yHalf(); y <= yHalf(); y += 0.045) {
        var x = xLeft(), yy = y, n = 0;
        ctx.beginPath(); ctx.moveTo(sx(x), sy(yy));
        while (x < xRight() && n++ < 2000) {
          if (!velWind(x, yy)) break;
          var u = VEL[0], v = VEL[1], s = 0.01 / Math.max(0.2, Math.hypot(u, v));
          if (!velWind(x + u * s * 0.5, yy + v * s * 0.5)) break;
          x += VEL[0] * s; yy += VEL[1] * s;
          ctx.lineTo(sx(x), sy(yy));
        }
        ctx.strokeStyle = colors.line; ctx.stroke();
      }
      drawFoil();
    }

    function setAoa(v) {
      target = Math.max(AMIN, Math.min(AMAX, v));
      input.value = target;
      var s = Math.abs(target).toFixed(1);
      out.textContent = (target < 0 ? '−' : '+') + (s.length < 4 ? '0' + s : s);
      var stall = target > STALL;
      view.classList.toggle('stall', stall);
      stateEl.textContent = stall ? 'SEPARATION' : 'ATTACHED';
      stateEl.classList.toggle('bad', stall);
      if (reduce) { aoa = target; drawStatic(); }
    }
    input.addEventListener('input', function () { setAoa(parseFloat(input.value)); });
    cv.addEventListener('pointermove', function (e) {
      var y = e.clientY - cv.getBoundingClientRect().top;
      setAoa(Math.round((AMIN + (H - 12 - y) / (H - 24) * (AMAX - AMIN)) * 2) / 2);
    });

    var visible = true, last = 0;
    if ('IntersectionObserver' in window) new IntersectionObserver(function (es) { visible = es[0].isIntersecting; }).observe(cv);
    function loop(t) {
      var dt = last ? Math.min((t - last) / 1000, 0.05) : 0.016; last = t;
      if (visible) frame(dt);
      requestAnimationFrame(loop);
    }
    size();
    for (i = 0; i < 280; i++) P.push(spawn({}, true));
    window.addEventListener('resize', function () { size(); if (reduce) drawStatic(); });
    if (document.fonts && document.fonts.ready) document.fonts.ready.then(function () { if (reduce) drawStatic(); });
    var q0 = parseFloat(params.get('aoa'));
    setAoa(isNaN(q0) ? 4 : q0);
    aoa = target;
    if (reduce) drawStatic(); else requestAnimationFrame(loop);
  })();

  /* ---------- shared helpers for the 3D viewers ---------- */
  function fetchBytes(name) {
    return fetch(name).then(function (r) { if (!r.ok) throw new Error(name + ' ' + r.status); return r.text(); }).then(function (b64) {
      var bin = atob(b64.trim()), bytes = new Uint8Array(bin.length);
      for (var k = 0; k < bin.length; k++) bytes[k] = bin.charCodeAt(k);
      return bytes;
    });
  }
  function whenNear(el, fn) {
    if (!('IntersectionObserver' in window)) { fn(); return; }
    var o = new IntersectionObserver(function (es) { if (es[0].isIntersecting) { o.disconnect(); fn(); } }, { rootMargin: '400px 0px' });
    o.observe(el);
  }
  function makeRenderer(canvas, msgEl) {
    if (!window.THREE || !THREE.OrbitControls) { msgEl.textContent = '3D 檢視器載入失敗'; return null; }
    try { var r = new THREE.WebGLRenderer({ canvas: canvas, antialias: true, alpha: true }); }
    catch (e) { msgEl.textContent = '這個瀏覽器不支援 3D 顯示'; return null; }
    r.setPixelRatio(Math.min(window.devicePixelRatio || 1, 2));
    return r;
  }
  function runLoop(canvas, renderer, scene, cam, ctl, extra) {
    function fit() {
      var w = canvas.clientWidth, h = canvas.clientHeight;
      renderer.setSize(w, h, false); cam.aspect = w / h; cam.updateProjectionMatrix();
    }
    fit();
    if ('ResizeObserver' in window) new ResizeObserver(fit).observe(canvas);
    var on = true, lastT = 0;
    new IntersectionObserver(function (es) { on = es[0].isIntersecting; }).observe(canvas);
    (function tick(t) {
      var dt = lastT ? Math.min((t - lastT) / 1000, 0.05) : 0.016; lastT = t || 0;
      if (on) { if (extra) extra(dt); ctl.update(); renderer.render(scene, cam); }
      requestAnimationFrame(tick);
    })(0);
  }
  // colormaps: [position 0..1, r, g, b]. Speed is mapped over 0..1.4 so free stream (1.0) sits at 0.714 as a quiet green
  var CP_STOPS = [[0, 111, 140, 255], [0.25, 70, 84, 154], [0.5, 107, 101, 96], [0.75, 208, 106, 38], [1, 255, 58, 31]];
  var SP_STOPS = [[0, 18, 30, 74], [0.45, 24, 70, 58], [0.714, 36, 120, 52], [0.86, 255, 204, 80], [1, 255, 90, 31]];
  function cmap(stops, t, out) {
    t = Math.min(1, Math.max(0, t));
    var i = 0;
    while (i < stops.length - 2 && t > stops[i + 1][0]) i++;
    var a = stops[i], b = stops[i + 1], f = (t - a[0]) / (b[0] - a[0]);
    out[0] = a[1] + (b[1] - a[1]) * f; out[1] = a[2] + (b[2] - a[2]) * f; out[2] = a[3] + (b[3] - a[3]) * f;
    return out;
  }
  function cssGrad(stops) { return 'linear-gradient(90deg,' + stops.map(function (s) { return 'rgb(' + s.slice(1).join(',') + ') ' + (s[0] * 100) + '%'; }).join(',') + ')'; }

  /* ---------- Hrist CFD viewer: test.STL body, flow from C7 at 15 m/s, alpha 4 deg, props off (OpenFOAM simpleFoam k-omega SST on lena) ---------- */
  (function () {
    var hcv = document.getElementById('hr-3d');
    if (!hcv) return;
    var hload = document.getElementById('hr-load');
    document.getElementById('leg-cp').style.background = cssGrad(CP_STOPS);
    document.getElementById('leg-sp').style.background = cssGrad(SP_STOPS);
    // CFD frame: flow +X, up +Y, span Z (m); nose at x = -0.25
    var STREAM_PTS = 120, STREAM_COUNTS = [41, 41, 41, 41, 41, 41, 41, 9];
    // seed layers: below wing 22 mm, above wing 8 / 35 / 65 / 100 / 150 / 220 mm, cooling inlet; only the ones near the wing are drawn
    var SHOW_LAYER = [true, true, true, false, false, false, false, true];
    // body is test.STL moved into the CFD frame; its Cp is taken from the nearest C7 surface faces (median gap 0.7 mm)
    var FILES = { geom: 'media/hrist-geom.b64.txt', cp: 'media/hrist-cp.b64.txt', stream: 'media/c7-stream-a4.b64.txt' };

    whenNear(hcv, function () {
      var renderer = makeRenderer(hcv, hload);
      if (!renderer) return;
      var scene = new THREE.Scene();
      var cam = new THREE.PerspectiveCamera(33, 1, 0.005, 20);
      scene.add(new THREE.HemisphereLight(0xffffff, 0x332a22, 0.85));
      var sun = new THREE.DirectionalLight(0xffffff, 0.5); sun.position.set(-1, 2, 1.5); scene.add(sun);
      var ctl = new THREE.OrbitControls(cam, hcv);
      ctl.enableDamping = true; ctl.enablePan = false;
      ctl.minDistance = 0.25; ctl.maxDistance = 2.4;
      cam.position.set(-0.6, 0.38, 0.7); ctl.target.set(0.04, 0.02, 0); ctl.update();

      var bodyGeo = new THREE.BufferGeometry();
      scene.add(new THREE.Mesh(bodyGeo, new THREE.MeshLambertMaterial({ vertexColors: true, side: THREE.DoubleSide })));
      var streams = new THREE.Group(); scene.add(streams);
      var partGeo = new THREE.BufferGeometry();
      var parts = new THREE.Points(partGeo, new THREE.PointsMaterial({ size: 0.007, color: 0xffcc50, transparent: true, opacity: 0.95, depthWrite: false }));
      parts.frustumCulled = false; scene.add(parts);
      var a4 = 4 * Math.PI / 180;
      scene.add(new THREE.ArrowHelper(new THREE.Vector3(Math.cos(a4), Math.sin(a4), 0), new THREE.Vector3(-0.46, -0.03, 0), 0.12, 0xff9830, 0.03, 0.018));

      var cp = null, rgb = [0, 0, 0], lines = [], PER = 6, pt = null;
      function paintBody(showCp) {
        var n = cp.length, cols = new Float32Array(n * 9);
        for (var t = 0; t < n; t++) {
          if (showCp) cmap(CP_STOPS, (cp[t] + 1) / 2, rgb); else { rgb[0] = 170; rgb[1] = 166; rgb[2] = 158; }
          for (var k = 0; k < 3; k++) { cols[t * 9 + k * 3] = rgb[0] / 255; cols[t * 9 + k * 3 + 1] = rgb[1] / 255; cols[t * 9 + k * 3 + 2] = rgb[2] / 255; }
        }
        bodyGeo.setAttribute('color', new THREE.BufferAttribute(cols, 3));
      }
      function buildStreams(data) {
        var line = 0;
        STREAM_COUNTS.forEach(function (count, gi) {
          if (!SHOW_LAYER[gi]) { line += count; return; }
          for (var l = 0; l < count; l++, line++) {
            var pos = new Float32Array(STREAM_PTS * 3), col = new Float32Array(STREAM_PTS * 3), pts = [];
            for (var i = 0; i < STREAM_PTS; i++) {
              var o = (line * STREAM_PTS + i) * 4;
              pos[i * 3] = data[o]; pos[i * 3 + 1] = data[o + 1]; pos[i * 3 + 2] = data[o + 2];
              cmap(SP_STOPS, data[o + 3] / 1.4, rgb);
              col[i * 3] = rgb[0] / 255; col[i * 3 + 1] = rgb[1] / 255; col[i * 3 + 2] = rgb[2] / 255;
              pts.push([data[o], data[o + 1], data[o + 2], data[o + 3]]);
            }
            var g = new THREE.BufferGeometry();
            g.setAttribute('position', new THREE.BufferAttribute(pos, 3));
            g.setAttribute('color', new THREE.BufferAttribute(col, 3));
            streams.add(new THREE.Line(g, new THREE.LineBasicMaterial({ vertexColors: true, transparent: true, opacity: 0.7 })));
            lines.push(pts);
          }
        });
        pt = new Float32Array(lines.length * PER).map(function (_, i) { return (i % PER) / PER + Math.random() / PER; });
        partGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(lines.length * PER * 3), 3));
      }
      function stepParticles(dt) {
        if (!pt || !parts.visible) return;
        var pos = partGeo.attributes.position.array, n = STREAM_PTS;
        for (var li = 0; li < lines.length; li++) {
          var ln = lines[li];
          for (var k = 0; k < PER; k++) {
            var idx = li * PER + k, t = pt[idx], i = Math.min(n - 2, Math.floor(t * (n - 1)));
            if (!reduce) t += dt * 0.1 * Math.max(0.1, ln[i][3]);
            if (t >= 1) t -= 1;
            pt[idx] = t;
            var f = t * (n - 1), j = Math.min(n - 2, Math.floor(f)), w = f - j, a = ln[j], b = ln[j + 1];
            pos[idx * 3] = a[0] + (b[0] - a[0]) * w; pos[idx * 3 + 1] = a[1] + (b[1] - a[1]) * w; pos[idx * 3 + 2] = a[2] + (b[2] - a[2]) * w;
          }
        }
        partGeo.attributes.position.needsUpdate = true;
      }

      function toggle(id, fn) {
        var b = document.getElementById(id);
        b.addEventListener('click', function () { var v = b.getAttribute('aria-pressed') !== 'true'; b.setAttribute('aria-pressed', v); fn(v); });
      }
      toggle('hr-cp', function (v) { if (cp) paintBody(v); });
      toggle('hr-st', function (v) { streams.visible = v; parts.visible = v; });

      hload.textContent = 'LOADING · 2.0 MB';
      Promise.all([fetchBytes(FILES.geom), fetchBytes(FILES.cp), fetchBytes(FILES.stream)]).then(function (r) {
        bodyGeo.setAttribute('position', new THREE.BufferAttribute(new Float32Array(r[0].buffer), 3));
        bodyGeo.computeVertexNormals();
        bodyGeo.computeBoundingSphere();   // the first frames ran before the mesh arrived and cached an empty sphere
        cp = new Float32Array(r[1].buffer);
        paintBody(true);
        buildStreams(new Float32Array(r[2].buffer));
        stepParticles(0);
        hload.textContent = '';
      }).catch(function () { hload.textContent = 'CFD 資料載入失敗'; });

      runLoop(hcv, renderer, scene, cam, ctl, stepParticles);
    });
  })();

  /* ---------- Geiravör 3D viewer ---------- */
  (function () {
    var gcv = document.getElementById('gei-3d');
    if (!gcv) return;
    var gload = document.getElementById('gei-load');
    whenNear(gcv, function () {
      var renderer = makeRenderer(gcv, gload);
      if (!renderer) return;
      var scene = new THREE.Scene();
      var cam = new THREE.PerspectiveCamera(32, 1, 1, 20000);
      scene.add(new THREE.HemisphereLight(0xffffff, 0x222222, 0.6));
      var sun = new THREE.DirectionalLight(0xffffff, 0.8); sun.position.set(1, 2, 1.5); scene.add(sun);
      var rim = new THREE.DirectionalLight(0xff9830, 0.7); rim.position.set(-1.5, 0.5, -2); scene.add(rim);
      var mat = new THREE.MeshStandardMaterial({ roughness: 0.75, metalness: 0.05, flatShading: true });
      mat.color.set(getComputedStyle(document.documentElement).getPropertyValue('--model').trim());
      var ctl = new THREE.OrbitControls(cam, gcv);
      ctl.enableDamping = true; ctl.enablePan = false;
      ctl.autoRotate = !reduce; ctl.autoRotateSpeed = 1.1;
      ctl.addEventListener('start', function () { ctl.autoRotate = false; });
      fetchBytes('media/geiravor.b64.txt').then(function (bytes) {
        var buf = bytes.buffer, n = new Uint32Array(buf, 0, 2), g = new THREE.BufferGeometry();
        g.setAttribute('position', new THREE.BufferAttribute(new Float32Array(buf, 8, n[0] * 3), 3));
        g.setIndex(new THREE.BufferAttribute(new Uint32Array(buf, 8 + n[0] * 12, n[1] * 3), 1));
        g.computeBoundingSphere();
        scene.add(new THREE.Mesh(g, mat));
        var r = g.boundingSphere.radius;
        cam.position.set(r * 1.5, r * 0.9, r * 1.9); cam.near = r / 50; cam.far = r * 20; cam.updateProjectionMatrix();
        ctl.minDistance = r * 0.8; ctl.maxDistance = r * 6;
        ctl.target.set(0, 0, 0); ctl.update();
        gload.textContent = '';
      }).catch(function () { gload.textContent = '模型載入失敗'; });
      runLoop(gcv, renderer, scene, cam, ctl);
    });
  })();
})();
