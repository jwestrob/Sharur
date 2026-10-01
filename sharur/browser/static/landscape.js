/* Functional landscape: genomes on a canvas, coloured by clade, with hover, click, box selection and highlight. */
(function () {
  "use strict";
  const wrap = document.getElementById("land");
  if (!wrap) return;
  const canvas = wrap.querySelector("canvas");
  const tip = wrap.querySelector(".land-tip");
  const ctx = canvas.getContext("2d");
  const $ = (sel) => document.querySelector(sel);
  const status = $("[data-land-status]"), search = $("[data-land-search]"), rankSelect = $("[data-land-rank]");
  const kindSelect = $("[data-land-kind]"), matches = $("[data-land-matches]");
  const selBar = $("[data-land-selection]"), selLabel = $("[data-land-selected]"), selLimit = $("[data-land-limit]");
  const css = getComputedStyle(document.documentElement);
  const tok = (name, fallback) => (css.getPropertyValue(name) || fallback).trim();
  const C = { ink2: tok("--ink-2", "#47505c"), ink3: tok("--ink-3", "#7c8591"), line: tok("--line", "#e3e0d8"),
              accent: tok("--accent", "#0f7c6e"), surface: tok("--surface", "#fff"), grey: "#8b939c" };
  const RANKS = ["phylum", "class", "order", "family", "genus"];
  const TOP = 12;

  let D = null, rank = wrap.dataset.color || "", colours = [], colourOf = null, order = [];
  let highlight = null, selection = new Set(), hover = -1, view = null, drag = null;

  kindSelect.addEventListener("change", () => {
    const u = new URL(location.href);
    u.searchParams.set("kind", kindSelect.value);
    location.href = u.toString();
  });
  rankSelect.addEventListener("change", () => {
    rank = rankSelect.value || D.default_rank;
    syncUrl();
    recolour();
    draw();
  });

  function syncUrl() {
    const u = new URL(location.href);
    u.searchParams.set("color", rankSelect.value);
    if (search.value.trim()) u.searchParams.set("q", search.value.trim()); else u.searchParams.delete("q");
    history.replaceState(null, "", u.toString());
  }

  function load() {
    fetch("/api/landscape?kind=" + encodeURIComponent(wrap.dataset.kind)).then((r) => r.json().then((j) => [r.status, j]))
      .then(([code, j]) => {
        if (code === 202) { setTimeout(load, 1200); return; }
        if (j.status !== "ready") { status.textContent = j.message || "The map could not be computed."; return; }
        D = j;
        status.hidden = true;
        rank = rank || D.default_rank;
        if (!rankSelect.value) rankSelect.querySelector('option[value=""]').textContent = "auto (" + D.default_rank + ")";
        recolour();
        axes();
        if (search.value.trim()) applySearch();
        layout();
      }).catch(() => { status.textContent = "The map could not be loaded."; });
  }

  // ---- colours ------------------------------------------------------------
  function recolour() {
    const table = D.ranks[rank];
    const counts = new Map();
    table.idx.forEach((k) => counts.set(k, (counts.get(k) || 0) + 1));
    const named = [...counts.entries()].filter(([k]) => table.names[k] !== "Unclassified").sort((a, b) => b[1] - a[1]);
    colours = new Array(table.names.length).fill(null);
    named.slice(0, TOP).forEach(([k], i) => { colours[k] = D.palette[i % D.palette.length]; });
    colourOf = (i) => colours[table.idx[i]];
    // grey first, then coloured clades from largest to smallest so small ones stay visible
    const rankOf = new Map(named.map(([k], i) => [k, i]));
    order = D.ids.map((_, i) => i).sort((a, b) => {
      const ca = colours[table.idx[a]] ? 1 : 0, cb = colours[table.idx[b]] ? 1 : 0;
      if (ca !== cb) return ca - cb;
      return (rankOf.get(table.idx[a]) ?? 0) - (rankOf.get(table.idx[b]) ?? 0);
    });
    $("[data-land-rank-label]").textContent = "by " + rank;
    legend(named, counts);
  }

  function legend(named, counts) {
    const box = $("[data-land-clades]");
    box.innerHTML = "";
    const table = D.ranks[rank];
    named.slice(0, TOP).forEach(([k, n]) => box.appendChild(chip(table.names[k], n, colours[k])));
    const rest = named.slice(TOP).reduce((s, [, n]) => s + n, 0) + (counts.get(table.names.indexOf("Unclassified")) || 0);
    if (rest) {
      const more = document.createElement("span");
      more.className = "count land-rest";
      more.innerHTML = '<i class="land-swatch" style="background:' + C.grey + '"></i>';
      more.append(document.createTextNode(rest.toLocaleString() + " in other or unclassified " + rank + " groups"));
      box.appendChild(more);
    }
  }

  function chip(name, n, colour) {
    const a = document.createElement("button");
    a.type = "button";
    a.className = "chip land-chip";
    a.innerHTML = '<i style="--dot:' + colour + '"></i>';
    const label = document.createElement("span");
    label.textContent = name;
    const count = document.createElement("span");
    count.className = "count";
    count.textContent = " " + n.toLocaleString();
    a.append(label, count);
    a.addEventListener("click", () => { search.value = search.value === name ? "" : name; applySearch(); });
    return a;
  }

  // ---- axes summary -------------------------------------------------------
  function axes() {
    const box = $("[data-land-axes]");
    box.innerHTML = "";
    const pct = (x) => (100 * x).toFixed(1) + "%";
    const r = (x) => (x === null || x === undefined ? "–" : x.toFixed(2));
    ["Axis 1 (horizontal)", "Axis 2 (vertical)"].forEach((title, i) => {
      const sec = document.createElement("div");
      sec.className = "land-axis";
      const h = document.createElement("h3");
      h.textContent = title + " · " + pct(D.explained[i]) + " of variance";
      const p = document.createElement("p");
      p.className = "count";
      p.textContent = "r = " + r(D.size_r[i]) + " with features per genome" +
        (D.completeness_r[i] === null ? "" : ", r = " + r(D.completeness_r[i]) + " with completeness");
      sec.append(h, p);
      [["high", "→ / ↑"], ["low", "← / ↓"]].forEach(([end, arrow]) => {
        const row = document.createElement("div");
        row.className = "land-load";
        const tag = document.createElement("span");
        tag.className = "count land-end";
        tag.textContent = (i === 0 ? (end === "high" ? "right" : "left") : (end === "high" ? "top" : "bottom"));
        row.appendChild(tag);
        D.loadings[i][end].forEach((f) => {
          const a = document.createElement("a");
          a.className = "chip";
          a.href = D.feature_url.replace("{id}", encodeURIComponent(f.id));
          a.title = f.name;
          a.textContent = f.label;
          row.appendChild(a);
        });
        sec.appendChild(row);
      });
      box.appendChild(sec);
    });
    $("[data-land-method]").textContent = "Principal components of " + D.label + " presence (" +
      D.used_features.toLocaleString() + " features carried by two or more genomes, centred per feature; truncated SVD). " +
      "Axes that follow features per genome partly reflect genome size and annotation depth." +
      (D.left_out ? " " + D.left_out.toLocaleString() + " genomes with no " + D.label + " are left off." : "") +
      (D.note ? " " + D.note : "");
  }

  // ---- layout and drawing -------------------------------------------------
  function layout() {
    const width = Math.max(320, wrap.clientWidth);
    const height = Math.round(Math.min(760, Math.max(380, width * 0.64)));
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(width * dpr);
    canvas.height = Math.round(height * dpr);
    canvas.style.width = width + "px";
    canvas.style.height = height + "px";
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    let x0 = Infinity, x1 = -Infinity, y0 = Infinity, y1 = -Infinity;
    for (let i = 0; i < D.x.length; i++) {
      x0 = Math.min(x0, D.x[i]); x1 = Math.max(x1, D.x[i]);
      y0 = Math.min(y0, D.y[i]); y1 = Math.max(y1, D.y[i]);
    }
    const pad = 28;
    // one scale for both axes, so distances read the same in every direction
    const s = Math.min((width - 2 * pad) / ((x1 - x0) || 1), (height - 2 * pad) / ((y1 - y0) || 1));
    const ox = (width - s * (x1 - x0)) / 2, oy = (height - s * (y1 - y0)) / 2;
    view = { width, height, sx: (x) => ox + (x - x0) * s, sy: (y) => height - (oy + (y - y0) * s) };
    D.px = Float32Array.from(D.x, view.sx);
    D.py = Float32Array.from(D.y, view.sy);
    draw();
  }

  function radius(i) {
    const c = D.completeness[i];
    if (c === null) return 2.6;
    return 1.6 + 2.6 * Math.max(0, Math.min(1, (c - 50) / 50));
  }

  function draw() {
    if (!view) return;
    const { width, height } = view;
    ctx.clearRect(0, 0, width, height);
    // origin crosshair
    ctx.strokeStyle = C.line;
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(view.sx(0), 0); ctx.lineTo(view.sx(0), height);
    ctx.moveTo(0, view.sy(0)); ctx.lineTo(width, view.sy(0));
    ctx.stroke();
    ctx.fillStyle = C.ink3;
    ctx.font = "11px ui-sans-serif, system-ui, sans-serif";
    ctx.fillText("axis 1 · " + (100 * D.explained[0]).toFixed(1) + "%", width - 110, view.sy(0) - 6);
    ctx.save();
    ctx.translate(view.sx(0) + 12, 14);
    ctx.fillText("axis 2 · " + (100 * D.explained[1]).toFixed(1) + "%", 0, 0);
    ctx.restore();

    for (const i of order) dot(i, highlight && !highlight.has(i) ? 0.12 : (highlight ? 0.95 : 0.78));
    if (highlight) for (const i of order) if (highlight.has(i)) dot(i, 1, true);
    ctx.lineWidth = 1.6;
    ctx.strokeStyle = C.accent;
    for (const i of selection) { ctx.beginPath(); ctx.arc(D.px[i], D.py[i], radius(i) + 2.2, 0, 2 * Math.PI); ctx.stroke(); }
    if (hover >= 0) {
      ctx.lineWidth = 2;
      ctx.strokeStyle = C.ink2;
      ctx.beginPath(); ctx.arc(D.px[hover], D.py[hover], radius(hover) + 3.5, 0, 2 * Math.PI); ctx.stroke();
    }
    if (drag && drag.moved) {
      ctx.fillStyle = "rgba(127,127,127,0.10)";
      ctx.strokeStyle = C.accent;
      ctx.lineWidth = 1;
      const x = Math.min(drag.x0, drag.x1), y = Math.min(drag.y0, drag.y1);
      ctx.fillRect(x, y, Math.abs(drag.x1 - drag.x0), Math.abs(drag.y1 - drag.y0));
      ctx.strokeRect(x, y, Math.abs(drag.x1 - drag.x0), Math.abs(drag.y1 - drag.y0));
    }
  }

  function dot(i, alpha, emphasise) {
    const colour = colourOf(i) || C.grey;
    const r = radius(i);
    const hollow = D.completeness[i] !== null && D.completeness[i] < 70;
    ctx.globalAlpha = colourOf(i) ? alpha : alpha * 0.6;
    ctx.beginPath();
    ctx.arc(D.px[i], D.py[i], r, 0, 2 * Math.PI);
    if (hollow) { ctx.strokeStyle = colour; ctx.lineWidth = 1.2; ctx.stroke(); }
    else { ctx.fillStyle = colour; ctx.fill(); }
    if (emphasise) { ctx.strokeStyle = C.surface; ctx.lineWidth = 0.8; ctx.stroke(); }
    ctx.globalAlpha = 1;
  }

  // ---- interaction --------------------------------------------------------
  function nearest(x, y) {
    let best = -1, bd = 100;   // within 10 px
    for (let i = 0; i < D.px.length; i++) {
      if (highlight && !highlight.has(i)) continue;
      const dx = D.px[i] - x, dy = D.py[i] - y, d = dx * dx + dy * dy;
      if (d < bd) { bd = d; best = i; }
    }
    return best;
  }

  function pos(ev) {
    const r = canvas.getBoundingClientRect();
    return [ev.clientX - r.left, ev.clientY - r.top];
  }

  function lineage(i) {
    // GTDB reuses placeholder names across ranks; show each name once
    const names = RANKS.map((r) => D.ranks[r].names[D.ranks[r].idx[i]]).filter((n) => n && n !== "Unclassified");
    return names.filter((n, k) => n !== names[k - 1]).join("; ");
  }

  function showTip(i, x, y) {
    if (i < 0) { tip.hidden = true; return; }
    tip.innerHTML = "<b></b><div class='count l'></div><div class='count c'></div>";
    tip.querySelector("b").textContent = D.ids[i];
    tip.querySelector(".l").textContent = lineage(i);
    const c = D.completeness[i];
    tip.querySelector(".c").textContent = (c === null ? "completeness unknown" : c.toFixed(1) + "% complete") +
      " · " + D.features[i].toLocaleString() + " " + D.label;
    tip.hidden = false;
    const left = x + 16 + 260 > view.width ? x - 270 : x + 16;
    tip.style.left = left + "px";
    tip.style.top = Math.max(0, y - 10) + "px";
  }

  canvas.addEventListener("mousedown", (ev) => {
    if (!D) return;
    const [x, y] = pos(ev);
    drag = { x0: x, y0: y, x1: x, y1: y, moved: false, add: ev.shiftKey };
  });
  window.addEventListener("mousemove", (ev) => {
    if (!D || !view) return;
    const [x, y] = pos(ev);
    if (drag) {
      drag.x1 = x; drag.y1 = y;
      if (Math.abs(x - drag.x0) + Math.abs(y - drag.y0) > 4) drag.moved = true;
      if (drag.moved) { tip.hidden = true; draw(); return; }
    }
    if (ev.target !== canvas) return;
    const i = nearest(x, y);
    if (i !== hover) { hover = i; draw(); }
    canvas.style.cursor = i >= 0 ? "pointer" : "crosshair";
    showTip(i, x, y);
  });
  window.addEventListener("mouseup", (ev) => {
    if (!drag) return;
    const d = drag;
    drag = null;
    const [ux, uy] = pos(ev);   // the release point counts even when no mousemove arrived
    d.x1 = ux; d.y1 = uy;
    if (Math.abs(ux - d.x0) + Math.abs(uy - d.y0) > 4) d.moved = true;
    if (d.moved) {
      const xa = Math.min(d.x0, d.x1), xb = Math.max(d.x0, d.x1), ya = Math.min(d.y0, d.y1), yb = Math.max(d.y0, d.y1);
      if (!d.add) selection = new Set();
      for (let i = 0; i < D.px.length; i++) {
        if (highlight && !highlight.has(i)) continue;
        if (D.px[i] >= xa && D.px[i] <= xb && D.py[i] >= ya && D.py[i] <= yb) selection.add(i);
      }
      updateSelection();
      draw();
    } else if (ev.target === canvas && hover >= 0) {
      const href = "/genome/" + encodeURIComponent(D.ids[hover]);
      if (ev.metaKey || ev.ctrlKey) window.open(href, "_blank"); else location.href = href;
    }
  });
  canvas.addEventListener("mouseleave", () => { if (!drag) { hover = -1; tip.hidden = true; draw(); } });

  function updateSelection() {
    const n = selection.size;
    selBar.hidden = n === 0;
    if (!n) return;
    selLabel.textContent = n.toLocaleString() + " genome" + (n === 1 ? "" : "s") + " selected";
    const ids = [...selection].map((i) => D.ids[i]);
    const inline = "genomes:" + ids.join(",");
    const links = selBar.querySelectorAll("[data-land-open]");
    const point = (token) => links.forEach((a) => {
      a.classList.remove("disabled");
      a.href = a.dataset.landOpen === "matrix"
        ? "/matrix?a=" + encodeURIComponent(token) + "&order=cluster"
        : "/compare?a=" + encodeURIComponent(token);
    });
    selLimit.textContent = "";
    if (inline.length <= D.inline_selection) point(inline);
    else {
      // long selections travel as a short key the server keeps until it restarts
      links.forEach((a) => { a.classList.add("disabled"); a.href = "#"; });
      const asked = ids.join(",");
      fetch("/api/selection", { method: "POST", headers: { "Content-Type": "application/json" },
                                body: JSON.stringify({ genomes: ids }) })
        .then((r) => r.json()).then((j) => { if (j.token && asked === [...selection].map((i) => D.ids[i]).join(",")) point(j.token); })
        .catch(() => { selLimit.textContent = "The selection could not be stored."; });
    }
    // a selection that is mostly one clade can open as that clade, at any size
    const table = D.ranks[rank], tally = new Map();
    for (const i of selection) tally.set(table.idx[i], (tally.get(table.idx[i]) || 0) + 1);
    const [top, count] = [...tally.entries()].sort((a, b) => b[1] - a[1])[0];
    const name = table.names[top];
    if (n >= 5 && count / n >= 0.9 && name !== "Unclassified") {
      const clade = rank + ":" + name;
      const span = document.createElement("span");
      span.className = "count";
      span.append(document.createTextNode(Math.round(100 * count / n) + "% " + name + ": "));
      const m = document.createElement("a");
      m.href = "/matrix?a=" + encodeURIComponent(clade);
      m.textContent = "matrix of the " + rank;
      const t = document.createElement("a");
      t.href = "/taxa/" + rank + "/" + encodeURIComponent(name);
      t.textContent = "clade page";
      span.append(m, document.createTextNode(" · "), t);
      selLimit.appendChild(span);
    }
  }
  $("[data-land-clear]").addEventListener("click", () => { selection = new Set(); updateSelection(); draw(); });

  function applySearch() {
    const q = search.value.trim().toLowerCase();
    syncUrl();
    if (!D) return;
    if (!q) { highlight = null; matches.textContent = ""; draw(); return; }
    highlight = new Set();
    const hitNames = {};
    RANKS.forEach((r) => { hitNames[r] = new Set(D.ranks[r].names.map((n, k) => (n.toLowerCase() === q ? k : -1)).filter((k) => k >= 0)); });
    const exact = RANKS.some((r) => hitNames[r].size);
    for (let i = 0; i < D.ids.length; i++) {
      const byClade = exact ? RANKS.some((r) => hitNames[r].has(D.ranks[r].idx[i]))
                            : RANKS.some((r) => D.ranks[r].names[D.ranks[r].idx[i]].toLowerCase().includes(q));
      if (byClade || D.ids[i].toLowerCase().includes(q)) highlight.add(i);
    }
    matches.textContent = highlight.size.toLocaleString() + " genome" + (highlight.size === 1 ? "" : "s") + " match" +
      (highlight.size ? " · Enter selects them" : "");
    draw();
  }
  let timer = null;
  search.addEventListener("input", () => { clearTimeout(timer); timer = setTimeout(applySearch, 150); });
  search.addEventListener("keydown", (ev) => {
    if (ev.key !== "Enter") return;
    ev.preventDefault();
    if (highlight && highlight.size) { selection = new Set(highlight); updateSelection(); draw(); }
  });

  let resizeTimer = null;
  window.addEventListener("resize", () => { clearTimeout(resizeTimer); resizeTimer = setTimeout(() => D && layout(), 120); });
  load();
})();
