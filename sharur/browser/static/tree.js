/* Taxonomy tree: radial or rectangular cladogram with feature rings, collapsible clades. */
(function () {
  "use strict";
  const source = document.getElementById("tree-data");
  const svg = document.getElementById("tree-svg");
  if (!source || !svg) return;
  const data = JSON.parse(source.textContent);
  const NS = "http://www.w3.org/2000/svg";
  const host = document.getElementById("tree-host");
  const tip = host.querySelector(".tree-tip");
  const detail = document.getElementById("tree-detail");
  const F = data.features;
  const ROOT = data.tree;
  let layout = data.layout;

  // ---- tree bookkeeping ------------------------------------------------------
  const byId = new Map();
  (function index(node, parent, depth, id) {
    node.id = id; node.parent = parent; node.depth = depth; node.open = false;
    byId.set(id, node);
    node.children.forEach((c, i) => index(c, node, depth + 1, id + "/" + i));
  })(ROOT, null, 0, "r");
  const unclassified = (n) => n.name === "Unclassified";

  function tips() {
    const out = [];
    (function walk(n) {
      if (n.open && n.children.length) n.children.forEach(walk); else out.push(n);
    })(ROOT);
    return out;
  }
  function openTo(rank) {
    const target = data.ranks.indexOf(rank);
    byId.forEach((n) => { n.open = n.children.length > 0 && data.ranks.indexOf(n.rank) < target; });
    ROOT.open = ROOT.children.length > 0;
  }
  function autoOpen(limit) {
    byId.forEach((n) => { n.open = false; });
    ROOT.open = true;
    let frontier = [ROOT];
    for (;;) {
      const next = [];
      frontier.forEach((n) => n.children.forEach((c) => next.push(c)));
      const expandable = next.filter((c) => c.children.length);
      if (!expandable.length) break;
      let count = 0;
      next.forEach((c) => { count += c.children.length || 1; });
      if (count > limit) break;
      expandable.forEach((c) => { c.open = true; });
      frontier = expandable;
    }
  }
  if (data.open) openTo(data.open); else autoOpen(layout === "radial" ? 140 : 90);

  // ---- colors ------------------------------------------------------------------
  const share = (n, i) => (n.n ? n.k[i] / n.n : 0);
  const heat = (color, s) => (s <= 0 ? "var(--surface-3)" :
    "color-mix(in srgb, " + color + " " + Math.round(18 + 82 * s) + "%, var(--surface-3))");
  const weight = (n) => (n.open || !n.children.length ? 1 : Math.min(3.2, 1 + Math.log10(n.n)));

  // ---- svg helpers -----------------------------------------------------------
  function el(name, attrs, text) {
    const e = document.createElementNS(NS, name);
    for (const k in attrs) if (attrs[k] !== undefined && attrs[k] !== null) e.setAttribute(k, attrs[k]);
    if (text !== undefined) e.textContent = text;
    return e;
  }
  const polar = (cx, cy, r, a) => [cx + r * Math.cos(a), cy + r * Math.sin(a)];
  function sector(cx, cy, r0, r1, a0, a1) {
    const large = a1 - a0 > Math.PI ? 1 : 0;
    const [x0, y0] = polar(cx, cy, r1, a0), [x1, y1] = polar(cx, cy, r1, a1);
    const [x2, y2] = polar(cx, cy, r0, a1), [x3, y3] = polar(cx, cy, r0, a0);
    return `M${x0.toFixed(1)},${y0.toFixed(1)} A${r1},${r1} 0 ${large} 1 ${x1.toFixed(1)},${y1.toFixed(1)} ` +
      `L${x2.toFixed(1)},${y2.toFixed(1)} A${r0},${r0} 0 ${large} 0 ${x3.toFixed(1)},${y3.toFixed(1)} Z`;
  }
  function label(n) {
    const collapsed = !n.open && n.children.length;
    return n.name + (collapsed ? " (" + n.n + ")" : "");
  }

  // ---- radial ------------------------------------------------------------------
  function drawRadial() {
    const S = 1000, cx = S / 2, cy = S / 2;
    const rw = F.length <= 3 ? 16 : F.length <= 5 ? 12 : 9, gap = 2;
    const visible = tips();
    const labelSpace = visible.length > 260 ? 70 : 150;
    const rTip = S / 2 - 18 - labelSpace - F.length * (rw + gap) - 8;
    const scale = depthScale(visible);
    const radius = (d) => (d / scale.levels) * rTip;
    const wedgeEnd = (n) => radius(n.depth) + (rTip / scale.levels) * scale.length(n);
    const total = visible.reduce((s, n) => s + weight(n), 0);
    const span = 2 * Math.PI * (visible.length > 1 ? 1 : 0.999);
    const start = -Math.PI / 2;
    let acc = 0;
    visible.forEach((n) => {
      const w = weight(n) / total * span;
      n.a0 = start + acc; n.a1 = start + acc + w; n.a = start + acc + w / 2; acc += w;
    });
    (function place(n) {
      if (n.open && n.children.length) {
        n.children.forEach(place);
        n.a = (n.children[0].a + n.children[n.children.length - 1].a) / 2;
        n.a0 = n.children[0].a0; n.a1 = n.children[n.children.length - 1].a1;
      }
      n.r = radius(n.depth);
    })(ROOT);

    svg.setAttribute("viewBox", `0 0 ${S} ${S}`);
    const gEdges = el("g", { class: "t-edges" }), gNodes = el("g"), gRings = el("g"), gLabels = el("g", { class: "t-labels" });
    (function edges(n) {
      if (!(n.open && n.children.length)) return;
      const first = n.children[0], last = n.children[n.children.length - 1];
      if (n.children.length > 1 && n.r > 0) {
        const [x0, y0] = polar(cx, cy, n.r, first.a), [x1, y1] = polar(cx, cy, n.r, last.a);
        const large = last.a - first.a > Math.PI ? 1 : 0;
        gEdges.append(el("path", { d: `M${x0},${y0} A${n.r},${n.r} 0 ${large} 1 ${x1},${y1}` }));
      }
      n.children.forEach((c) => {
        const [x0, y0] = polar(cx, cy, n.r, c.a), [x1, y1] = polar(cx, cy, c.r, c.a);
        gEdges.append(el("line", { x1: x0, y1: y0, x2: x1, y2: y1 }));
        edges(c);
      });
    })(ROOT);
    visible.forEach((n) => {
      const collapsed = n.children.length && !n.open;
      if (collapsed) {
        const pad = Math.min(0.012, (n.a1 - n.a0) * 0.12);
        const rEnd = wedgeEnd(n);
        const [px, py] = polar(cx, cy, n.r, n.a);
        const [ax, ay] = polar(cx, cy, rEnd, n.a0 + pad), [bx, by] = polar(cx, cy, rEnd, n.a1 - pad);
        const large = n.a1 - n.a0 - 2 * pad > Math.PI ? 1 : 0;
        gNodes.append(el("path", { d: `M${px},${py} L${ax},${ay} A${rEnd},${rEnd} 0 ${large} 1 ${bx},${by} Z`,
          class: "t-wedge", "data-id": n.id }));
        if (rEnd < rTip - 1) {
          const [x0, y0] = polar(cx, cy, rEnd, n.a), [x1, y1] = polar(cx, cy, rTip, n.a);
          gEdges.append(el("line", { x1: x0, y1: y0, x2: x1, y2: y1, class: "t-align" }));
        }
      } else if (n.r < rTip) {
        const [x0, y0] = polar(cx, cy, n.r, n.a), [x1, y1] = polar(cx, cy, rTip, n.a);
        gEdges.append(el("line", { x1: x0, y1: y0, x2: x1, y2: y1, class: "t-align" }));
      }
      F.forEach((f, i) => {
        const r0 = rTip + 6 + i * (rw + gap);
        const pad = (n.a1 - n.a0) > 0.01 ? 0.0015 : 0;
        gRings.append(el("path", { d: sector(cx, cy, r0, r0 + rw, n.a0 + pad, n.a1 - pad), class: "t-cell",
          style: "fill:" + heat(f.color, share(n, i)), "data-id": n.id, "data-f": i }));
      });
      const rLab = rTip + 10 + F.length * (rw + gap);
      const room = (n.a1 - n.a0) * rLab;
      if (room >= 8.5) {
        const deg = n.a * 180 / Math.PI;
        const flip = deg > 90 && deg < 270;
        const [lx, ly] = polar(cx, cy, rLab, n.a);
        const t = el("text", { x: lx, y: ly, class: "t-label" + (n.children.length && !n.open ? " t-collapsed" : ""),
          "text-anchor": flip ? "end" : "start", "dominant-baseline": "middle", "data-id": n.id,
          transform: `rotate(${flip ? deg + 180 : deg} ${lx} ${ly})` }, clip(label(n), labelSpace / 6.2));
        gLabels.append(t);
      }
    });
    byId.forEach((n) => {
      if (!(n.open && n.children.length) || n === ROOT) return;
      const [x, y] = polar(cx, cy, n.r, n.a);
      gNodes.append(el("circle", { cx: x, cy: y, r: 2.5 + 3 * Math.sqrt(n.n / ROOT.n), class: "t-node", "data-id": n.id,
        style: F.length ? "fill:" + heat(F[0].color, share(n, 0)) : null }));
    });
    gNodes.append(el("circle", { cx, cy, r: 5, class: "t-node t-root", "data-id": ROOT.id }));
    svg.replaceChildren(gEdges, gRings, gNodes, gLabels);
  }

  // ---- rectangular ---------------------------------------------------------------
  function drawRect() {
    const W = 1000, cw = F.length <= 4 ? 24 : 18, gap = 3;
    const visible = tips();
    const head = F.length ? 120 : 16;
    const longest = Math.max(...visible.map((n) => Math.min(34, label(n).length)));
    const labelW = Math.max(90, Math.min(230, longest * 6.4 + 16));
    const heatW = F.length * (cw + gap);
    const treeW = Math.max(280, W - 24 - labelW - heatW - 10);
    const scale = depthScale(visible);
    const x = (d) => 12 + (d / scale.levels) * treeW;
    const xTip = x(scale.levels);
    const wedgeEnd = (n) => x(n.depth) + (treeW / scale.levels) * scale.length(n);
    let y = head;
    visible.forEach((n) => {
      const h = 17 * weight(n);
      n.y0 = y; n.y1 = y + h; n.y = y + h / 2; y += h;
    });
    const H = y + 14;
    (function place(n) {
      if (n.open && n.children.length) {
        n.children.forEach(place);
        n.y = (n.children[0].y + n.children[n.children.length - 1].y) / 2;
      }
      n.x = x(n.depth);
    })(ROOT);
    svg.setAttribute("viewBox", `0 0 ${W} ${H}`);
    const gEdges = el("g", { class: "t-edges" }), gNodes = el("g"), gRings = el("g"), gLabels = el("g", { class: "t-labels" });
    (function edges(n) {
      if (!(n.open && n.children.length)) return;
      const first = n.children[0], last = n.children[n.children.length - 1];
      gEdges.append(el("line", { x1: n.x, y1: first.y, x2: n.x, y2: last.y }));
      n.children.forEach((c) => { gEdges.append(el("line", { x1: n.x, y1: c.y, x2: c.x, y2: c.y })); edges(c); });
    })(ROOT);
    const heatX = xTip + labelW;
    visible.forEach((n) => {
      if (n.children.length && !n.open) {
        const xEnd = wedgeEnd(n);
        gNodes.append(el("path", { d: `M${n.x},${n.y} L${xEnd},${n.y0 + 2} L${xEnd},${n.y1 - 2} Z`, class: "t-wedge", "data-id": n.id }));
        if (xEnd < xTip - 1) gEdges.append(el("line", { x1: xEnd, y1: n.y, x2: xTip, y2: n.y, class: "t-align" }));
      } else if (n.x < xTip) {
        gEdges.append(el("line", { x1: n.x, y1: n.y, x2: xTip, y2: n.y, class: "t-align" }));
      }
      gLabels.append(el("text", { x: xTip + 6, y: n.y, "dominant-baseline": "middle", "data-id": n.id,
        class: "t-label" + (n.children.length && !n.open ? " t-collapsed" : "") }, clip(label(n), labelW / 6.4)));
      F.forEach((f, i) => {
        gRings.append(el("rect", { x: heatX + i * (cw + gap), y: n.y0 + 1, width: cw, height: Math.max(2, n.y1 - n.y0 - 2), rx: 2,
          class: "t-cell", style: "fill:" + heat(f.color, share(n, i)), "data-id": n.id, "data-f": i }));
      });
    });
    F.forEach((f, i) => {
      const hx = heatX + i * (cw + gap) + cw / 2 + 4;
      gLabels.append(el("text", { x: hx, y: head - 8, class: "t-head", transform: `rotate(-55 ${hx} ${head - 8})` }, clip(f.label, 22)));
      gLabels.append(el("rect", { x: heatX + i * (cw + gap), y: head - 5, width: cw, height: 3, rx: 1.5, style: "fill:" + f.color }));
    });
    byId.forEach((n) => {
      if (!(n.open && n.children.length)) return;
      gNodes.append(el("circle", { cx: n.x, cy: n.y, r: n === ROOT ? 4.5 : 2.5 + 3 * Math.sqrt(n.n / ROOT.n),
        class: "t-node" + (n === ROOT ? " t-root" : ""), "data-id": n.id,
        style: F.length && n !== ROOT ? "fill:" + heat(F[0].color, share(n, 0)) : null }));
    });
    svg.replaceChildren(gEdges, gRings, gNodes, gLabels);
  }

  // Depth axis for the clades on screen: their deepest level, plus one level for closed clades,
  // whose wedge length grows with genome count (log scale).
  function depthScale(visible) {
    const collapsed = visible.filter((n) => n.children.length && !n.open);
    const deepest = Math.max(1, ...visible.map((n) => n.depth));
    const levels = Math.max(1, deepest + (collapsed.length ? 1 : 0));
    const maxN = Math.max(2, ...collapsed.map((n) => n.n));
    return { levels, length: (n) => 0.3 + 0.7 * Math.log(n.n + 1) / Math.log(maxN + 1) };
  }

  function clip(text, chars) {
    const n = Math.max(6, Math.floor(chars));
    return text.length > n ? text.slice(0, n - 1) + "…" : text;
  }
  function draw() {
    svg.classList.toggle("t-rect", layout === "rect");
    if (layout === "rect") drawRect(); else drawRadial();
  }

  // ---- interaction ---------------------------------------------------------------------
  let selected = null;
  function pct(x) { return (x * 100).toFixed(x > 0 && x < 0.1 ? 1 : 0) + "%"; }
  function tooltip(ev) {
    const target = ev.target.closest("[data-id]");
    if (!target) { tip.hidden = true; return; }
    const n = byId.get(target.dataset.id);
    tip.replaceChildren();
    const head = document.createElement("b"); head.textContent = n.name;
    const sub = document.createElement("div"); sub.className = "count";
    sub.textContent = (n.rank === "root" ? "" : n.rank + " · ") + n.n.toLocaleString() + " genomes" +
      (n.children.length ? " · " + n.children.length + " " + (n.children[0].rank || "") + (n.children.length > 1 ? " clades" : " clade") : "");
    tip.append(head, sub);
    F.forEach((f, i) => {
      const row = document.createElement("div"); row.className = "tree-tip-row" + (target.dataset.f === String(i) ? " on" : "");
      const sw = document.createElement("i"); sw.style.background = f.color;
      row.append(sw, document.createTextNode(f.label + ": " + n.k[i] + " / " + n.n + " (" + pct(share(n, i)) + ")"));
      tip.append(row);
    });
    tip.hidden = false;
    const r = host.getBoundingClientRect();
    let left = ev.clientX - r.left + host.scrollLeft + 14, top = ev.clientY - r.top + host.scrollTop + 14;
    if (left + 300 > host.scrollLeft + host.clientWidth) left -= 320;
    tip.style.left = left + "px"; tip.style.top = top + "px";
  }
  function url(path, params) {
    const q = new URLSearchParams(params);
    return path + "?" + q.toString();
  }
  function featureTokens() { return F.map((f) => f.token).join(","); }
  function showDetail(n) {
    selected = n;
    detail.replaceChildren();
    const h = document.createElement("h2"); h.textContent = n.name;
    const sub = document.createElement("p"); sub.className = "count";
    sub.textContent = (n.rank === "root" ? "dataset" : n.rank) + " · " + n.n.toLocaleString() + " genomes";
    detail.append(h, sub);
    if (F.length) {
      const list = document.createElement("div"); list.className = "tree-shares";
      F.forEach((f, i) => {
        const row = document.createElement("div"); row.className = "tree-share";
        const name = document.createElement("span"); name.className = "tree-share-name"; name.textContent = f.label; name.title = f.name || f.label;
        const bar = document.createElement("div"); bar.className = "bar"; bar.style.setProperty("--bar", f.color);
        const fill = document.createElement("i"); fill.style.width = (share(n, i) * 100).toFixed(1) + "%"; bar.append(fill);
        const val = document.createElement("span"); val.className = "count"; val.textContent = n.k[i] + "/" + n.n;
        row.append(name, bar, val); list.append(row);
      });
      detail.append(list);
    }
    const actions = document.createElement("p"); actions.className = "actions";
    const real = n.rank !== "root" && !unclassified(n);
    if (n.children.length) {
      const b = document.createElement("button"); b.type = "button"; b.textContent = n.open ? "Close" : "Open";
      b.addEventListener("click", () => { toggle(n); showDetail(n); }); actions.append(b);
    }
    if (real) {
      const zoom = document.createElement("a"); zoom.className = "btn";
      zoom.href = url("/tree", { root: n.rank + ":" + n.name, layout, features: featureTokens() });
      zoom.textContent = "Zoom in"; actions.append(zoom);
      const page = document.createElement("a"); page.className = "btn";
      page.href = "/taxa/" + encodeURIComponent(n.rank) + "/" + encodeURIComponent(n.name); page.textContent = "Clade page";
      actions.append(page);
      const kinds = new Set(F.map((f) => f.kind));
      if (F.length && kinds.size === 1) {
        const kind = F[0].kind;
        const m = document.createElement("a"); m.className = "btn";
        m.href = url("/matrix", { a: n.rank + ":" + n.name, kind, features: F.map((f) => f.id).join(",") });
        m.textContent = "Genome matrix"; actions.append(m);
      }
    }
    detail.append(actions);
    svg.querySelectorAll(".t-selected").forEach((e) => e.classList.remove("t-selected"));
    svg.querySelectorAll(`[data-id="${n.id}"]`).forEach((e) => e.classList.add("t-selected"));
  }
  function toggle(n) {
    if (!n.children.length) return;
    n.open = !n.open;
    if (!n.open) byId.forEach((m) => { if (m.id.startsWith(n.id + "/")) m.open = false; });
    draw();
  }
  svg.addEventListener("mousemove", tooltip);
  svg.addEventListener("mouseleave", () => { tip.hidden = true; });
  svg.addEventListener("click", (ev) => {
    const target = ev.target.closest("[data-id]");
    if (!target) return;
    const n = byId.get(target.dataset.id);
    const shape = target.classList.contains("t-wedge") || target.classList.contains("t-node");
    if (shape && n !== ROOT) toggle(n);
    showDetail(n);
  });

  document.querySelectorAll("[data-layout]").forEach((b) => b.addEventListener("click", () => {
    layout = b.dataset.layout;
    document.querySelectorAll("[data-layout]").forEach((x) => x.classList.toggle("on", x === b));
    const u = new URL(location.href); u.searchParams.set("layout", layout); history.replaceState(null, "", u);
    draw();
    if (selected) showDetail(selected);
  }));
  const openSel = document.getElementById("tree-open");
  if (openSel) openSel.addEventListener("change", () => {
    if (openSel.value) openTo(openSel.value); else autoOpen(layout === "radial" ? 140 : 90);
    const u = new URL(location.href);
    if (openSel.value) u.searchParams.set("open", openSel.value); else u.searchParams.delete("open");
    history.replaceState(null, "", u);
    draw();
  });

  // feature picker
  const form = document.querySelector(".tree-add");
  if (form) {
    const input = form.querySelector(".picker-input"), list = form.querySelector(".picker-list");
    let timer = null;
    form.addEventListener("submit", (e) => e.preventDefault());
    input.addEventListener("input", () => {
      clearTimeout(timer);
      const q = input.value.trim();
      if (q.length < 2) { list.classList.remove("open"); return; }
      timer = setTimeout(() => {
        fetch("/api/tree/features?q=" + encodeURIComponent(q)).then((r) => r.json()).then((rows) => {
          list.replaceChildren();
          rows.forEach((r) => {
            const a = document.createElement("a"); a.href = "#";
            const name = document.createElement("span"); name.textContent = r.label + (r.id !== r.label ? "  " + r.id : "");
            const kind = document.createElement("span"); kind.className = "kind"; kind.textContent = r.kind;
            a.append(name, kind);
            a.addEventListener("click", (e) => {
              e.preventDefault();
              const current = form.dataset.current ? form.dataset.current.split(",") : [];
              if (!current.includes(r.token)) current.push(r.token);
              location.href = url("/tree", { root: form.dataset.root, layout, features: current.join(",") });
            });
            list.append(a);
          });
          list.classList.toggle("open", rows.length > 0);
        });
      }, 160);
    });
    input.addEventListener("blur", () => setTimeout(() => list.classList.remove("open"), 200));
  }

  draw();
})();
