/* Presence/absence matrix: canvas renderer with hover details and click-through. */
(function () {
  "use strict";

  // group pickers: taxon or genome suggestions fill the field with rank:name or the genome id
  document.querySelectorAll(".matrix-form .group-input").forEach((input) => {
    const list = input.parentElement.querySelector(".picker-list");
    let timer = null;
    input.addEventListener("input", () => {
      clearTimeout(timer);
      const q = input.value.trim();
      if (q.length < 2 || q.includes(":")) { list.classList.remove("open"); return; }
      timer = setTimeout(() => {
        fetch("/api/suggest?kind=taxon,genome&q=" + encodeURIComponent(q)).then((r) => r.json()).then((rows) => {
          list.innerHTML = "";
          rows.forEach((r) => {
            const a = document.createElement("a");
            a.href = "#";
            const name = document.createElement("span");
            name.textContent = r.label;
            const kind = document.createElement("span");
            kind.className = "kind";
            kind.textContent = r.kind === "taxon" ? r.sub : "genome";
            a.append(name, kind);
            a.addEventListener("click", (e) => {
              e.preventDefault();
              input.value = r.kind === "taxon" ? r.sub + ":" + r.label : r.label;
              list.classList.remove("open");
            });
            list.appendChild(a);
          });
          list.classList.toggle("open", rows.length > 0);
        });
      }, 150);
    });
    input.addEventListener("blur", () => setTimeout(() => list.classList.remove("open"), 200));
  });

  const host = document.getElementById("matrix");
  const source = document.getElementById("matrix-data");
  if (!host || !source) return;
  const data = JSON.parse(source.textContent);
  const canvas = host.querySelector("canvas");
  const tip = host.querySelector(".matrix-tip");
  const ctx = canvas.getContext("2d");
  const css = getComputedStyle(document.documentElement);
  const color = (name, fallback) => (css.getPropertyValue(name) || fallback).trim();
  const C = {
    ink: color("--ink", "#161b22"), ink2: color("--ink-2", "#47505c"), ink3: color("--ink-3", "#7c8591"),
    line: color("--line", "#e3e0d8"), bg: color("--surface", "#fff"),
    empty: color("--mx-empty", "#f1efe9"), c1: color("--mx-1", "#8fd3c7"), c2: color("--mx-2", "#2fa894"),
    c3: color("--mx-3", "#0b5c52"), comp: color("--accent", "#0f7c6e"), hatch: color("--mx-hatch", "#d9d4c8"),
    warn: color("--warn", "#c2410c"),
  };

  const G = data.genomes, F = data.features, V = data.values;
  const nG = G.length, nF = F.length;
  const module = data.kind === "module";   // cells hold module completeness
  const stepped = F.some((f) => f.step !== null && f.step !== undefined);   // KOs laid out by module step

  // layout
  const rowH = Math.max(2, Math.min(14, Math.floor(1600 / Math.max(nG, 1))));
  const showRowLabels = rowH >= 9;
  const labelW = showRowLabels ? 170 : 0;
  const stripW = 8, compW = 26, gap = 6, groupGap = 10;
  const leftW = labelW + stripW + 4 + compW + gap;
  const headerH = 150;
  const width = Math.max(host.clientWidth, 320);
  const cellW = Math.max(6, Math.min(24, Math.floor((width - leftW - 20) / Math.max(nF, 1))));
  const gridW = cellW * nF;
  const groupBreaks = [];
  const rowY = new Array(nG);
  let y = headerH;
  for (let i = 0; i < nG; i++) {
    if (i > 0 && G[i].group !== G[i - 1].group) { y += groupGap; groupBreaks.push(y - groupGap / 2); }
    rowY[i] = y;
    y += rowH;
  }
  const height = y + 4;
  const totalW = leftW + gridW + 8;

  const dpr = window.devicePixelRatio || 1;
  canvas.width = Math.ceil(totalW * dpr);
  canvas.height = Math.ceil(height * dpr);
  canvas.style.width = totalW + "px";
  canvas.style.height = height + "px";
  ctx.scale(dpr, dpr);

  // hatch for absences in incomplete genomes
  const pat = document.createElement("canvas");
  pat.width = pat.height = 6;
  const pctx = pat.getContext("2d");
  pctx.fillStyle = C.empty; pctx.fillRect(0, 0, 6, 6);
  pctx.strokeStyle = C.hatch; pctx.lineWidth = 1.2;
  pctx.beginPath(); pctx.moveTo(0, 6); pctx.lineTo(6, 0); pctx.stroke();
  const hatch = ctx.createPattern(pat, "repeat");

  function cellFill(v, completeness) {
    if (module) {
      if (v <= 0) return completeness !== null && completeness < 70 ? hatch : C.empty;
      ctx.globalAlpha = 0.25 + 0.75 * v;
      return v >= data.threshold ? C.c3 : C.c2;
    }
    if (v < data.threshold) return completeness !== null && completeness < 70 ? hatch : C.empty;
    return v >= 3 ? C.c3 : v >= 2 ? C.c2 : C.c1;
  }

  function draw() {
    ctx.clearRect(0, 0, totalW, height);
    ctx.font = "11px ui-sans-serif, system-ui, sans-serif";
    // module steps: alternate header bands
    if (stepped) {
      for (let j = 0; j < nF; j++) {
        if ((F[j].step || 0) % 2 === 1) {
          ctx.fillStyle = C.line;
          ctx.fillRect(leftW + j * cellW, headerH - 10, cellW, 8);
        }
      }
    }
    // column labels, rotated
    ctx.fillStyle = C.ink2;
    for (let j = 0; j < nF; j++) {
      ctx.save();
      ctx.translate(leftW + j * cellW + cellW / 2 + 3, headerH - 14);
      ctx.rotate(-Math.PI / 3);
      const label = F[j].label.length > 26 ? F[j].label.slice(0, 25) + "…" : F[j].label;
      ctx.fillText(label, 0, 0);
      ctx.restore();
    }
    // rows
    for (let i = 0; i < nG; i++) {
      const g = G[i], yy = rowY[i];
      if (showRowLabels) {
        ctx.fillStyle = C.ink2;
        const text = g.id.length > 26 ? g.id.slice(0, 25) + "…" : g.id;
        ctx.fillText(text, 2, yy + rowH - 2);
      }
      ctx.fillStyle = g.color || C.ink3;
      ctx.fillRect(labelW, yy, stripW, rowH - (rowH > 4 ? 1 : 0));
      if (g.completeness !== null) {
        ctx.fillStyle = C.line;
        ctx.fillRect(labelW + stripW + 4, yy, compW, rowH - (rowH > 4 ? 1 : 0));
        ctx.fillStyle = g.deficit ? C.warn : g.completeness < 70 ? C.hatch : C.comp;
        ctx.fillRect(labelW + stripW + 4, yy, compW * Math.min(g.completeness, 100) / 100, rowH - (rowH > 4 ? 1 : 0));
      }
      const row = V[i];
      for (let j = 0; j < nF; j++) {
        ctx.fillStyle = cellFill(row[j], g.completeness);
        ctx.fillRect(leftW + j * cellW, yy, cellW - (cellW > 7 ? 1 : 0), rowH - (rowH > 4 ? 1 : 0));
        ctx.globalAlpha = 1;
      }
    }
    // group separators and step separators
    ctx.strokeStyle = C.ink3; ctx.lineWidth = 1;
    groupBreaks.forEach((yy) => { ctx.beginPath(); ctx.moveTo(0, yy); ctx.lineTo(totalW, yy); ctx.stroke(); });
    if (stepped) {
      ctx.strokeStyle = C.ink3;
      for (let j = 1; j < nF; j++) {
        if (F[j].step !== F[j - 1].step) {
          const x = leftW + j * cellW - 0.5;
          ctx.beginPath(); ctx.moveTo(x, headerH - 10); ctx.lineTo(x, height); ctx.stroke();
        }
      }
    }
  }
  draw();

  // hit testing
  function locate(ev) {
    const r = canvas.getBoundingClientRect();
    const x = ev.clientX - r.left, yy = ev.clientY - r.top;
    const j = Math.floor((x - leftW) / cellW);
    let i = -1;
    if (yy >= headerH) {
      // binary search over row starts
      let lo = 0, hi = nG - 1;
      while (lo <= hi) {
        const mid = (lo + hi) >> 1;
        if (rowY[mid] <= yy) { i = mid; lo = mid + 1; } else hi = mid - 1;
      }
      if (i >= 0 && yy >= rowY[i] + rowH) i = -1;
    }
    return { x, y: yy, i, j: x >= leftW && j < nF ? j : -1, header: yy < headerH };
  }

  function fmt(v) {
    if (module) return v > 0 ? Math.round(v * 100) + "% complete" : "no steps found";
    if (v <= 0) return "not detected";
    return v + (v === 1 ? " copy" : " copies");
  }

  function show(ev) {
    const at = locate(ev);
    let html = "";
    if (at.i >= 0) {
      const g = G[at.i];
      const comp = (g.completeness === null ? "completeness unknown" :
        g.completeness.toFixed(1) + "% complete" + (g.contamination !== null ? ", " + g.contamination.toFixed(1) + "% contamination" : "")) +
        " · " + g.proteins + " proteins" + (g.deficit ? " (gene calls likely missing)" : "");
      html = "<b></b><div class='count lineage'></div><div class='count comp'></div>";
      tip.innerHTML = html;
      tip.querySelector("b").textContent = g.id;
      tip.querySelector(".lineage").textContent = g.lineage;
      tip.querySelector(".comp").textContent = comp;
      if (at.j >= 0) {
        const f = F[at.j], div = document.createElement("div");
        div.className = "feat";
        div.textContent = f.id + " · " + f.name + ": " + fmt(V[at.i][at.j]);
        tip.appendChild(div);
      }
    } else if (at.j >= 0 && at.header) {
      const f = F[at.j];
      tip.innerHTML = "<b></b><div class='count'></div>";
      tip.querySelector("b").textContent = f.id;
      tip.querySelector(".count").textContent = f.name + (f.step !== null && f.step !== undefined ? " · step " + (f.step + 1) : "");
    } else { tip.hidden = true; canvas.style.cursor = "default"; return; }
    tip.hidden = false;
    canvas.style.cursor = "pointer";
    const hostR = host.getBoundingClientRect();
    let left = ev.clientX - hostR.left + host.scrollLeft + 14, top = ev.clientY - hostR.top + host.scrollTop + 14;
    if (left + 280 > host.scrollLeft + host.clientWidth) left -= 300;
    tip.style.left = left + "px";
    tip.style.top = top + "px";
  }

  canvas.addEventListener("mousemove", show);
  canvas.addEventListener("mouseleave", () => { tip.hidden = true; });
  canvas.addEventListener("click", (ev) => {
    const at = locate(ev);
    let href = null;
    if (at.i >= 0 && at.j >= 0) {
      href = data.cell_url.replace("{g}", encodeURIComponent(G[at.i].id)).replace("{f}", encodeURIComponent(F[at.j].id));
    } else if (at.i >= 0 && at.x < leftW) {
      href = "/genome/" + encodeURIComponent(G[at.i].id);
    } else if (at.j >= 0 && at.header) {
      href = F[at.j].url;
    }
    if (href) {
      if (ev.metaKey || ev.ctrlKey) window.open(href, "_blank"); else window.location = href;
    }
  });
})();
