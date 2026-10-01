/* Theme switch and the ⌘K command palette. */
(function () {
  "use strict";
  const read = (k, d) => { try { return JSON.parse(localStorage.getItem(k)) ?? d; } catch (e) { return d; } };

  // ---- theme: system → light → dark ---------------------------------------
  const root = document.documentElement;
  const themeBtn = document.getElementById("theme-toggle");
  const names = { auto: "follows the system", light: "light", dark: "dark" };
  function applyTheme(t) {
    if (t === "light" || t === "dark") root.dataset.theme = t; else delete root.dataset.theme;
    try { if (t === "auto") localStorage.removeItem("sharur-theme"); else localStorage.setItem("sharur-theme", t); } catch (e) { /* private mode */ }
    if (themeBtn) { themeBtn.dataset.mode = t; themeBtn.title = "Theme: " + names[t]; }
  }
  const current = () => root.dataset.theme || "auto";
  const nextTheme = () => ({ auto: "light", light: "dark", dark: "auto" })[current()];
  if (themeBtn) {
    themeBtn.dataset.mode = current();
    themeBtn.title = "Theme: " + names[current()];
    themeBtn.addEventListener("click", () => applyTheme(nextTheme()));
  }

  // ---- command palette ------------------------------------------------------
  const pal = document.getElementById("palette");
  if (!pal) return;
  const input = pal.querySelector(".palette-input");
  const list = pal.querySelector(".palette-list");
  const PAGES = [
    ["Overview", "/", "home"], ["Tree of life", "/taxa", "taxa clades treemap"],
    ["Taxonomy tree", "/tree", "cladogram map features rings"], ["Functional landscape", "/landscape", "map genomes pca"],
    ["Functions", "/functions", "labels predicates"], ["Pathways", "/pathways", "kegg modules"],
    ["Families", "/domains", "pfam vog domains"], ["Systems", "/systems", "defense secretion crispr"],
    ["CRISPR arrays", "/crispr", "minced repeats spacers"], ["Discover", "/discover", "giants fusions unannotated"],
    ["Presence/absence matrix", "/matrix", "genomes features ko"], ["Compare", "/compare", "two genomes clades"],
    ["Function heatmap", "/heatmap", "prevalence clades"], ["Architecture search", "/architecture", "domain pattern"],
    ["Data health", "/health", "checks quality"], ["Flags & notes", "/flags", "curation"],
    ["Triage", "/triage", "curation queue"], ["Collection", "/collection", "saved"], ["Findings", "/findings", "agents"],
    ["JSON API", "/api", "json export scripts programmatic"],
  ];
  const ACTIONS = [["Switch theme", () => applyTheme(nextTheme()), "dark light mode"],
                   ["Surprise me", () => { window.location = "/discover/random"; }, "random"]];
  let rows = [], sel = 0, timer = null, seq = 0;

  function fuzzy(text, q) {
    text = text.toLowerCase(); q = q.toLowerCase();
    if (text.includes(q)) return 2 + (text.startsWith(q) ? 1 : 0);
    let i = 0;
    for (const c of text) if (c === q[i]) i++;
    return i === q.length ? 1 : 0;
  }
  function draw() {
    list.innerHTML = "";
    rows.forEach((r, i) => {
      const a = document.createElement("a");
      a.href = r.url || "#";
      a.className = i === sel ? "on" : "";
      const label = document.createElement("span");
      label.textContent = r.label;
      if (r.sub) { const s = document.createElement("span"); s.className = "sub"; s.textContent = " · " + r.sub; label.appendChild(s); }
      const kind = document.createElement("span");
      kind.className = "kind";
      kind.textContent = r.kind;
      a.append(label, kind);
      a.addEventListener("mouseenter", () => { sel = i; mark(); });
      a.addEventListener("click", (e) => { e.preventDefault(); go(r); });
      list.appendChild(a);
    });
  }
  function mark() { Array.from(list.children).forEach((a, i) => a.classList.toggle("on", i === sel)); }
  function go(r) {
    if (!r) return;
    if (r.run) { close(); r.run(); return; }
    window.location = r.url;
  }
  function local(q) {
    if (!q) {
      const recent = read("sharur_recent", []).slice(0, 6).map((r) => ({ label: r.title, url: r.url, kind: "recent" }));
      return recent.concat(PAGES.slice(0, 8).map(([label, url]) => ({ label, url, kind: "page" })));
    }
    const pages = PAGES.map(([label, url, words]) => ({ label, url, kind: "page", score: Math.max(fuzzy(label, q), fuzzy(words, q) - 1) }));
    const actions = ACTIONS.map(([label, run, words]) => ({ label, run, kind: "action", score: Math.max(fuzzy(label, q), fuzzy(words, q) - 1) }));
    return pages.concat(actions).filter((r) => r.score > 0).sort((a, b) => b.score - a.score).slice(0, 6);
  }
  function update() {
    const q = input.value.trim();
    rows = local(q); sel = 0; draw();
    clearTimeout(timer);
    if (q.length < 2) return;
    const mine = ++seq;
    timer = setTimeout(() => {
      fetch("/api/suggest?q=" + encodeURIComponent(q)).then((r) => r.json()).then((found) => {
        if (mine !== seq) return;
        rows = local(q).concat(found.map((f) => ({ label: f.label, sub: f.sub, url: f.url, kind: f.kind })));
        sel = Math.min(sel, Math.max(rows.length - 1, 0)); draw();
      }).catch(() => {});
    }, 120);
  }
  function open() { pal.hidden = false; input.value = ""; update(); input.focus(); }
  function close() { pal.hidden = true; }
  document.getElementById("palette-open")?.addEventListener("click", open);
  document.addEventListener("keydown", (e) => {
    if ((e.metaKey || e.ctrlKey) && e.key.toLowerCase() === "k") { e.preventDefault(); pal.hidden ? open() : close(); }
  });
  pal.addEventListener("click", (e) => { if (e.target === pal) close(); });
  input.addEventListener("input", update);
  input.addEventListener("keydown", (e) => {
    if (e.key === "Escape") { e.preventDefault(); close(); }
    else if (e.key === "ArrowDown") { e.preventDefault(); sel = Math.min(sel + 1, rows.length - 1); mark(); }
    else if (e.key === "ArrowUp") { e.preventDefault(); sel = Math.max(sel - 1, 0); mark(); }
    else if (e.key === "Enter") { e.preventDefault(); go(rows[sel]); }
  });
})();
