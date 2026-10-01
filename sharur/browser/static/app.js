// Sharur browser: search suggestions, sortable and filterable tables.
(function () {
  "use strict";

  // ---- search suggestions -------------------------------------------------
  const input = document.querySelector(".search input");
  const box = document.querySelector(".suggest");
  let timer = null, items = [], selected = -1;

  function render(results) {
    items = results;
    selected = -1;
    box.innerHTML = "";
    results.forEach((r) => {
      const a = document.createElement("a");
      a.href = r.url;
      const left = document.createElement("span");
      left.textContent = r.label;
      if (r.sub) {
        const sub = document.createElement("span");
        sub.className = "sub";
        sub.textContent = " · " + r.sub;
        left.appendChild(sub);
      }
      const kind = document.createElement("span");
      kind.className = "kind";
      kind.textContent = r.kind;
      a.append(left, kind);
      box.appendChild(a);
    });
    box.classList.toggle("open", results.length > 0);
  }

  if (input && box) {
    input.addEventListener("input", () => {
      clearTimeout(timer);
      const q = input.value.trim();
      if (q.length < 2) { render([]); return; }
      timer = setTimeout(() => {
        fetch("/api/suggest?q=" + encodeURIComponent(q))
          .then((r) => r.json()).then(render).catch(() => render([]));
      }, 120);
    });
    input.addEventListener("keydown", (e) => {
      const links = box.querySelectorAll("a");
      if (e.key === "ArrowDown" || e.key === "ArrowUp") {
        if (!links.length) return;
        e.preventDefault();
        selected = (selected + (e.key === "ArrowDown" ? 1 : -1) + links.length) % links.length;
        links.forEach((l, i) => l.classList.toggle("sel", i === selected));
      } else if (e.key === "Enter" && selected >= 0 && links[selected]) {
        e.preventDefault();
        window.location = links[selected].href;
      } else if (e.key === "Escape") {
        render([]);
      }
    });
    document.addEventListener("click", (e) => { if (!box.contains(e.target) && e.target !== input) render([]); });
    document.addEventListener("keydown", (e) => {
      if (e.key === "/" && document.activeElement.tagName !== "INPUT") { e.preventDefault(); input.focus(); }
    });
  }

  // ---- copy buttons --------------------------------------------------------
  document.querySelectorAll("[data-copy]").forEach((button) => {
    button.addEventListener("click", async () => {
      const source = document.getElementById(button.dataset.copy);
      if (!source) return;
      const text = source.value;
      let ok = false;
      if (navigator.clipboard && window.isSecureContext) {
        // a pending permission prompt can leave writeText unsettled; never let the button hang
        const timeout = new Promise((_, reject) => setTimeout(() => reject(new Error("timeout")), 700));
        try { await Promise.race([navigator.clipboard.writeText(text), timeout]); ok = true; } catch (e) { ok = false; }
      }
      if (!ok) {  // plain-http sharing: fall back to a selection copy
        source.classList.remove("visually-hidden");
        source.select();
        try { ok = document.execCommand("copy"); } catch (e) { ok = false; }
        source.classList.add("visually-hidden");
        window.getSelection().removeAllRanges();
      }
      button.textContent = ok ? "Copied ✓" : "Copy failed";
      button.classList.toggle("done", ok);
      setTimeout(() => { button.textContent = button.dataset.label; button.classList.remove("done"); }, 1600);
    });
  });

  // ---- tables -------------------------------------------------------------
  document.querySelectorAll("table[data-table]").forEach((table) => {
    const body = table.tBodies[0];
    if (!body) return;
    table.querySelectorAll("th[data-sort]").forEach((th, _, all) => {
      th.addEventListener("click", () => {
        const index = Array.from(th.parentNode.children).indexOf(th);
        const numeric = th.dataset.sort === "num";
        const dir = th.classList.contains("desc") ? 1 : -1;
        all.forEach((h) => h.classList.remove("asc", "desc"));
        th.classList.add(dir === 1 ? "asc" : "desc");
        const rows = Array.from(body.rows);
        rows.sort((a, b) => {
          const x = a.cells[index].dataset.v ?? a.cells[index].textContent.trim();
          const y = b.cells[index].dataset.v ?? b.cells[index].textContent.trim();
          if (numeric) return dir * ((parseFloat(x) || 0) - (parseFloat(y) || 0));
          return dir * x.localeCompare(y);
        });
        rows.forEach((r) => body.appendChild(r));
        table.dispatchEvent(new Event("rows-sorted"));
      });
    });
    const filter = document.querySelector(`input[data-filter="${table.id}"]`);
    const counter = document.querySelector(`[data-count="${table.id}"]`);
    const CAP = 40;
    let expanded = false;
    // long tables show their first rows; a toggle reveals the rest, and filtering always searches everything
    const toggle = document.createElement("button");
    toggle.type = "button";
    toggle.className = "show-all";
    toggle.addEventListener("click", () => { expanded = !expanded; update(); });
    (table.closest(".scroll") || table).insertAdjacentElement("afterend", toggle);
    function update() {
      const q = filter ? filter.value.trim().toLowerCase() : "";
      let matches = 0, shown = 0;
      Array.from(body.rows).forEach((r) => {
        const hit = !q || r.textContent.toLowerCase().includes(q);
        if (hit) matches++;
        const visible = hit && (expanded || matches <= CAP);
        r.style.display = visible ? "" : "none";
        if (visible) shown++;
      });
      toggle.hidden = matches <= CAP;
      toggle.textContent = expanded ? "Show the first " + CAP : "Show all " + matches.toLocaleString() + " rows";
      if (counter) counter.textContent = (shown < matches ? shown.toLocaleString() + " shown · " : "") +
        matches.toLocaleString() + " of " + body.rows.length.toLocaleString();
    }
    table.addEventListener("rows-sorted", update);
    if (filter) filter.addEventListener("input", update);
    update();
  });
})();

// ======================================================================
// Curation, collection, navigation, previews, export
// ======================================================================
(function () {
  "use strict";
  const FLAGS = [["verified", "Verified", "✓"], ["suspicious", "Suspicious", "!"],
                 ["interesting", "Interesting", "★"], ["follow_up", "Follow up", "↻"]];
  const store = {
    get(key, fallback) { try { return JSON.parse(localStorage.getItem(key)) ?? fallback; } catch (e) { return fallback; } },
    set(key, value) { try { localStorage.setItem(key, JSON.stringify(value)); } catch (e) { /* private mode */ } },
  };
  const cookie = (name) => (document.cookie.split("; ").find((c) => c.startsWith(name + "=")) || "").split("=").slice(1).join("=");
  const el = (tag, attrs = {}, text) => { const n = document.createElement(tag); Object.assign(n, attrs); if (text !== undefined) n.textContent = text; return n; };
  const typing = () => ["INPUT", "TEXTAREA", "SELECT"].includes(document.activeElement.tagName) || document.activeElement.isContentEditable;

  // ---- who am I -----------------------------------------------------------
  const nameEl = document.getElementById("whoami-name");
  const myName = () => decodeURIComponent(cookie("sharur_user")) || "anonymous";
  if (nameEl) nameEl.textContent = myName();
  const who = document.getElementById("whoami");
  if (who) who.addEventListener("click", (e) => {
    e.preventDefault();
    const input = el("input", { value: myName() === "anonymous" ? "" : myName(), placeholder: "Your name", className: "inline-input" });
    nameEl.replaceWith(input);
    input.focus();
    const save = () => fetch("/api/me", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ name: input.value }) })
      .then(() => location.reload());
    input.addEventListener("keydown", (k) => { if (k.key === "Enter") save(); if (k.key === "Escape") location.reload(); });
    input.addEventListener("blur", save);
  });

  // ---- notes panels -------------------------------------------------------
  async function api(body) {
    const r = await fetch("/api/notes", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify(body) });
    if (!r.ok) throw new Error(await r.text());
    return r.json();
  }
  function renderNotes(panel, state) {
    const me = myName();
    const buttons = panel.querySelector(".flag-buttons");
    buttons.innerHTML = "";
    FLAGS.forEach(([key, label, icon], idx) => {
      const who = state.flags[key] || [];
      const b = el("button", { type: "button", className: "flag-btn flag-" + key + (who.includes(me) ? " on" : ""), title: label + " (" + (idx + 1) + ")" });
      b.append(el("span", { className: "flag-icon" }, icon), document.createTextNode(" " + label));
      if (who.length) b.append(el("span", { className: "flag-count" }, String(who.length)));
      b.addEventListener("click", async () => renderNotes(panel, await api({ kind: state.kind, id: state.entity, flag: key })));
      buttons.append(b);
    });
    const list = panel.querySelector(".note-list");
    list.innerHTML = "";
    state.notes.slice().reverse().forEach((n) => {
      const item = el("div", { className: "note-item" });
      item.append(el("div", { className: "note-text" }, n.text));
      const meta = el("div", { className: "sub" }, n.author + " · " + n.created_at.replace("T", " ").slice(0, 16));
      if (n.author === me) {
        const del = el("a", { href: "#", className: "note-del" }, "delete");
        del.addEventListener("click", async (e) => {
          e.preventDefault();
          const r = await fetch("/api/notes/" + n.id + "/delete", { method: "POST" });
          if (r.ok) renderNotes(panel, await r.json());
        });
        meta.append(document.createTextNode(" · "), del);
      }
      item.append(meta);
      list.append(item);
    });
    panel.querySelector(".who").textContent = "Signed as " + me + ". Change the name in the left rail.";
    panel._state = state;
  }
  document.querySelectorAll("[data-notes-kind]").forEach(async (panel) => {
    const kind = panel.dataset.notesKind, id = panel.dataset.notesId;
    const form = panel.querySelector(".note-form"), area = form.querySelector("textarea");
    const submit = async (e) => {
      if (e) e.preventDefault();
      if (!area.value.trim()) return;
      renderNotes(panel, await api({ kind, id, text: area.value }));
      area.value = "";
    };
    form.addEventListener("submit", submit);
    area.addEventListener("keydown", (e) => { if (e.key === "Enter" && !e.shiftKey) { e.preventDefault(); submit(); } });
    const r = await fetch("/api/notes?kind=" + encodeURIComponent(kind) + "&id=" + encodeURIComponent(id));
    if (r.ok) renderNotes(panel, await r.json());
  });

  // compact flags on stacked loci rows
  document.querySelectorAll(".stack-item[data-entity]").forEach(async (row) => {
    const kind = row.dataset.entity, id = row.dataset.id;
    const holder = el("div", { className: "row-flags" });
    row.querySelector(".stack-label").append(holder);
    const draw = (state) => {
      holder.innerHTML = "";
      FLAGS.forEach(([key, label, icon]) => {
        const n = (state.flags[key] || []).length;
        const b = el("button", { type: "button", className: "mini-flag flag-" + key + ((state.flags[key] || []).includes(myName()) ? " on" : ""), title: label + (n ? " · " + n : "") }, icon);
        b.addEventListener("click", async () => draw(await api({ kind, id, flag: key })));
        holder.append(b);
      });
    };
    const r = await fetch("/api/notes?kind=" + kind + "&id=" + encodeURIComponent(id));
    if (r.ok) draw(await r.json());
  });

  // ---- collection ---------------------------------------------------------
  const cart = () => store.get("sharur_collection", { protein: [], genome: [] });
  const badge = () => {
    const c = cart(), n = c.protein.length + c.genome.length;
    const b = document.getElementById("collection-badge");
    if (b) b.textContent = n ? String(n) : "";
  };
  badge();
  document.querySelectorAll("[data-collect]").forEach((b) => {
    const kind = b.dataset.collect, id = b.dataset.id;
    const sync = () => { const has = cart()[kind].includes(id); b.textContent = has ? "✓ Collected" : "＋ Collect"; b.classList.toggle("done", has); };
    sync();
    b.addEventListener("click", () => {
      const c = cart();
      c[kind] = c[kind].includes(id) ? c[kind].filter((x) => x !== id) : c[kind].concat([id]);
      store.set("sharur_collection", c);
      sync(); badge();
    });
  });
  const collectionTable = document.getElementById("collection");
  if (collectionTable) {
    const c = cart();
    const total = c.protein.length + c.genome.length;
    document.getElementById("collection-count").textContent = total ? total + " items" : "";
    document.getElementById("collection-empty").hidden = total > 0;
    document.querySelector("#fasta-form input[name=ids]").value = c.protein.join("\n");
    document.getElementById("fasta-form").hidden = !c.protein.length;
    const prompt = document.getElementById("agent-prompt");
    prompt.value = "Here is a set of items from the Sharur dataset I am browsing.\n" +
      (c.protein.length ? "\nProteins:\n" + c.protein.map((p) => "- " + p).join("\n") + "\n" : "") +
      (c.genome.length ? "\nGenomes:\n" + c.genome.map((g) => "- " + g).join("\n") + "\n" : "") +
      "\nStart with `sharur card <protein>` and `sharur why <protein> <label>` for proteins, and `sharur describe` for the dataset. " +
      "Report observed domains separately from named functions.\n";
    document.getElementById("copy-prompt").addEventListener("click", async (e) => {
      const b = e.currentTarget;
      try { await Promise.race([navigator.clipboard.writeText(prompt.value), new Promise((_, r) => setTimeout(() => r(new Error("t")), 700))]); b.textContent = "Copied ✓"; }
      catch (err) { prompt.classList.remove("visually-hidden"); prompt.select(); const ok = document.execCommand("copy"); prompt.classList.add("visually-hidden"); b.textContent = ok ? "Copied ✓" : "Copy failed"; }
      setTimeout(() => { b.textContent = b.dataset.label; }, 1600);
    });
    // collected genomes as one group for the matrix or a comparison
    const asGroup = (go) => fetch("/api/selection", { method: "POST", headers: { "Content-Type": "application/json" },
      body: JSON.stringify({ genomes: c.genome }) }).then((r) => r.json()).then((d) => { if (d.token) go(encodeURIComponent(d.token)); });
    const mx = document.getElementById("collection-matrix"), cmp = document.getElementById("collection-compare");
    if (mx && c.genome.length) { mx.hidden = false; mx.addEventListener("click", () => asGroup((t) => { location.href = "/matrix?a=" + t; })); }
    if (cmp && c.genome.length > 1) { cmp.hidden = false; cmp.addEventListener("click", () => asGroup((t) => { location.href = "/compare?a=" + t; })); }
    document.getElementById("clear-collection").addEventListener("click", () => { store.set("sharur_collection", { protein: [], genome: [] }); location.reload(); });
    if (total) fetch("/api/collection", { method: "POST", headers: { "Content-Type": "application/json" }, body: JSON.stringify({ proteins: c.protein, genomes: c.genome }) })
      .then((r) => r.json()).then((rows) => {
        const body = collectionTable.tBodies[0];
        rows.forEach((r) => {
          const tr = el("tr");
          const a = el("a", { href: r.url, className: "mono" }, r.id.split("|").pop());
          const td1 = el("td"); td1.append(el("span", { className: "tag" }, r.kind), document.createTextNode(" "), a);
          const td2 = el("td"); td2.append(el("span", {}, r.lineage));
          const td3 = el("td", { className: "num" }, r.kind === "protein" ? (r.length || 0).toLocaleString() + " aa" : (r.length / 1e6).toFixed(2) + " Mb");
          const td4 = el("td", {}, r.best);
          const td5 = el("td"); const x = el("a", { href: "#" }, "remove");
          x.addEventListener("click", (e) => { e.preventDefault(); const cc = cart(); cc[r.kind] = cc[r.kind].filter((i) => i !== r.id); store.set("sharur_collection", cc); location.reload(); });
          td5.append(x);
          tr.append(td1, td2, td3, td4, td5);
          body.append(tr);
        });
        collectionTable.dispatchEvent(new Event("rows-sorted"));
        document.dispatchEvent(new Event("tables-changed"));
      });
  }

  // ---- recently viewed ----------------------------------------------------
  const recent = store.get("sharur_recent", []);
  const here = { url: location.pathname + location.search, title: document.title.replace(/ · Sharur$/, "") };
  if (!/^\/(search|api|flags|collection|triage)/.test(location.pathname)) {
    store.set("sharur_recent", [here].concat(recent.filter((r) => r.url !== here.url)).slice(0, 8));
  }
  const recentList = document.getElementById("recent-list");
  if (recentList) {
    recent.filter((r) => r.url !== here.url).slice(0, 6).forEach((r) => recentList.append(el("a", { href: r.url, title: r.title }, r.title)));
    recentList.parentElement.hidden = !recentList.children.length;
  }

  // ---- hover previews -----------------------------------------------------
  const card = el("div", { className: "preview" });
  card.hidden = true;
  document.body.append(card);
  const cache = new Map();
  let timer = null, current = null;
  document.addEventListener("mouseover", (e) => {
    const a = e.target.closest && e.target.closest('a[href^="/protein/"], a[href^="/genome/"]');
    if (!a || a === current || a.closest(".preview") || a.closest("svg")) return;
    current = a;
    clearTimeout(timer);
    timer = setTimeout(async () => {
      const href = a.getAttribute("href");
      if (!cache.has(href)) {
        const r = await fetch("/api/preview?href=" + encodeURIComponent(href));
        cache.set(href, r.ok ? await r.text() : "");
      }
      if (current !== a || !cache.get(href)) return;
      card.innerHTML = cache.get(href);
      const box = a.getBoundingClientRect();
      card.hidden = false;
      const left = Math.min(window.innerWidth - card.offsetWidth - 12, Math.max(12, box.left));
      const top = box.bottom + 8 + card.offsetHeight > window.innerHeight ? box.top - card.offsetHeight - 8 : box.bottom + 8;
      card.style.left = left + window.scrollX + "px";
      card.style.top = top + window.scrollY + "px";
    }, 350);
  });
  document.addEventListener("mouseout", (e) => {
    if (current && (!e.relatedTarget || !current.contains(e.relatedTarget))) { clearTimeout(timer); current = null; card.hidden = true; }
  });

  // ---- tables: TSV export and j/k selection --------------------------------
  function decorateTables() {
    document.querySelectorAll("table[data-table]").forEach((table) => {
      if (table.dataset.tsv || table.dataset.noTsv !== undefined) return;
      table.dataset.tsv = "1";
      const tools = table.closest(".panel") && table.closest(".panel").querySelector(".table-tools");
      if (!tools) return;
      const b = el("button", { type: "button", className: "tsv-btn", title: "Download the rows shown as TSV" }, "TSV");
      b.addEventListener("click", () => {
        const head = Array.from(table.tHead ? table.tHead.rows[0].cells : []).map((c) => c.textContent.trim());
        const q = (document.querySelector(`input[data-filter="${table.id}"]`) || { value: "" }).value.trim().toLowerCase();
        const rows = Array.from(table.tBodies[0].rows).filter((r) => !q || r.textContent.toLowerCase().includes(q))
          .map((r) => Array.from(r.cells).map((c) => (c.dataset.v ?? c.textContent).trim().replace(/\s+/g, " ")).join("\t"));
        const blob = new Blob([[head.join("\t")].concat(rows).join("\n") + "\n"], { type: "text/tab-separated-values" });
        const link = el("a", { href: URL.createObjectURL(blob), download: (document.title.split(" · ")[0] || "table").replace(/[^\w.-]+/g, "_") + ".tsv" });
        document.body.append(link); link.click(); link.remove();
      });
      tools.append(b);
    });
  }
  decorateTables();
  document.addEventListener("tables-changed", decorateTables);
  let selected = -1;
  function visibleRows() {
    const t = document.querySelector("table[data-table]");
    return t ? Array.from(t.tBodies[0].rows).filter((r) => r.style.display !== "none") : [];
  }

  // ---- keyboard ------------------------------------------------------------
  const notesPanel = document.querySelector("[data-notes-kind]");
  document.addEventListener("keydown", (e) => {
    if (typing() || e.metaKey || e.ctrlKey || e.altKey) return;
    const go = (key) => { const a = document.querySelector('[data-key="' + key + '"]'); if (a && !a.classList.contains("off")) { location.href = a.href; return true; } return false; };
    if (e.key === "?") { const h = document.getElementById("help"); h.hidden = !h.hidden; return; }
    if (e.key === "Escape") { document.getElementById("help").hidden = true; return; }
    if (e.key === "[") { go("prev-gene"); return; }
    if (e.key === "]") { go("next-gene"); return; }
    if (e.key === "c") { const b = document.querySelector("[data-collect]"); if (b) b.click(); return; }
    if (notesPanel && /^[1-4]$/.test(e.key)) { const b = notesPanel.querySelectorAll(".flag-btn")[Number(e.key) - 1]; if (b) b.click(); return; }
    if (notesPanel && e.key === "n") { e.preventDefault(); notesPanel.querySelector("textarea").focus(); return; }
    if (e.key === "j" || e.key === "k") {
      if (go(e.key === "j" ? "next" : "prev")) return;
      const rows = visibleRows();
      if (!rows.length) return;
      rows.forEach((r) => r.classList.remove("kbd-sel"));
      selected = Math.max(0, Math.min(rows.length - 1, selected + (e.key === "j" ? 1 : -1)));
      rows[selected].classList.add("kbd-sel");
      rows[selected].scrollIntoView({ block: "nearest" });
      return;
    }
    if (e.key === "Enter" || e.key === "o") {
      const open = document.querySelector(".triage-bar") && document.querySelector(".panel h2 a");
      if (e.key === "o" && open) { location.href = open.href; return; }
      const rows = visibleRows();
      if (selected >= 0 && rows[selected]) { const a = rows[selected].querySelector("a"); if (a) location.href = a.href; }
    }
  });

  // ---- figure export (SVG / PNG, light theme) --------------------------------
  const STYLE_PROPS = ["fill", "fill-opacity", "stroke", "stroke-width", "stroke-dasharray", "stroke-opacity", "opacity",
                       "font-family", "font-size", "font-weight", "text-anchor", "paint-order", "letter-spacing"];
  function inlineCopy(svg) {
    document.documentElement.classList.add("export-light");
    const clone = svg.cloneNode(true);
    const src = svg.querySelectorAll("*"), dst = clone.querySelectorAll("*");
    src.forEach((node, i) => {
      const cs = getComputedStyle(node);
      dst[i].setAttribute("style", STYLE_PROPS.map((p) => p + ":" + cs.getPropertyValue(p)).join(";"));
    });
    document.documentElement.classList.remove("export-light");
    clone.querySelectorAll("a").forEach((a) => { while (a.firstChild) a.parentNode.insertBefore(a.firstChild, a); a.remove(); });
    const vb = svg.viewBox.baseVal;
    clone.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    clone.setAttribute("width", vb.width); clone.setAttribute("height", vb.height);
    clone.removeAttribute("class");
    return clone;
  }
  function stackSvg(stack) {
    // rows become one figure: genome label column + the shared-scale tracks
    const rows = Array.from(stack.querySelectorAll(".stack-item"));
    const labelW = 260, rowH = 44, width = 1040 + labelW;
    const out = document.createElementNS("http://www.w3.org/2000/svg", "svg");
    out.setAttribute("xmlns", "http://www.w3.org/2000/svg");
    out.setAttribute("viewBox", "0 0 " + width + " " + rows.length * rowH);
    out.setAttribute("width", width); out.setAttribute("height", rows.length * rowH);
    rows.forEach((row, i) => {
      const t = document.createElementNS("http://www.w3.org/2000/svg", "text");
      t.setAttribute("x", 4); t.setAttribute("y", i * rowH + 26);
      t.setAttribute("style", "font-family:Helvetica,Arial,sans-serif;font-size:12px;fill:#161b22");
      t.textContent = (row.querySelector(".stack-label a").textContent + "  " + (row.querySelector(".stack-label .sub") || { textContent: "" }).textContent).slice(0, 48);
      out.append(t);
      const g = document.createElementNS("http://www.w3.org/2000/svg", "g");
      g.setAttribute("transform", "translate(" + labelW + "," + i * rowH + ")");
      Array.from(inlineCopy(row.querySelector("svg")).childNodes).forEach((n) => g.append(n));
      out.append(g);
    });
    return out;
  }
  function download(svgNode, name, png) {
    const text = new XMLSerializer().serializeToString(svgNode);
    const fname = (document.title.split(" · ")[0] + "_" + name).replace(/[^\w.-]+/g, "_");
    if (!png) {
      const link = el("a", { href: URL.createObjectURL(new Blob([text], { type: "image/svg+xml" })), download: fname + ".svg" });
      document.body.append(link); link.click(); link.remove(); return;
    }
    const img = new Image();
    img.onload = () => {
      const scale = 3, canvas = el("canvas");
      canvas.width = svgNode.getAttribute("width") * scale; canvas.height = svgNode.getAttribute("height") * scale;
      const ctx = canvas.getContext("2d");
      ctx.fillStyle = "#ffffff"; ctx.fillRect(0, 0, canvas.width, canvas.height);
      ctx.drawImage(img, 0, 0, canvas.width, canvas.height);
      canvas.toBlob((blob) => { const link = el("a", { href: URL.createObjectURL(blob), download: fname + ".png" }); document.body.append(link); link.click(); link.remove(); });
    };
    img.src = "data:image/svg+xml;charset=utf-8," + encodeURIComponent(text);
  }
  function addExport(container, make, name) {
    const bar = el("div", { className: "fig-export" });
    ["SVG", "PNG"].forEach((fmt) => {
      const b = el("button", { type: "button", title: "Download this figure as " + fmt + " (light theme)" }, "⤓ " + fmt);
      b.addEventListener("click", () => download(make(), name, fmt === "PNG"));
      bar.append(b);
    });
    container.prepend(bar);
  }
  document.querySelectorAll("svg.track, svg.hood, svg.contig-track, svg.hist, svg.strip, svg.genome-ring, svg.cc-spectrum, svg.tree-fig, svg.dotplot, svg.pw-diagram, svg.pw-heat").forEach((svg) => {
    const panel = svg.closest(".panel");
    if (panel && !panel.querySelector(":scope > .fig-export")) addExport(panel, () => inlineCopy(svg), svg.getAttribute("aria-label") || "figure");
  });
  document.querySelectorAll(".stack").forEach((stack) => {
    if (stack.querySelector(".stack-item")) addExport(stack.closest(".panel"), () => stackSvg(stack), "loci");
  });
})();
