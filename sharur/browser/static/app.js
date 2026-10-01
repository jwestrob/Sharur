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
      });
    });
    const filter = document.querySelector(`input[data-filter="${table.id}"]`);
    const counter = document.querySelector(`[data-count="${table.id}"]`);
    function update() {
      const q = filter ? filter.value.trim().toLowerCase() : "";
      let shown = 0;
      Array.from(body.rows).forEach((r) => {
        const hit = !q || r.textContent.toLowerCase().includes(q);
        r.style.display = hit ? "" : "none";
        if (hit) shown++;
      });
      if (counter) counter.textContent = shown.toLocaleString() + " of " + body.rows.length.toLocaleString();
    }
    if (filter) filter.addEventListener("input", update);
    update();
  });
})();
