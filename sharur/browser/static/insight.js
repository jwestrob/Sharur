// Sharur browser: similar proteins, structure viewer, finding verification.
(function () {
  "use strict";

  function el(tag, attrs, text) {
    const node = document.createElement(tag);
    Object.entries(attrs || {}).forEach(([k, v]) => node.setAttribute(k, v));
    if (text !== undefined && text !== null) node.textContent = text;
    return node;
  }

  // ---- similar proteins ---------------------------------------------------
  const similar = document.querySelector("[data-similar-url]");
  const status = similar && similar.querySelector("[data-similar-status]");
  if (similar && status) {
    fetch(similar.dataset.similarUrl).then((r) => r.json()).then((data) => {
      const table = similar.querySelector("[data-similar-table]");
      if (!data.neighbors || !data.neighbors.length) {
        status.textContent = data.state === "available" ? "No neighbors found for this protein." :
          "Similarity search is unavailable (" + (data.detail || data.state) + ").";
        return;
      }
      const body = table.tBodies[0];
      data.neighbors.forEach((n) => {
        const row = el("tr");
        const p = el("td");
        p.appendChild(el("a", {href: n.url, class: "mono", title: n.protein_id}, n.protein_id.split("|").pop()));
        const g = el("td");
        if (n.genome_url) g.appendChild(el("a", {href: n.genome_url}, n.lineage || n.bin_id));
        const len = el("td", {class: "num"}, n.length ? n.length.toLocaleString() + " aa" : "–");
        const what = el("td");
        what.appendChild(el("div", {}, n.architecture || "no placed domains"));
        if (n.best_hit) what.appendChild(el("div", {class: "sub"}, n.best_hit));
        const sim = el("td", {class: "num"});
        const bar = el("div", {class: "cell-bar"});
        const track = el("div", {class: "bar"});
        const fill = el("i");
        fill.style.width = Math.max(0, Math.min(1, n.similarity)) * 100 + "%";
        track.appendChild(fill);
        bar.append(track, el("span", {class: "count"}, n.similarity.toFixed(3)));
        sim.appendChild(bar);
        row.append(p, g, len, what, sim);
        body.appendChild(row);
      });
      table.hidden = false;
      status.textContent = "Cosine similarity of mean-pooled embeddings (values crowd near 1, so rank matters more " +
        "than the absolute score); search took " + data.ms + " ms.";
      status.className = "note";
    }).catch(() => { status.textContent = "Similarity search failed to load."; });
  }

  // ---- structure viewer ---------------------------------------------------
  const mol = document.getElementById("mol");
  if (mol) {
    const molStatus = mol.querySelector("[data-mol-status]");
    let viewer = null;
    let bScale = 1;  // pLDDT in the B-factor column: 0-1 (ESM3) or 0-100 (AlphaFold)
    function plddtColor(atom) {
      const b = atom.b * bScale;
      return b > 90 ? "#0053d6" : b > 70 ? "#65cbf3" : b > 50 ? "#ffdb13" : "#ff7d45";
    }
    function show(key) {
      fetch("/structure-file/" + encodeURIComponent(key)).then((r) => r.text()).then((pdb) => {
        if (!viewer) {
          if (molStatus) molStatus.remove();
          viewer = window.$3Dmol.createViewer(mol, {backgroundColor: "rgba(0,0,0,0)"});
        }
        viewer.clear();
        const model = viewer.addModel(pdb, "pdb");
        const maxB = model.selectedAtoms({}).reduce((m, atom) => Math.max(m, atom.b), 0);
        bScale = maxB <= 1.0 ? 100 : 1;
        viewer.setStyle({}, {cartoon: {colorfunc: plddtColor}});
        viewer.zoomTo();
        viewer.render();
      }).catch(() => { if (molStatus) molStatus.textContent = "Could not load the structure file."; });
    }
    const script = el("script", {src: "https://cdn.jsdelivr.net/npm/3dmol@2.5.5/build/3Dmol-min.js"});
    script.onload = () => show(mol.dataset.first);
    script.onerror = () => { if (molStatus) molStatus.textContent = "The 3D viewer library could not be loaded (offline?)."; };
    document.head.appendChild(script);
    document.querySelectorAll("[data-structure]").forEach((button) => {
      button.addEventListener("click", () => { if (window.$3Dmol) show(button.dataset.structure); });
    });
  }

  // ---- finding verification -----------------------------------------------
  const verify = document.querySelector("[data-verify-url]");
  const button = verify && verify.querySelector("[data-verify]");
  if (verify && button) {
    button.addEventListener("click", () => {
      button.disabled = true;
      button.textContent = "Running…";
      fetch(verify.dataset.verifyUrl, {method: "POST"}).then((r) => r.json()).then((data) => {
        let pass = 0, fail = 0;
        (data.results || []).forEach((result, i) => {
          const row = verify.querySelector(`[data-check="${i}"] .check-result`);
          if (!row) return;
          row.textContent = "";
          const cls = {pass: "ok", fail: "warn", error: "warn"}[result.status] || "";
          row.appendChild(el("span", {class: "tag " + cls}, result.status));
          const detail = result.detail || (result.observed !== undefined ? "observed " + JSON.stringify(result.observed) : "");
          row.appendChild(el("span", {class: "count"}, " " + detail + (result.ms !== undefined ? " · " + result.ms + " ms" : "")));
          if (result.status === "pass") pass++;
          if (result.status === "fail" || result.status === "error") fail++;
        });
        button.textContent = `Re-run checks (${pass} passed, ${fail} failed)`;
      }).catch(() => { button.textContent = "Check run failed"; }).finally(() => { button.disabled = false; });
    });
  }
})();
