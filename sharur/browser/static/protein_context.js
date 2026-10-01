/* Protein page: fetch the context panels (architecture elsewhere, copies, usual neighbours) after the page renders. */
(function () {
  "use strict";
  const host = document.getElementById("protein-context");
  if (!host) return;
  fetch(host.dataset.contextUrl, { credentials: "same-origin" })
    .then((r) => (r.ok ? r.text() : Promise.reject(r.status)))
    .then((html) => { host.innerHTML = html; })
    .catch(() => {
      host.innerHTML = '<div class="panel"><p class="empty">This protein\'s context could not be loaded.</p></div>';
    });
})();
