'use strict';
(() => {
  const data = JSON.parse(document.getElementById('results-data').textContent);
  const items = data.items;
  const byId = new Map(items.map(item => [item.id, item]));
  const $ = id => document.getElementById(id);
  const number = value => Number(value || 0).toLocaleString('en-US');
  const labels = {ready: 'Ready', reference: 'Reference data only', partial: 'Needs attention', unavailable: 'No current plots'};
  const state = {view: 'analyses', page: 0, measurement: ''};
  const PAGE_SIZE = 48;
  const url = path => path.split('/').map(encodeURIComponent).join('/');
  const el = (tag, text, cls) => {const node = document.createElement(tag); if (text !== undefined && text !== null) node.textContent = text; if (cls) node.className = cls; return node;};
  const button = (text, action, cls = 'secondary') => {const node = el('button', text, cls); node.type = 'button'; node.addEventListener('click', action); return node;};
  const link = (text, href, external = false) => {const node = el('a', text); node.href = external ? href : url(href); if (external) {node.target = '_blank'; node.rel = 'noopener noreferrer';} return node;};
  const pretty = value => JSON.stringify(value, null, 2);
  const date = value => value ? value.replace('T', ' ').replace(/\+00:00$/, ' UTC').replace(/Z$/, ' UTC') : 'Not recorded';
  const badge = item => el('span', labels[item.availability], 'badge ' + item.availability);

  $('created').textContent = date(data.created_at);
  const plotTotal = items.reduce((sum, i) => sum + i.plots.length, 0);
  const ready = items.filter(i => i.availability === 'ready').length;
  [[items.length, 'registered analyses'], [plotTotal, 'available plots'], [ready, 'ready to browse'], [items.length - ready, 'need attention']].forEach(([value, label], index) => {
    const node = el('div', null, 'metric' + (index === 3 ? ' attention' : ''));
    node.append(el('strong', number(value)), el('span', label)); $('metrics').append(node);
  });
  [...new Set(items.map(i => i.experiment))].sort().forEach(experiment => {const option = el('option', experiment); option.value = experiment; $('experiment').append(option);});
  const allPlots = items.flatMap(item => item.plots.map(plot => ({item, plot})));
  const searchText = item => [item.id, item.title, item.data_type, item.authority, item.selected?.tag, item.latest?.tag].join(' ').toLowerCase();
  const matches = item => (!$('experiment').value || item.experiment === $('experiment').value) && (!$('availability').value || item.availability === $('availability').value) && (!state.measurement || item.id === state.measurement);
  const search = () => $('search').value.trim().toLowerCase();
  const setView = view => {state.view = view; state.page = 0; render();};
  function showPlots(id) {state.measurement = id; $('search').value = ''; $('experiment').value = ''; $('availability').value = ''; $('info-dialog').close(); setView('plots'); $('result-title').scrollIntoView({behavior:'smooth', block:'start'});}

  function card(item) {
    const node = el('article', null, 'analysis-card');
    const top = el('div', null, 'card-top'); top.append(el('span', item.experiment, 'experiment'), badge(item));
    node.append(top, el('h3', item.data_type || item.title), el('div', item.id, 'identifier'));
    const ref = item.data_available ? `${item.authority} · ${number(item.reference_entries)} reference entries` : 'Generator diagnostic · no experimental reference';
    node.append(el('p', ref, 'card-summary'));
    const run = item.selected || item.latest;
    const details = el('div', null, 'run-line');
    if (run) {
      details.append(el('strong', item.selected ? 'Displayed campaign' : 'Latest recorded campaign'), el('span', run.tag, 'tag'));
      details.append(el('span', `${number(run.successful_shards)} / ${number(run.total_shards)} shards complete · ${number(item.plots.length)} plots`));
    } else if (item.availability === 'reference') {
      details.append(el('strong', `${number(item.plots.length)} published-data plots`), el('span', 'Simulation is not implemented for this measurement.'));
    } else details.append(el('span', 'No production campaign found in the scanned folders.'));
    if (item.latest && (!item.selected || item.latest.path !== item.selected.path)) details.append(el('div', item.latest.activity, 'small-status'));
    if (!item.selected && item.latest && !item.latest.compatible) details.append(el('div', 'Source definitions differ; historical plots are withheld.', 'small-status'));
    if (item.diagnostics.total_masked_bins !== undefined && item.diagnostics.total_masked_bins !== null) details.append(el('div', `${number(item.diagnostics.total_masked_bins)} masked bins`, 'identifier'));
    const actions = el('div', null, 'card-actions');
    const plots = button(`View ${number(item.plots.length)} plots`, () => showPlots(item.id), 'primary'); plots.disabled = item.plots.length === 0;
    actions.append(plots, button('Analysis info', () => info(item.id)));
    node.append(details, actions); return node;
  }

  function plotCard({item, plot}) {
    const node = el('figure', null, 'plot-card');
    if (plot.files.png) {
      const anchor = link('', plot.files.png); anchor.className = 'image-link'; anchor.target = '_blank'; anchor.rel = 'noopener';
      const image = el('img'); image.src = url(plot.files.png); image.loading = 'lazy'; image.alt = `${item.id}: ${plot.name}`;
      anchor.append(image); node.append(anchor);
    } else node.append(el('p', 'PDF plot — use the download link below.', 'empty'));
    const caption = el('figcaption'); caption.append(el('span', item.experiment + ' · ' + item.id, 'experiment'));
    caption.append(el('div', plot.name.split('/').pop(), 'plot-name'));
    const description = plot.labels.Title || [plot.labels.XLabel, plot.labels.YLabel].filter(Boolean).join(' → ');
    if (description) caption.append(el('div', description, 'plot-label'));
    const downloads = el('div', null, 'plot-links');
    Object.entries(plot.files).forEach(([format, path]) => {const anchor = link(format.toUpperCase(), path); anchor.target = '_blank'; anchor.rel = 'noopener'; downloads.append(anchor);});
    downloads.append(button('Analysis & run info', () => info(item.id), 'text-button'));
    caption.append(downloads); node.append(caption); return node;
  }

  function render() {
    $('analyses-tab').setAttribute('aria-pressed', String(state.view === 'analyses'));
    $('plots-tab').setAttribute('aria-pressed', String(state.view === 'plots'));
    $('results').replaceChildren(); $('pagination').replaceChildren(); $('active-filter').replaceChildren();
    if (state.measurement) {const chip = el('div', state.measurement, 'filter-chip'); chip.append(button('×', () => {state.measurement = ''; state.page = 0; render();}, 'text-button')); $('active-filter').append(chip);}
    const query = search();
    let filtered;
    if (state.view === 'analyses') {
      filtered = items.filter(item => matches(item) && (!query || searchText(item).includes(query)));
      $('results').className = 'analysis-grid'; $('result-title').textContent = 'Analysis overview';
      filtered.forEach(item => $('results').append(card(item)));
      $('result-count').textContent = `${number(filtered.length)} of ${number(items.length)} analyses`;
    } else {
      filtered = allPlots.filter(({item, plot}) => matches(item) && (!query || (searchText(item) + ' ' + plot.name + ' ' + Object.values(plot.labels).join(' ')).toLowerCase().includes(query)));
      const pages = Math.max(1, Math.ceil(filtered.length / PAGE_SIZE)); state.page = Math.min(state.page, pages - 1);
      $('results').className = 'plot-grid'; $('result-title').textContent = state.measurement ? byId.get(state.measurement).data_type : 'All available plots';
      filtered.slice(state.page * PAGE_SIZE, (state.page + 1) * PAGE_SIZE).forEach(value => $('results').append(plotCard(value)));
      $('result-count').textContent = `${number(filtered.length)} matching plots`;
      if (pages > 1) {
        const previous = button('← Previous', () => {state.page--; render(); $('result-title').scrollIntoView();}); previous.disabled = state.page === 0;
        const next = button('Next →', () => {state.page++; render(); $('result-title').scrollIntoView();}); next.disabled = state.page === pages - 1;
        $('pagination').append(previous, el('span', `Page ${state.page + 1} of ${pages}`), next);
      }
    }
    if (!filtered.length) $('results').append(el('div', 'No results match these filters. Try another search or reset the filters.', 'empty'));
  }

  function section(title, parent) {const node = el('section', null, 'info-section'); node.append(el('h3', title)); parent.append(node); return node;}
  function details(title, value, parent) {if (value === undefined || value === null || (typeof value === 'object' && !Object.keys(value).length)) return; const node = el('details'); node.append(el('summary', title), el('pre', pretty(value))); parent.append(node);}
  function info(id) {
    const item = byId.get(id); const body = $('info-content'); body.replaceChildren();
    const title = el('h2', item.title); title.id = 'info-title'; body.append(title, el('p', item.id, 'identifier'), badge(item));
    item.issues.forEach(message => body.append(el('p', message, 'notice')));
    const basics = section('Measurement and numerical source', body);
    basics.append(el('p', item.description || item.data_type), el('p', item.authority_note || item.authority));
    basics.append(el('p', item.data_available ? `${number(item.reference_entries)} reference entries in ${number(item.dataset_count)} groups. Overlapping projections are not independent datasets.` : 'Generator-only diagnostic: numerical experimental data are not supplied.'));
    const sources = el('div', null, 'info-links'); item.links.forEach(([label, href]) => sources.append(link(label + ' ↗', href, true))); Object.entries(item.downloads).forEach(([label, href]) => sources.append(link(label, href))); basics.append(sources);
    const run = item.selected || item.latest;
    const status = section('Campaign and result status', body);
    if (run) {
      const facts = el('div', null, 'info-facts');
      [[item.selected ? 'Displayed campaign' : 'Latest campaign', run.tag], ['Created', date(run.created_at)], ['Successful shards', `${number(run.successful_shards)} / ${number(run.total_shards)}`], ['Last recorded update', date(run.updated_at)], ['Latest campaign activity', item.latest?.activity], ['Source identity', run.compatibility_note]].forEach(([label, value]) => {const cell = el('div'); cell.append(el('strong', label), el('span', value || 'Not recorded')); facts.append(cell);});
      status.append(facts);
      details('Run configuration: families, helicities, PDFs, scales and statistics', run.configuration, status);
      details('Generator versions and source provenance', {versions:run.versions, source:run.source_control, campaign:run.path, manifest_sha256:run.manifest_sha256}, status);
    } else status.append(el('p', item.availability === 'reference' ? 'Published experimental points only. No Herwig prediction or Rivet event implementation is enabled.' : 'No production manifest was found for this analysis.'));
    details('Simulation status and requirements', item.simulation, status);
    const physics = section('Cuts, beams and physics definitions', body);
    const cuts = item.reference_metadata.selection;
    if (typeof cuts === 'string') physics.append(el('p', cuts)); else details('Implemented selection', cuts, physics);
    if (item.run_info) physics.append(el('p', item.run_info));
    details('Beam / target model and interpretation', item.physics, physics);
    details('Generator cuts and family definitions', {generator_cuts:item.generator_cuts, families:item.families}, physics);
    details('PDF ensembles and scale definitions', {pdf_ensembles:item.pdf_ensembles, scales:item.scales}, physics);
    details('Reference definition, binning, estimators and limitations', item.reference_metadata, physics);
    if (Object.keys(item.diagnostics).length) {
      const quality = section('Fit diagnostics and masks', body);
      if (item.diagnostics.total_masked_bins !== null) quality.append(el('p', `${number(item.diagnostics.total_masked_bins)} masked bins across the recorded datasets.`));
      details('Masked bins by dataset', item.diagnostics.masked_bins_by_dataset, quality);
      details('Recorded goodness of fit and uncertainty conventions', item.diagnostics.details, quality);
      quality.append(el('p', 'Values are read from the campaign summary. Large arrays are explicitly abbreviated; no new fit or covariance is calculated.'));
    }
    const history = section(`Campaign history · ${item.history.length} candidates`, body);
    if (item.history.length) {
      const wrap = el('div', null, 'table-scroll'), table = el('table', null, 'history-table'); const head = el('tr'); ['Campaign', 'Created', 'Shards', 'Definition'].forEach(value => head.append(el('th', value))); const thead = el('thead'); thead.append(head); table.append(thead); const tbody = el('tbody');
      item.history.forEach(c => {const row = el('tr'); [c.tag, date(c.created_at), `${number(c.successful_shards)} / ${number(c.total_shards)}`, c.compatible ? 'Compatible' : 'Different source'].forEach(value => row.append(el('td', value))); tbody.append(row);}); table.append(tbody); wrap.append(table); history.append(wrap);
    }
    details('Skipped candidates and scan issues', item.rejected, history);
    body.append(el('p', item.availability === 'reference' ? 'Reference-only registration; no compiled Rivet analysis is claimed.' : `Rivet metadata status: ${item.rivet_status}. Successful generation does not establish experimental fidelity or statistical convergence.`, 'method-note'));
    if (item.plots.length) body.append(button(`Browse ${number(item.plots.length)} plots`, () => showPlots(id), 'primary'));
    $('info-dialog').showModal(); $('info-dialog').scrollTop = 0;
  }

  $('close-info').addEventListener('click', () => $('info-dialog').close());
  $('info-dialog').addEventListener('click', event => {if (event.target === $('info-dialog')) {const r = event.target.getBoundingClientRect(); if (event.clientX < r.left || event.clientX > r.right || event.clientY < r.top || event.clientY > r.bottom) event.target.close();}});
  $('analyses-tab').addEventListener('click', () => setView('analyses'));
  $('plots-tab').addEventListener('click', () => setView('plots'));
  ['search', 'experiment', 'availability'].forEach(id => $(id).addEventListener(id === 'search' ? 'input' : 'change', () => {state.page = 0; render();}));
  $('reset').addEventListener('click', () => {$('search').value = ''; $('experiment').value = ''; $('availability').value = ''; state.measurement = ''; state.page = 0; render();});
  $('scan-policy').append(el('p', data.selection_policy), el('p', data.signature_policy), el('p', data.include_smoke ? 'Smoke tests were included in this snapshot.' : 'Recorded smoke tests were excluded.'), el('p', 'Campaign roots: ' + data.campaign_roots.join(' · ')));
  render();
})();
