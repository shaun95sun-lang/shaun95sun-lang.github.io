/* Archive side card: click (or Enter/Space on) a list item to open its detail card on the right.
   The detail lives in each item's hidden .more block, so the page still reads fine without JS. */
(function () {
  var drawer = document.getElementById('drawer');
  if (!drawer) return;
  var items = document.querySelectorAll('.list li');
  var xBtn = document.getElementById('drawer-x');
  var hEl = document.getElementById('drawer-h'), catEl = document.getElementById('drawer-cat');
  var stampEl = document.getElementById('drawer-stamp'), moreEl = document.getElementById('drawer-more'), tEl = document.getElementById('drawer-t');
  var current = null;

  function open(li, focus) {
    if (current) current.setAttribute('aria-expanded', 'false');
    current = li;
    li.setAttribute('aria-expanded', 'true');
    var group = li.closest('.group');
    catEl.textContent = group ? group.querySelector('h3').textContent : '專案';
    hEl.textContent = li.querySelector('b').textContent;
    var st = li.querySelector(':scope > .stamp');
    stampEl.innerHTML = st ? st.outerHTML : '';
    var more = li.querySelector('.more');
    moreEl.innerHTML = more ? more.innerHTML : '<p>' + li.querySelector('.d').textContent + '</p>';
    var t = li.querySelector('.t');
    tEl.textContent = t ? t.textContent : '';
    drawer.hidden = false;
    document.body.classList.add('drawer-on');
    drawer.querySelector('.drawer-body').scrollTop = 0;
    if (focus) xBtn.focus();
    try { history.replaceState(null, '', '#' + li.id); } catch (e) {}
  }

  function close() {
    if (drawer.hidden) return;
    drawer.hidden = true;
    document.body.classList.remove('drawer-on');
    if (current) { current.setAttribute('aria-expanded', 'false'); current.focus(); }
    current = null;
    try { history.replaceState(null, '', location.pathname + location.search); } catch (e) {}
  }

  items.forEach(function (li) {
    if (!li.id) return;
    li.tabIndex = 0;
    li.setAttribute('role', 'button');
    li.setAttribute('aria-controls', 'drawer');
    li.setAttribute('aria-expanded', 'false');
    li.addEventListener('click', function (e) {
      if (e.target.closest('a')) return;
      if (current === li && !drawer.hidden) close(); else open(li, false);
    });
    li.addEventListener('keydown', function (e) {
      if (e.key === 'Enter' || e.key === ' ') { e.preventDefault(); open(li, true); }
    });
  });

  xBtn.addEventListener('click', close);
  document.addEventListener('keydown', function (e) { if (e.key === 'Escape') close(); });
  document.addEventListener('click', function (e) {
    if (!drawer.hidden && !drawer.contains(e.target) && !e.target.closest('.list li')) close();
  });

  // links like work.html#omni (from the About page) open that card directly
  var target = location.hash && document.getElementById(location.hash.slice(1));
  if (target && target.matches('.list li')) open(target, false);
})();
