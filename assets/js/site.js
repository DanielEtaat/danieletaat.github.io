(function () {
  var menuButton = document.querySelector('[data-menu-toggle]');
  var menu = document.querySelector('[data-menu]');

  if (menuButton && menu) {
    var compactMenu = window.matchMedia('(max-width: 820px)');
    var setMenuOpen = function (open) {
      menuButton.setAttribute('aria-expanded', String(open));
      menu.toggleAttribute('data-open', open);
      menu.inert = compactMenu.matches && !open;
    };
    setMenuOpen(false);
    menuButton.addEventListener('click', function () {
      var open = menuButton.getAttribute('aria-expanded') === 'true';
      setMenuOpen(!open);
    });
    menu.addEventListener('click', function (event) {
      if (event.target.closest('a')) setMenuOpen(false);
    });
    compactMenu.addEventListener('change', function () {
      setMenuOpen(false);
    });
    document.addEventListener('keydown', function (event) {
      if (event.key === 'Escape' && menuButton.getAttribute('aria-expanded') === 'true') {
        setMenuOpen(false);
        menuButton.focus();
      }
    });
  }

  var progress = document.querySelector('[data-reading-progress]');
  if (progress) {
    var updateProgress = function () {
      var doc = document.documentElement;
      var available = doc.scrollHeight - doc.clientHeight;
      var value = available > 0 ? (doc.scrollTop / available) * 100 : 0;
      progress.style.setProperty('--reading-progress', value + '%');
    };
    updateProgress();
    document.addEventListener('scroll', updateProgress, { passive: true });
  }

  var themeButton = document.querySelector('[data-theme-toggle]');
  if (themeButton) {
    var themeTransitionTimer;
    var themeColor = document.querySelector('meta[name="theme-color"]');
    var syncThemeColor = function () {
      if (!themeColor) return;
      var background = window.getComputedStyle(document.documentElement).getPropertyValue('--bg').trim();
      if (background) themeColor.setAttribute('content', background);
    };
    var storedTheme = null;
    try { storedTheme = localStorage.getItem('daniel-theme'); } catch (error) {}
    if (storedTheme) document.documentElement.dataset.theme = storedTheme;
    var initialTheme = document.documentElement.dataset.theme || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
    themeButton.setAttribute('aria-pressed', String(initialTheme === 'dark'));
    themeButton.setAttribute('aria-label', 'Use ' + (initialTheme === 'light' ? 'dark' : 'light') + ' theme');
    syncThemeColor();
    themeButton.addEventListener('click', function () {
      var current = document.documentElement.dataset.theme || (window.matchMedia('(prefers-color-scheme: dark)').matches ? 'dark' : 'light');
      var next = current === 'light' ? 'dark' : 'light';
      window.clearTimeout(themeTransitionTimer);
      if (!window.matchMedia('(prefers-reduced-motion: reduce)').matches) {
        document.documentElement.classList.add('theme-changing');
        themeTransitionTimer = window.setTimeout(function () {
          document.documentElement.classList.remove('theme-changing');
        }, 450);
      }
      document.documentElement.dataset.theme = next;
      try { localStorage.setItem('daniel-theme', next); } catch (error) {}
      themeButton.setAttribute('aria-label', 'Use ' + (next === 'light' ? 'dark' : 'light') + ' theme');
      themeButton.setAttribute('aria-pressed', String(next === 'dark'));
      syncThemeColor();
    });
  }

  var toc = document.querySelector('[data-toc]');
  var article = document.querySelector('[data-article-body]');
  if (toc && article) {
    var headings = Array.prototype.slice.call(article.querySelectorAll('h2'));
    if (!headings.length) {
      toc.hidden = true;
    } else {
      var list = document.createElement('ol');
      headings.forEach(function (heading) {
        if (!heading.id) heading.id = heading.textContent.toLowerCase().replace(/[^a-z0-9]+/g, '-').replace(/(^-|-$)/g, '');
        var item = document.createElement('li');
        var link = document.createElement('a');
        link.href = '#' + heading.id;
        link.textContent = heading.textContent;
        item.appendChild(link);
        list.appendChild(item);
      });
      toc.appendChild(list);
    }
  }

  var responsiveToc = document.querySelector('[data-responsive-toc]');
  if (responsiveToc && window.matchMedia('(max-width: 820px)').matches) {
    responsiveToc.removeAttribute('open');
  }
})();
