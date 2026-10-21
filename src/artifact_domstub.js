/* A very small DOM, good enough to execute the atlas page under jsc.
   It is not a browser. It exists so that a wrong element id, a renamed data
   field or a runtime error is caught by a check rather than by someone opening
   the published page. Nesting and descendant selectors are tracked because the
   page relies on "nav.tabs button" meaning only the tab buttons. */
var ALL = [];
var VOID = {br:1, img:1, input:1, link:1, meta:1, hr:1, source:1, col:1};

function El(tag, cls, ds, id){
  this.tag = tag || "div";
  this.cls = cls || [];
  this.dataset = ds || {};
  this.id = id || "";
  this.parent = null; this.owner = null;
  this._html = ""; this.textContent = ""; this.attrs = {}; this.style = {};
  this.children = []; this.onclick = null; this.disabled = false;
  var self = this;
  this.classList = {
    toggle: function(c, on){
      var i = self.cls.indexOf(c);
      if (on === undefined) on = i < 0;
      if (on && i < 0) self.cls.push(c);
      if (!on && i >= 0) self.cls.splice(i, 1);
    },
    add: function(c){ if (self.cls.indexOf(c) < 0) self.cls.push(c); },
    remove: function(c){ var i = self.cls.indexOf(c); if (i >= 0) self.cls.splice(i, 1); },
    contains: function(c){ return self.cls.indexOf(c) >= 0; }
  };
}
Object.defineProperty(El.prototype, "innerHTML", {
  get: function(){ return this._html; },
  set: function(v){
    this._html = String(v);
    var self = this;
    ALL = ALL.filter(function(e){ return !isWithin(e, self); });
    parseInto(this._html, this);
  }
});
El.prototype.setAttribute = function(k, v){ this.attrs[k] = String(v); };
El.prototype.getAttribute = function(k){ return this.attrs[k]; };
El.prototype.appendChild = function(c){ this.children.push(c); c.parent = this; return c; };
El.prototype.click = function(){ if (this.onclick) this.onclick(); };
El.prototype.matches = function(sel){
  // one compound selector: tag, .class, or tag.class
  var m = /^([a-zA-Z0-9]*)((?:\.[A-Za-z0-9_-]+)*)$/.exec(sel);
  if (!m) return false;
  if (m[1] && this.tag !== m[1]) return false;
  var classes = m[2] ? m[2].split(".").filter(Boolean) : [];
  for (var i = 0; i < classes.length; i++)
    if (this.cls.indexOf(classes[i]) < 0) return false;
  return true;
};

function isWithin(el, ancestor){
  for (var p = el.parent; p; p = p.parent) if (p === ancestor) return true;
  return false;
}

/* "nav.tabs button" -> the element matches "button" AND has an ancestor
   matching "nav.tabs". Without this, every button on the page matches. */
function selectorMatch(el, sel){
  var parts = sel.trim().split(/\s+/);
  if (!el.matches(parts[parts.length - 1])) return false;
  var i = parts.length - 2, p = el.parent;
  while (i >= 0){
    while (p && !p.matches(parts[i])) p = p.parent;
    if (!p) return false;
    p = p.parent; i--;
  }
  return true;
}

function parseInto(html, owner){
  var re = /<(\/?)([a-zA-Z][a-zA-Z0-9]*)\b([^>]*?)(\/?)>/g, m;
  var stack = [];
  var top = function(){ return stack.length ? stack[stack.length - 1] : owner; };
  while ((m = re.exec(html)) !== null){
    var closing = m[1] === "/", tag = m[2].toLowerCase(), attrs = m[3],
        selfClose = m[4] === "/";
    if (closing){ if (stack.length) stack.pop(); continue; }
    var cm = /class="([^"]*)"/.exec(attrs);
    var im = /\bid="([^"]*)"/.exec(attrs);
    var cls = cm ? cm[1].split(/\s+/).filter(Boolean) : [];
    var ds = {}, dm, dre = /data-([a-zA-Z0-9-]+)="([^"]*)"/g;
    while ((dm = dre.exec(attrs)) !== null)
      ds[dm[1].replace(/-([a-z])/g, function(_, c){ return c.toUpperCase(); })] = dm[2];
    var el = new El(tag, cls, ds, im ? im[1] : "");
    el.parent = top(); el.owner = owner;
    var am, are = /([a-zA-Z-]+)="([^"]*)"/g;
    while ((am = are.exec(attrs)) !== null) el.attrs[am[1]] = am[2];
    ALL.push(el);
    if (!selfClose && !VOID[tag]) stack.push(el);
  }
}

var byId = {};
var document = {
  getElementById: function(id){
    if (byId[id]) return byId[id];
    var hit = ALL.filter(function(e){ return e.id === id; })[0];
    if (!hit) hit = new El("div", [], {}, id);
    byId[id] = hit;
    return hit;
  },
  querySelectorAll: function(sel){
    return ALL.filter(function(e){ return selectorMatch(e, sel); });
  },
  querySelector: function(sel){ return document.querySelectorAll(sel)[0] || null; },
  createElement: function(t){ return new El(t); }
};
El.prototype.querySelectorAll = function(sel){
  var self = this;
  return ALL.filter(function(e){ return isWithin(e, self) && selectorMatch(e, sel); });
};
El.prototype.querySelector = function(sel){ return this.querySelectorAll(sel)[0] || null; };

var window = { open: function(){} };
