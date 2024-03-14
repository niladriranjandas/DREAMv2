#!/usr/bin/env pyXplor

"""Generate Python documentation in HTML or text for interactive use using
modifications to the Python pydoc module.
"""

from os import environ
import collections
xplorDir=environ['XPLOR_DIR']
arch=environ['ARCH']
xplorPaths={}
xplorPaths['high'] = xplorDir+'/python/'
xplorPaths['low'] = xplorDir+'/bin.'+arch+'/'
includedModules=[]
pythonModules = ['__builtin__','xplor','repr','exceptions','glob']
# base for non-Xplor-NIH modules
import sys
pyDocURL="https://docs.python.org/%d/library/" % sys.version_info[0]
help_xplorDoc_dontCall=1

import pydoc
import inspect

# --------------------------------------------------------- common routines

initSysPath=sys.path # sys.path sometimes changed by loading modules
def pathdirs():
    """Convert sys.path into a list of absolute, existing, unique paths."""
    dirs = []
    normdirs = []
    for dir in initSysPath:
        dir = os.path.abspath(dir or '.')
        normdir = os.path.normcase(dir)
        if normdir not in normdirs and os.path.isdir(dir):
            dirs.append(dir)
            normdirs.append(normdir)
    return dirs

def getdoc(thing):
    """Get the doc string or comments for an object."""
    result=''
    if hasattr(thing,"pyXplorHelp") and not hasattr(thing,
                                                    "help_xplorDoc_dontCall"):
#        try:
        result = thing.help()
#        except:
#            pass
        pass
    result2 = inspect.getdoc(thing) or inspect.getcomments(thing) or ''
    result +=  '\n\n' + result2
    return result
    

import re
_re_stripid = re.compile(r' at 0x[0-9a-f]{6,16}(>+)$', re.IGNORECASE)
def stripid(text):
    """Remove the hexadecimal id from a Python object representation."""
    # The behaviour of %p is implementation-dependent in terms of case.
    return _re_stripid.sub(r'\1', text)


moduleBlacklist=[r'numarray\..*',
                 r'numpy\..*',
                 r'matplotlib\..*',
                 r'pyx\..*',
                 r'test\..*',
                 r'bzrlib\..*test.*',
                 r'version',
                 r'IPython\..*',
                 r'win32clip',
                 r'siteconfig',
                 r'bzr.*',
                 ]

import re
def xplorSynopsis(modname, cache={}, scanContents=False):
    """Get the one-line summary out of a Xplor-NIH module file."""
    for pattern in moduleBlacklist:
        if re.match(pattern,modname): return None
    try:
        module = safeimport(modname,cache)
    except pydoc.ErrorDuringImport:
        return None
    if scanContents:
        result = getdoc(module).strip()
    else:
        result = getdoc(module).strip().split('\n')[0]
        pass
    return result


import pkgutil, os

class HTMLDoc(pydoc.HTMLDoc):
    """Formatter class for Xplor-NIH HTML documentation."""

    # ------------------------------------------- HTML formatting utilities

    _repr_instance = pydoc.HTMLRepr()
    repr = _repr_instance.repr
    escape = _repr_instance.escape

    def classlink(self, object, modname):
        """Make a link for a class."""
        name, module = object.__name__, sys.modules.get(object.__module__)
        if hasattr(module, name) and getattr(module, name) is object:
            moduleLink=module.__name__
            if module.__name__ in pythonModules:
                moduleLink=pyDocURL + module.__name__
            return '<a href="%s.html#%s">%s</a>' % (
                module.__name__, name, pydoc.classname(object, modname))
        return pydoc.classname(object, modname)

    def modulelink(self, object):
        """Make a link for a module."""
        name = object.__name__
        url = name
#        if includedModules and name in pythonModules:
        if name in pythonModules:
            url = pyDocURL + name
            pass
        return '<a href="%s.html">%s</a>' % (url, name)

    def markup(self, text, escape=None, funcs={}, classes={}, methods={}):
        """Mark up some plain text, given a context of symbols to look for.
        some simple tags are defined:
           <\m NAME> - creates a link for Xplor-NIH module named NAME
           <\s NAME> - creates a link for Python module named NAME
           <\l LINK NAME> - creates a link LINK with tag NAME
        [In each case, please do not include the backslash character.]
        Each context dictionary maps object names to anchor names."""
        escape = escape or self.escape
        results = []
        here = 0
        pattern = re.compile(r'(<l\s+(\S+)\s+([^>]+)>|'
                             r'\b(http|ftp)://\S+[\w/]|'
                             r'\bRFC[- ]?(\d+)|'
                             r'\bPEP[- ]?(\d+)|'
                             r'\<m\s+(\S+\w)\>|'
                             r'\<s\s+(\S+\w)\>|'
                             r'\b(self\.)?(\w+))')
        while True:
            match = pattern.search(text, here)
            if not match: break
            start, end = match.span()
            results.append(escape(text[here:start]))

            (all, link, linkName, scheme, rfc, pep,
             modName, pyModName, selfdot, name) = match.groups()
            if link:
                results.append('<a href="%s">%s</a>' % (link, linkName))
            elif scheme:
                url = escape(all).replace('"', '&quot;')
                results.append('<a href="%s">%s</a>' % (url, url))
            elif rfc:
                url = 'http://www.rfc-editor.org/rfc/rfc%d.txt' % int(rfc)
                results.append('<a href="%s">%s</a>' % (url, escape(all)))
            elif pep:
                url = 'http://www.python.org/dev/peps/pep-%04d/' % int(pep)
                results.append('<a href="%s">%s</a>' % (url, escape(all)))
            elif selfdot:
                results.append('self.<strong>%s</strong>' % name)
            elif modName:
                results.append('<a href="%s">%s</a>' % (modName+
                                                        '.html',modName))
            elif pyModName:
                results.append('<a href="%s/%s">%s</a>' % (pyDocURL,pyModName+
                                                        '.html',pyModName))
            else:
                results.append(self.namelink(name, []))
            #elif text[end:end+1] == '(':
            #    results.append(self.namelink(name, methods, funcs, classes))
            #else:
            #    results.append(self.namelink(name, classes))
            here = end
        results.append(escape(text[here:]))
        return ''.join(results)

    # ---------------------------------------------- type-specific routines

    def docmodule(self, object, name=None, mod=None, *ignored):
        """Produce HTML documentation for a module object."""
        name = object.__name__ # ignore the passed-in name
        try:
            all = object.__all__
        except AttributeError:
            all = None
        parts = name.split('.')
        links = []
        for i in range(len(parts)-1):
            links.append(
                '<a href="%s.html"><font color="#ffffff">%s</font></a>' %
                ('.'.join(parts[:i+1]), parts[i]))
        linkedname = '.'.join(links + parts[-1:])
        head = '<big><big><strong>%s</strong></big></big>' % linkedname
        try:
            path = inspect.getabsfile(object)
            import urllib.parse
            url = urllib.parse.quote(path)
            filelink = self.filelink(url, path)
        except TypeError:
            filelink = '(built-in)'
        info = []
        if hasattr(object, '__version__'):
            version = str(object.__version__)
            if version[:11] == '$' + 'Revision: ' and version[-1:] == '$':
                version = version[11:-1].strip()
            info.append('version %s' % self.escape(version))
        if hasattr(object, '__date__'):
            info.append(self.escape(str(object.__date__)))
        if info:
            head = head + ' (%s)' % ', '.join(info)
        docloc = self.getdocloc(object)
        if docloc is not None:
            docloc = '<br><a href="%(docloc)s">Module Reference</a>' % locals()
        else:
            docloc = ''
        result = self.heading(
            head, '#ffffff', '#7799ee',
            '<a href=".">index</a><br>' + filelink + docloc)

        modules0 = inspect.getmembers(object, inspect.ismodule)
        modules=[]
        for module in modules0:
            if module in includedModules: modules.append(module)
            pass

        classes, cdict = [], {}
        for key, value in inspect.getmembers(object, inspect.isclass):
            # if __all__ exists, believe it.  Otherwise use old heuristic.
            if (all is not None or
                (inspect.getmodule(value) or object) is object):
                if pydoc.visiblename(key, all, object):
                    classes.append((key, value))
                    cdict[key] = cdict[value] = '#' + key
        for key, value in classes:
            for base in value.__bases__:
                key, modname = base.__name__, base.__module__
                module = sys.modules.get(modname)
                if modname != name and module and hasattr(module, key):
                    if getattr(module, key) is base:
                        if not key in cdict:
                            if modname in pythonModules:
                                modname = pyDocURL + modname
                            cdict[key] = cdict[base] = modname + '.html#' + key
        funcs, fdict = [], {}
        for key, value in inspect.getmembers(object, inspect.isroutine):
            # if __all__ exists, believe it.  Otherwise use old heuristic.
            if (all is not None or
                inspect.isbuiltin(value) or inspect.getmodule(value) is object):
                if pydoc.visiblename(key, all, object):
                    funcs.append((key, value))
                    fdict[key] = '#-' + key
                    if inspect.isfunction(value): fdict[value] = fdict[key]
        data = []
        for key, value in inspect.getmembers(object, pydoc.isdata):
            if pydoc.visiblename(key, all, object):
                data.append((key, value))

        doc = self.markup(getdoc(object), self.preformat, fdict, cdict)
        doc = doc and '<tt>%s</tt>' % doc
        result = result + '<p>%s</p>\n' % doc

        if hasattr(object, '__path__'):
            modpkgs = []
            for importer, modname, ispkg in pkgutil.iter_modules(object.__path__):
                modpkgs.append((modname, name, ispkg, 0))
            modpkgs.sort()
            contents = self.multicolumn(modpkgs, self.modpkglink)
            result = result + self.bigsection(
                'Package Contents', '#ffffff', '#aa55cc', contents)
        elif modules:
            pass
#            contents = self.multicolumn(
#                modules, lambda t: self.modulelink(t[1]))
#            result = result + self.bigsection(
#                'Modules', '#ffffff', '#aa55cc', contents)

        if classes:
            classlist = [value for (key, value) in classes]
            contents = [
                self.formattree(inspect.getclasstree(classlist, 1), name)]
            for key, value in classes:
                contents.append(self.document(value, key, name, fdict, cdict))
            result = result + self.bigsection(
                'Classes', '#ffffff', '#ee77aa', ' '.join(contents))
        if funcs:
            contents = []
            for key, value in funcs:
                contents.append(self.document(value, key, name, fdict, cdict))
            result = result + self.bigsection(
                'Functions', '#ffffff', '#eeaa77', ' '.join(contents))
        if data:
            contents = []
            for key, value in data:
                contents.append(self.document(value, key))
            result = result + self.bigsection(
                'Data', '#ffffff', '#55aa55', '<br>\n'.join(contents))
        if hasattr(object, '__author__'):
            contents = self.markup(str(object.__author__), self.preformat)
            result = result + self.bigsection(
                'Author', '#ffffff', '#7799ee', contents)
        if hasattr(object, '__credits__'):
            contents = self.markup(str(object.__credits__), self.preformat)
            result = result + self.bigsection(
                'Credits', '#ffffff', '#7799ee', contents)

        return result

    def docdata(self, object, name=None, mod=None, cl=None):
        """Produce html documentation for a data descriptor."""
        return self._docdescriptor(name, object, mod)

    def index(self, dir, shadowed=None):
        """Generate an HTML index for a directory of modules."""
        modpkgs = []
        if shadowed is None: shadowed = {}
        for importer, name, ispkg in pkgutil.iter_modules([dir]):
            ##FIX: should this list only files ending in .py?
            if any((0xD800 <= ord(ch) <= 0xDFFF) for ch in name):
                # ignore a module if its name contains a surrogate character
                continue
            modpkgs.append((name, '', ispkg, name in shadowed))
            shadowed[name] = 1

        modpkgs.sort()
        contents = self.multicolumn(modpkgs, self.modpkglink)
        return self.bigsection(dir, '#ffffff', '#ee77aa', contents)

class TextDoc(pydoc.TextDoc):
    """Formatter class for text documentation."""

    # ------------------------------------------- text formatting utilities

    _repr_instance = pydoc.TextRepr()
    repr = _repr_instance.repr

    def section(self, title, contents):
        """Format a section with a given heading."""
        return pydoc.TextDoc.section(self,title, self.markup(contents))

    def markup(self,text):
        """remove markup tags."""
        results = []
        here = 0
        pattern = re.compile(r'(<l\s+(\S+)\s+([^>]+)>|'
                             r'\<m\s+(\S+\w)\>|'
                             r'\b(self\.)?(\w+))')
        while True:
            match = pattern.search(text, here)
            if not match: break
            start, end = match.span()
            results.append(text[here:start])

            all, link,linkName, modName, selfdot, name = match.groups()
            if link:
                results.append('%s' % linkName)
            elif modName:
                results.append('%s' % modName)
            else:
                results.append(name)
                pass
            here = end
            pass
        results.append(text[here:])
        return ''.join(results)

    def docmodule(self, object, name=None, mod=None):
        """Produce text documentation for a given module object."""
        name = object.__name__ # ignore the passed-in name
        synop, desc = pydoc.splitdoc(getdoc(object))
        result = self.section('NAME', name + (synop and ' - ' + synop))
        all = getattr(object, '__all__', None)
        docloc = self.getdocloc(object)
        if docloc is not None:
            result = result + self.section('MODULE REFERENCE', docloc + """

The following documentation is automatically generated from the Python
source files.  It may be incomplete, incorrect or include features that
are considered implementation detail and may vary between Python
implementations.  When in doubt, consult the module reference at the
location listed above.
""")

        if desc:
            result = result + self.section('DESCRIPTION', desc)

        classes = []
        for key, value in inspect.getmembers(object, inspect.isclass):
            # if __all__ exists, believe it.  Otherwise use old heuristic.
            if (all is not None
                or (inspect.getmodule(value) or object) is object):
                if pydoc.visiblename(key, all, object):
                    classes.append((key, value))
        funcs = []
        for key, value in inspect.getmembers(object, inspect.isroutine):
            # if __all__ exists, believe it.  Otherwise use old heuristic.
            if (all is not None or
                inspect.isbuiltin(value) or inspect.getmodule(value) is object):
                if pydoc.visiblename(key, all, object):
                    funcs.append((key, value))
        data = []
        for key, value in inspect.getmembers(object, pydoc.isdata):
            if pydoc.visiblename(key, all, object):
                data.append((key, value))

        modpkgs = []
        modpkgs_names = set()
        if hasattr(object, '__path__'):
            for importer, modname, ispkg in pkgutil.iter_modules(object.__path__):
                modpkgs_names.add(modname)
                if ispkg:
                    modpkgs.append(modname + ' (package)')
                else:
                    modpkgs.append(modname)

            modpkgs.sort()
            result = result + self.section(
                'PACKAGE CONTENTS', '\n'.join(modpkgs))

        # Detect submodules as sometimes created by C extensions
        submodules = []
        for key, value in inspect.getmembers(object, inspect.ismodule):
            if value.__name__.startswith(name + '.') and key not in modpkgs_names:
                submodules.append(key)
        if submodules:
            submodules.sort()
            result = result + self.section(
                'SUBMODULES', '\n'.join(submodules))

        if classes:
            classlist = [value for key, value in classes]
            contents = [self.formattree(
                inspect.getclasstree(classlist, 1), name)]
            for key, value in classes:
                contents.append(self.document(value, key, name))
            result = result + self.section('CLASSES', '\n'.join(contents))

        if funcs:
            contents = []
            for key, value in funcs:
                contents.append(self.document(value, key, name))
            result = result + self.section('FUNCTIONS', '\n'.join(contents))

        if data:
            contents = []
            for key, value in data:
                contents.append(self.docother(value, key, name, maxlen=70))
            result = result + self.section('DATA', '\n'.join(contents))

        if hasattr(object, '__version__'):
            version = str(object.__version__)
            if version[:11] == '$' + 'Revision: ' and version[-1:] == '$':
                version = version[11:-1].strip()
            result = result + self.section('VERSION', version)
        if hasattr(object, '__date__'):
            result = result + self.section('DATE', str(object.__date__))
        if hasattr(object, '__author__'):
            result = result + self.section('AUTHOR', str(object.__author__))
        if hasattr(object, '__credits__'):
            result = result + self.section('CREDITS', str(object.__credits__))
        try:
            file = inspect.getabsfile(object)
        except TypeError:
            file = '(built-in)'
        result = result + self.section('FILE', file)
        return result




# --------------------------------------- interactive interpreter interface

text = TextDoc()
#plaintext = _PlainTextDoc() FIX: ??
html = HTMLDoc()

#def doc(thing, title='Python Library Documentation: %s'):
#    """Display text documentation, given an object or a path to an object."""
#    try:
#        object, name = resolve(thing)
#        desc = describe(object)
#        module = inspect.getmodule(object)
#        if name and '.' in name:
#            desc += ' in ' + name[:name.rfind('.')]
#        elif module and module is not object:
#            desc += ' in module ' + module.__name__
#        pager(title % desc + '\n\n' + text.document(object, name))
#    except (ImportError, ErrorDuringImport) as value:
#        print(value)
#        pass
#    return

def render_doc(thing, title='Python Library Documentation: %s', forceload=0,
        renderer=None):
    """Render text documentation, given an object or a path to an object."""
    if renderer is None:
        renderer = text
    object, name = pydoc.resolve(thing, forceload)
    desc = pydoc.describe(object)
    module = inspect.getmodule(object)
    if name and '.' in name:
        desc += ' in ' + name[:name.rfind('.')]
    elif module and module is not object:
        desc += ' in module ' + module.__name__
        pass
    
    if not (inspect.ismodule(object) or
              inspect.isclass(object) or
              inspect.isroutine(object) or
              inspect.isgetsetdescriptor(object) or
              inspect.ismemberdescriptor(object) or
              isinstance(object, property)):
        # If the passed object is a piece of data or an instance,
        # document its available methods instead of its value.
        object = type(object)
        desc += ' object'
    return title % desc + '\n\n' + renderer.document(object, name)

def doc(thing, title='Python Library Documentation: %s', forceload=0,
        output=None):
    """Display text documentation, given an object or a path to an object."""
    try:
        if output is None:
            pydoc.pager(render_doc(thing, title, forceload))
        else:
            output.write(render_doc(thing, title, forceload, plaintext))
    except (ImportError, pydoc.ErrorDuringImport) as value:
        print(value)


modulesProcessed = []


def writedoc(thing, forceload=0):
    """Write HTML documentation to a file in the current directory."""
    try:
        object, name = pydoc.resolve(thing, forceload)
        page = html.page(pydoc.describe(object), html.document(object, name))
        with open(name + '.html', 'w', encoding='utf-8') as file:
            file.write(page)
        print('wrote', name + '.html')
        modulesProcessed.append(name)
    except (ImportError, pydoc.ErrorDuringImport) as value:
        print(value)

def writeIndex(moduleList):
    file = open('index.html','w')
    pkgs = [(x,'',0,0) for x in moduleList]
    pkgs.sort()
    contents = html.multicolumn(pkgs,html.modpkglink)
    contents += r'''
    <hr><p>Search this website for a keyword or phrase

'''

    contents += r'''
<form action="https://www.google.com/search" class="searchform" method="get" name="searchform" target="_blank">
<input name="sitesearch" type="hidden" value="nmr.cit.nih.gov">
<input autocomplete="on" class="form-controls search" name="q" placeholder="Search in nmr.cit.nih.gov" required="required"  type="text">
<button class="button" type="submit">Search</button>
</form>
'''

    contents += '</p>'
    file.write(html.page("Index of Xplor-NIH Modules",
                         html.bigsection("Index of Xplor-NIH Python modules",
                                         '#ffffff', '#ee77aa',contents) ))
    return

class Helper(pydoc.Helper):
    def help(self, request):
        if type(request) is type(''):
            if request == 'help': self.intro()
            elif request == 'keywords': self.listkeywords()
            elif request == 'topics': self.listtopics()
            elif request == 'modules': self.listmodules()
            elif request[:8] == 'modules ':
                self.listmodules(split(request)[1])
            elif request in self.keywords: self.showtopic(request)
            elif request in self.topics: self.showtopic(request)
            elif request: doc(request, 'Help on %s:')
        elif isinstance(request, Helper): self()
        else: doc(request, 'Help on %s:')
        self.output.write('\n')


help = Helper()
 
class ModuleScanner(pydoc.ModuleScanner):
    """An interruptible scanner that searches module synopses."""
    def run(self, callback, key=None, completer=None):
        if key: key = lower(key)
        self.quit = False
        seen = {}

        for modname in sys.builtin_module_names:
            if modname != '__main__':
                seen[modname] = 1
                if key is None:
                    callback(None, modname, '')
                else:
                    if self.scanContents:
                        desc = getdoc(__import__(modname)) or ''
                    else:
                        desc = split(getdoc(__import__(modname)) or '', '\n')[0]
                        pass
                    if find(lower(modname + ' - ' + desc), key) >= 0:
                        callback(None, modname, desc)

        while not self.quit:
            node = self.next()
            if not node: break
            path, package = node
            modname = inspect.getmodulename(path)
            if os.path.isfile(path) and modname:
                xplorModule=0
                for xplorPath in xplorPaths.values():
                    if os.path.dirname(path)==os.path.dirname(xplorPath):
                        xplorModule=1
                    pass
                modname = package + (package and '.') + modname
                if not modname in seen:
                    seen[modname] = 1 # if we see spam.py, skip spam.pyc
                    if key is None:
                        callback(path, modname, '')
                    else:
                        if xplorModule:
                            desc = xplorSynopsis(modname,{},
                                                 self.scanContents) or ''
                        else:
                            desc = synopsis(path,{},self.scanContents) or ''
                        if find(lower(modname + ' - ' + desc), key) >= 0:
                            callback(path, modname, desc)
        if completer: completer()


def apropos(key,scanContents=False):
    """Print all the one-line module summaries that contain a substring.
    If scanContents argument is True, scan full module documentation for
    the substring.
    """
    def callback(path, modname, desc):
        if modname[-9:] == '.__init__':
            modname = modname[:-9] + ' (package)'
        print(modname, desc and '- ' + desc)
    try: import warnings
    except ImportError: pass
    else: warnings.filterwarnings('ignore') # ignore problems during import
    ModuleScanner().run(callback, key, scanContents)

# --------------------------------------------------- web browser interface

def _url_handler(url, content_type="text/html"):
    """The pydoc url handler for use with the pydoc server.

    If the content_type is 'text/css', the _pydoc.css style
    sheet is read and returned if it exits.

    If the content_type is 'text/html', then the result of
    get_html_page(url) is returned.
    """
    class _HTMLDoc(HTMLDoc):

        def page(self, title, contents):
            """Format an HTML page."""
            css_path = "pydoc_data/bin.Linux_x86_64/lib/python3.7/_pydoc.css"
            import os
            css_path = "bin."+os.environ['ARCH']
            import sys
            css_path += "/lib/python{}.{}".format( *sys.version_info[:2] )
            css_path += "/pydoc_data/_pydoc.css"
            css_link = (
                '<link rel="stylesheet" type="text/css" href="%s">' %
                css_path)
            return '''\
<!DOCTYPE html PUBLIC "-//W3C//DTD HTML 4.0 Transitional//EN">
<html><head><title>Pydoc: %s</title>
<meta http-equiv="Content-Type" content="text/html; charset=utf-8">
%s</head><body bgcolor="#f0f0f8">%s<div style="clear:both;padding-top:.5em;">%s</div>
</body></html>''' % (title, css_link, html_navbar(), contents)

        def filelink(self, url, path):
            return '<a href="getfile?key=%s">%s</a>' % (url, path)


    html = _HTMLDoc()

    def html_navbar():
        import platform
        version = html.escape("%s [%s, %s]" % (platform.python_version(),
                                               platform.python_build()[0],
                                               platform.python_compiler()))
        return """
            <div style='float:left'>
                Python %s<br>%s
            </div>
            <div style='float:right'>
                <div style='text-align:center'>
                  <a href="index.html">Module Index</a>
                  : <a href="topics.html">Topics</a>
                  : <a href="keywords.html">Keywords</a>
                </div>
                <div>
                    <form action="get" style='display:inline;'>
                      <input type=text name=key size=15>
                      <input type=submit value="Get">
                    </form>&nbsp;
                    <form action="search" style='display:inline;'>
                      <input type=text name=key size=15>
                      <input type=submit value="Search">
                    </form>
                </div>
            </div>
            """ % (version, html.escape(platform.platform(terse=True)))

    def html_index():
        """Module Index page."""

        def bltinlink(name):
            return '<a href="%s.html">%s</a>' % (name, name)

        heading = html.heading(
            '<big><big><strong>Xplor-NIH: Index of Modules</strong></big></big>',
            '#ffffff', '#7799ee')
        names = [name for name in sys.builtin_module_names
                 if name != '__main__']
        names.append('xplor')
        contents = html.multicolumn(names, bltinlink)
        contents = [heading, '<p>' + html.bigsection(
            'Built-in Modules', '#ffffff', '#ee77aa', contents)]

        seen = {}
#        for dir in sys.path:
#            contents.append(html.index(dir, seen))
        indices = []

        #FIX: add name: 'High-level modules',
        indices.append(html.index(
                                  xplorPaths['high'], seen))
        #FIX: add name: 'Low-level modules',
        indices.append(html.index(
                                  xplorPaths['low'], seen))
        
        
        indices.append('<p>' + html.bigsection(
            'Built-in Modules', '#ffffff', '#ee77aa', contents))
        
        for dir in pathdirs():
            if dir in list(xplorPaths.values()): continue
            #FIX: add name: 'Python modules: '+dir,
            indices.append(html.index(
                                      dir, seen))
            pass
        
        contents = heading + ' '.join(indices) + '''<p align=right>
<font color="#909090" face="helvetica, arial"><strong>
based on pydoc</strong> by Ka-Ping Yee &lt;ping@lfw.org&gt;</font>'''
        
        return 'Index of Modules', ''.join(contents)

    def html_search(key):
        """Search results page."""
        # scan for modules
        search_result = []

        def callback(path, modname, desc):
            if modname[-9:] == '.__init__':
                modname = modname[:-9] + ' (package)'
            search_result.append((modname, desc and '- ' + desc))

        with warnings.catch_warnings():
            warnings.filterwarnings('ignore') # ignore problems during import
            def onerror(modname):
                pass
            ModuleScanner().run(callback, key, onerror=onerror)

        # format page
        def bltinlink(name):
            return '<a href="%s.html">%s</a>' % (name, name)

        results = []
        heading = html.heading(
            '<big><big><strong>Search Results</strong></big></big>',
            '#ffffff', '#7799ee')
        for name, desc in search_result:
            results.append(bltinlink(name) + desc)
        contents = heading + html.bigsection(
            'key = %s' % key, '#ffffff', '#ee77aa', '<br>'.join(results))
        return 'Search Results', contents

    def html_getfile(path):
        """Get and display a source file listing safely."""
        path = urllib.parse.unquote(path)
        with tokenize.open(path) as fp:
            lines = html.escape(fp.read())
        body = '<pre>%s</pre>' % lines
        heading = html.heading(
            '<big><big><strong>File Listing</strong></big></big>',
            '#ffffff', '#7799ee')
        contents = heading + html.bigsection(
            'File: %s' % path, '#ffffff', '#ee77aa', body)
        return 'getfile %s' % path, contents

    def html_topics():
        """Index of topic texts available."""

        def bltinlink(name):
            return '<a href="topic?key=%s">%s</a>' % (name, name)

        heading = html.heading(
            '<big><big><strong>INDEX</strong></big></big>',
            '#ffffff', '#7799ee')
        names = sorted(Helper.topics.keys())

        contents = html.multicolumn(names, bltinlink)
        contents = heading + html.bigsection(
            'Topics', '#ffffff', '#ee77aa', contents)
        return 'Topics', contents

    def html_keywords():
        """Index of keywords."""
        heading = html.heading(
            '<big><big><strong>INDEX</strong></big></big>',
            '#ffffff', '#7799ee')
        names = sorted(Helper.keywords.keys())

        def bltinlink(name):
            return '<a href="topic?key=%s">%s</a>' % (name, name)

        contents = html.multicolumn(names, bltinlink)
        contents = heading + html.bigsection(
            'Keywords', '#ffffff', '#ee77aa', contents)
        return 'Keywords', contents

    def html_topicpage(topic):
        """Topic or keyword help page."""
        buf = io.StringIO()
        htmlhelp = Helper(buf, buf)
        contents, xrefs = htmlhelp._gettopic(topic)
        if topic in htmlhelp.keywords:
            title = 'KEYWORD'
        else:
            title = 'TOPIC'
        heading = html.heading(
            '<big><big><strong>%s</strong></big></big>' % title,
            '#ffffff', '#7799ee')
        contents = '<pre>%s</pre>' % html.markup(contents)
        contents = html.bigsection(topic , '#ffffff','#ee77aa', contents)
        if xrefs:
            xrefs = sorted(xrefs.split())

            def bltinlink(name):
                return '<a href="topic?key=%s">%s</a>' % (name, name)

            xrefs = html.multicolumn(xrefs, bltinlink)
            xrefs = html.section('Related help topics: ',
                                 '#ffffff', '#ee77aa', xrefs)
        return ('%s %s' % (title, topic),
                ''.join((heading, contents, xrefs)))

    def html_getobj(url):
        obj = pydoc.locate(url, forceload=1)
        if obj is None and url != 'None':
            raise ValueError('could not find object')
        title = pydoc.describe(obj)
        content = html.document(obj, url)
        return title, content

    def html_error(url, exc):
        heading = html.heading(
            '<big><big><strong>Error</strong></big></big>',
            '#ffffff', '#7799ee')
        from traceback import format_exception_only
        contents = '<br>'.join(html.escape(line) for line in
                               format_exception_only(type(exc), exc))
        contents = heading + html.bigsection(url, '#ffffff', '#bb0000',
                                             contents)
        return "Error - %s" % url, contents

    def get_html_page(url):
        """Generate an HTML page for url."""
        complete_url = url
        if url.endswith('.html'):
            url = url[:-5]
        try:
            if url in ("", "index"):
                title, content = html_index()
            elif url == "topics":
                title, content = html_topics()
            elif url == "keywords":
                title, content = html_keywords()
            elif '=' in url:
                op, _, url = url.partition('=')
                if op == "search?key":
                    title, content = html_search(url)
                elif op == "getfile?key":
                    title, content = html_getfile(url)
                elif op == "topic?key":
                    # try topics first, then objects.
                    try:
                        title, content = html_topicpage(url)
                    except ValueError:
                        title, content = html_getobj(url)
                elif op == "get?key":
                    # try objects first, then topics.
                    if url in ("", "index"):
                        title, content = html_index()
                    else:
                        try:
                            title, content = html_getobj(url)
                        except ValueError:
                            title, content = html_topicpage(url)
                else:
                    raise ValueError('bad pydoc url')
            else:
                title, content = html_getobj(url)
        except Exception as exc:
            # Catch any errors and display them in an error page.
            title, content = html_error(complete_url, exc)
        return html.page(title, content)

    if url.startswith('/'):
        url = url[1:]
    if content_type == 'text/css':
        import xplorDoc
        path_here = os.path.dirname(os.path.realpath(xplorDoc.__file__))
        css_path = os.path.join(path_here, url)
        with open(css_path) as fp:
            return ''.join(fp.readlines())
    elif content_type == 'text/html':
        return get_html_page(url)
    # Errors outside the url handler are caught by the server.
    raise TypeError('unknown content type %r for url %s' % (content_type, url))


def browse(port=0, *, open_browser=True, hostname='localhost'):
    """Start the enhanced pydoc Web server and open a Web browser.

    Use port '0' to start the server on an arbitrary port.
    Set open_browser to False to suppress opening a browser.
    """
    import webbrowser
    serverthread = pydoc._start_server(_url_handler, hostname, port)
    if serverthread.error:
        print(serverthread.error)
        return
    if serverthread.serving:
        server_help_msg = 'Server commands: [b]rowser, [q]uit'
        if open_browser:
            webbrowser.open(serverthread.url)
        try:
            print('Server ready at', serverthread.url)
            print(server_help_msg)
            while serverthread.serving:
                cmd = input('server> ')
                cmd = cmd.lower()
                if cmd == 'q':
                    break
                elif cmd == 'b':
                    webbrowser.open(serverthread.url)
                else:
                    print(server_help_msg)
        except (KeyboardInterrupt, EOFError):
            print()
        finally:
            if serverthread.serving:
                serverthread.stop()
                print('Server stopped')

def serve(port, callback=None, completer=None):
    import http.server, mimetools, select

    # Patch up mimetools.Message so it doesn't break if rfc822 is reloaded.
    class Message(mimetools.Message):
        def __init__(self, fp, seekable=1):
            Message = self.__class__
            Message.__bases__[0].__bases__[0].__init__(self, fp, seekable)
            self.encodingheader = self.getheader('content-transfer-encoding')
            self.typeheader = self.getheader('content-type')
            self.parsetype()
            self.parseplist()

    class DocHandler(http.server.BaseHTTPRequestHandler):
        def send_document(self, title, contents):
            try:
                self.send_response(200)
                self.send_header('Content-Type', 'text/html')
                self.end_headers()
                self.wfile.write(html.page(title, contents))
            except IOError: pass

        def do_GET(self):
            path = self.path
            if path[-5:] == '.html': path = path[:-5]
            if path[:1] == '/': path = path[1:]
            if path and path != '.':
                try:
                    obj = locate(path)
                except pydoc.ErrorDuringImport as value:
                    self.send_document(path, html.escape(str(value)))
                    return
                if obj:
                    self.send_document(pydoc.describe(obj), html.document(obj, path))
                else:
                    self.send_document(path,
'no Python documentation found for %s' % repr(path))
            else:
                heading = html.heading(
'<big><big><strong>Xplor-NIH: Index of Python Modules</strong></big></big>',
'#ffffff', '#7799ee')
                def bltinlink(name):
                    return '<a href="%s.html">%s</a>' % (name, name)
                names = list([x for x in sys.builtin_module_names if x != '__main__'])
                names.append('xplor')
                contents = html.multicolumn(names, bltinlink)
                indices = []
                seen = {}
                
                indices.append(html.index('High-level modules',
                                          xplorPaths['high'], seen))
                indices.append(html.index('Low-level modules',
                                          xplorPaths['low'], seen))

                
                indices.append('<p>' + html.bigsection(
                    'Built-in Modules', '#ffffff', '#ee77aa', contents))

                for dir in pathdirs():
                    if dir in list(xplorPaths.values()): continue
                    indices.append(html.index('Python modules: '+dir,
                                              dir, seen))

                contents = heading + join(indices) + '''<p align=right>
<font color="#909090" face="helvetica, arial"><strong>
based on pydoc</strong> by Ka-Ping Yee &lt;ping@lfw.org&gt;</font>'''
                self.send_document('Index of Modules', contents)

        def log_message(self, *args): pass

    class DocServer(http.server.HTTPServer):
        def __init__(self, port, callback):
            host = (sys.platform == 'mac') and '127.0.0.1' or 'localhost'
            self.address = ('', port)
            self.url = 'http://%s:%d/' % (host, port)
            self.callback = callback
            self.base.__init__(self, self.address, self.handler)

        def serve_until_quit(self):
            import select
            self.quit = False
            while not self.quit:
                rd, wr, ex = select.select([self.socket.fileno()], [], [], 1)
                if rd: self.handle_request()

        def server_activate(self):
            self.base.server_activate(self)
            if self.callback: self.callback(self)

    DocServer.base = http.server.HTTPServer
    DocServer.handler = DocHandler
    DocHandler.MessageClass = Message
    try:
        try:
            DocServer(port, callback).serve_until_quit()
        except (KeyboardInterrupt, select.error):
            pass
    finally:
        if completer: completer()

# -------------------------------------------------- command-line interface

def includeDirModules(dir):
    for path in os.listdir(dir):
        if ispackage(path):
            includeDirModules(path)
        elif os.path.isfile(path):
            mod = importfile(path)
            object, name = pydoc.resolve(mod)
            includedModules.append(name)
            pass
        pass
    return

    
def cli():
    """Command-line interface (looks at sys.argv to decide what to do)."""
    import getopt
    class BadUsage(Exception): pass

    pydoc._adjust_cli_sys_path()

    import os
    try:
        opts, args = getopt.getopt(sys.argv[1:], 'bk:n:p:w')
        writing = False
        start_server = False
        open_browser = False
        port = 0
        hostname = 'localhost'
        for opt, val in opts:
            if opt == '-b':
                start_server = True
                open_browser = True
            if opt == '-k':
                apropos(val)
                return
            if opt == '-p':
                start_server = True
                port = val
            if opt == '-w':
                writing = True
            if opt == '-n':
                start_server = True
                hostname = val

        if start_server:
            browse(port, hostname=hostname, open_browser=open_browser)
            return

        if not args: raise BadUsage
        # first figure out which modules we'll be documenting-
        # this so that we don't make links to non-existant modules.
        if writing:
            for arg in args:
                if pydoc.ispath(arg) and not os.path.exists(arg):
                    break
                try:
                    if pydoc.ispath(arg) and os.path.isfile(arg):
                        arg = importfile(arg)
                        object, name = pydoc.resolve(arg)
                        includedModules.append(name)
                    elif pydoc.ispath(arg) and os.path.isdir(arg):
                        includeDirModules(arg)
                        pass
                    pass
                except pydoc.ErrorDuringImport as value:
                    pass
                pass
            pass
        for arg in args:
            if pydoc.ispath(arg) and not os.path.exists(arg):
                print('file %r does not exist' % arg)
                break
            try:
                if pydoc.ispath(arg) and os.path.isfile(arg):
                    arg = importfile(arg)
                if writing:
                    if pydoc.ispath(arg) and os.path.isdir(arg):
                        writedocs(arg)
                    else:
                        writedoc(arg)
                else:
                    help.help(arg)
            except pydoc.ErrorDuringImport as value:
                print(value)
        if writing and modulesProcessed:
            writeIndex(modulesProcessed)
            pass

    except (getopt.error, BadUsage):
        cmd = os.path.splitext(os.path.basename(sys.argv[0]))[0]
        print("""pydoc - the Python documentation tool

{cmd} <name> ...
    Show text documentation on something.  <name> may be the name of a
    Python keyword, topic, function, module, or package, or a dotted
    reference to a class or function within a module or module in a
    package.  If <name> contains a '{sep}', it is used as the path to a
    Python source file to document. If name is 'keywords', 'topics',
    or 'modules', a listing of these things is displayed.

{cmd} -k <keyword>
    Search for a keyword in the synopsis lines of all available modules.

{cmd} -n <hostname>
    Start an HTTP server with the given hostname (default: localhost).

{cmd} -p <port>
    Start an HTTP server on the given port on the local machine.  Port
    number 0 can be used to get an arbitrary unused port.

{cmd} -b
    Start an HTTP server on an arbitrary unused port and open a Web browser
    to interactively browse documentation.  This option can be used in
    combination with -n and/or -p.

{cmd} -w <name> ...
    Write out the HTML documentation for a module to a file in the current
    directory.  If <name> contains a '{sep}', it is treated as a filename; if
    it names a directory, documentation is written for all the contents.
""".format(cmd=cmd, sep=os.sep))

#
# fixup code in HTMLutil.py
#
string_re = re.compile(r"(('[^'\n]+')|" + r'("[^"\n]+"))')
markup_re = re.compile(r"[^\n]*((</FONT)|(</STRONG))")
def find_string_literal( s, begin=0 ):
    """find single-line strings which aren't in HTML tags, and aren't
    already marked up

    This won't work inside multiline regions, such as docstrings - strings
    declared within the docstrings will be marked up.
    """
    while 1:
        match = string_re.search(s, begin)
        if not match: break
        if markup_re.match(s,match.start(0)):
            begin = match.end(0)
        else:
            return (match.start(0),match.end(0))
        pass
    return (None, None)


def source2HTML(filename):
    """Colorize/markup Python source code.
    
    Pass filename.  writes output to filename.html
    
    Strings are marked as green, while comments are made red.
    doc strings are made Blue. Functions are made Bold and
    classes Copper.

    
    """
    import HTMLgen, HTMLcolors, HTMLutil
    
    Red = HTMLgen.Font(color=HTMLcolors.RED3)
    Blue = HTMLgen.Font(color=HTMLcolors.BLUE3)
    Green = HTMLgen.Font(color=HTMLcolors.GREEN6)
    Copper = HTMLgen.Font(color=HTMLcolors.COPPER)
    Bold = HTMLgen.Strong()
    def addModuleLinks( text):
        '''
        in the following, substitute NAME with
        <a href="module_dir/NAME">NAME</a>
        from NAME import
        import NAME
        if NAME is in the list XPLORNIH_MODULES
        '''
        from os import environ as env
        from sys import version
        versionStr=version.split()[0]
        try:
            xplorNIH_modules=env['XPLORNIH_MODULES'].split()
        except KeyError:
            xplorNIH_modules=[]
            pass
        try:
            docdir=env['XPLORNIH_DOCDIR']
        except KeyError:
            docdir=[]
            pass
        results = []
        here = 0
        pattern = re.compile(r'^(\s*)(from\s+(\S+)\s+import|'
                             r'import\s+(\S+))',re.MULTILINE)
        while True:
            match = pattern.search(text, here)
            if not match: break
            start, end = match.span()
            results.append(text[here:start])

            leadingSpace, clause, mod1, mod2 = match.groups()
            mod = mod1 or mod2

            if mod in xplorNIH_modules:
                url= docdir + "/" + mod + ".html"
            else:
                url='http://www.python.org/doc/%s/lib/module-%s.html' % \
                     (versionStr,mod)
                pass

            replacement = '<a href="%s">%s</a>' % (url,mod)
            
            if mod1:
                results.append('%sfrom %s import' % (leadingSpace,replacement))
            else:
                results.append('%simport %s' % (leadingSpace,replacement))
            here = end
            pass
        results.append(text[here:])
        return join(results, '')

    def markup(source):
        # it's important to mark the strings first as the HTML tags
        # will be used to detect if key characters (like #) occur in
        # string literals.
        source = HTMLutil.global_substitute(HTMLutil.find_comment,
                                            Red, source)
        source = HTMLutil.global_substitute(find_string_literal,
                                            Green, source)
        source = HTMLutil.global_substitute(HTMLutil.find_docstring,
                                            Green, source)
        source = HTMLutil.global_substitute(HTMLutil.find_function,
                                            Copper, source)
        source = HTMLutil.global_substitute(HTMLutil.find_class,
                                            Bold, source)
        source = addModuleLinks(source)
        return HTMLgen.Pre(source)
    source = open(filename).read()
    import os
    d = HTMLgen.SimpleDocument(bgcolor=HTMLcolors.GREY3,
                               title=os.path.split(filename)[1])
    d.append(markup(source))
    d.write(filename+'.html')


        

if __name__ == '__main__': 
    cli()


