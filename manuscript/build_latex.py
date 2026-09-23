#!/usr/bin/env python3
"""main.qmd -> main.pdf + latex/main.tex (PDF build).

House format (matches the reference paper.tex of 2026-09-16): 12pt letterpaper
article, 1in margins, \\setstretch{1.5} body, title \\thanks for acknowledgements
and funding, author \\thanks for affiliation and contact, Keywords line after the
abstract, colorlinks. Everything except the two \\thanks strings comes from the
qmd's `pdf` format block; the \\thanks are patched in here because they are
LaTeX-only and would leak into the HTML if put in the YAML.

Re-run after any edit to main.qmd. Requires quarto + pdflatex.
"""
import subprocess, shutil, os, re, sys

HERE = os.path.dirname(os.path.abspath(__file__))
LTX  = os.path.join(HERE, 'latex')

# --- the two footnotes on the title page -------------------------------------
ARCHIVE_REPO = 'https://github.com/sokubo/paper-refreshment-designs-replication'
ARCHIVE_TAG = 'paper-v0.6'
ARCHIVE_COMMIT = 'COMMIT7'          # filled in after the release check of the published snapshot
TITLE_THANKS = (
    r"\thanks{Code for every simulation and deterministic check in this paper is in the replication "
    r"archive at \url{%s} (fixed version: tag \texttt{%s}, commit \texttt{%s}). The empirical "
    r"illustration uses licensed JLPS microdata, which are not redistributed; the archive contains the "
    r"code that produces the reported aggregates and a synthetic example (see Data and code availability). "
    r"This work was supported by JSPS KAKENHI Grant Number 22K13525.}" % (ARCHIVE_REPO, ARCHIVE_TAG, ARCHIVE_COMMIT)
)
AUTHOR_THANKS = (
    r"\thanks{Department of Sociology, Toyo University, Tokyo, Japan. "
    r"Email: okubo080@toyo.jp. Website: sokubo.github.io.}"
)
DATE = 'September 23, 2026'
KEYWORDS = (r"\noindent\textbf{Keywords:} attrition; panel conditioning; partial identification; "
            r"refreshment samples; rotation panels; survey design")
# -----------------------------------------------------------------------------

_src = open(os.path.join(HERE, 'main.qmd')).read()
_ph = sorted(set(re.findall(r'\[\[[A-Z0-9_]+\]\]', _src)))
if _ph or 'COMMIT7' in ARCHIVE_COMMIT:
    print('WARNING: unfilled placeholders:', _ph + (['ARCHIVE_COMMIT'] if 'COMMIT7' in ARCHIVE_COMMIT else []))
subprocess.run(['quarto', 'render', 'main.qmd', '--to', 'pdf'], cwd=HERE, check=True)

tex_src = os.path.join(HERE, 'main.tex')          # keep-tex output
tex = open(tex_src).read()
if 'CSLReferences' not in tex:
    sys.exit('bibliography not resolved into the tex — check citeproc')

# global theorem numbering (Theorem 1, 2, ...) so that the PDF matches the HTML and the prose cross-references
tex = re.sub(r'(\\newtheorem\{(theorem|corollary|proposition|lemma|definition|refremark|refsolution)\}\{[^}]*\})\[section\]', r'\1', tex)

# title \thanks: attach to the end of the \title{...} argument
m = re.search(r'\\title\{', tex)
if not m:
    sys.exit('no \\title in generated tex')
i, depth = m.end(), 1
while depth:
    if tex[i] == '{': depth += 1
    elif tex[i] == '}': depth -= 1
    i += 1
tex = tex[:i-1] + TITLE_THANKS + tex[i-1:]

# author \thanks
m = re.search(r'\\author\{([^}]*)\}', tex)
if not m:
    sys.exit('no \\author in generated tex')
tex = tex[:m.end()-1] + AUTHOR_THANKS + tex[m.end()-1:]

# visible Keywords line after the abstract (house format)
if r'\textbf{Keywords:}' not in tex:
    tex = tex.replace(r'\end{abstract}', '\\end{abstract}\n\n' + KEYWORDS + '\n', 1)

# date, spelled out as in the house format
tex = re.sub(r'\\date\{[^}]*\}', r'\\date{' + DATE + '}', tex, count=1)

# amsthm's proof environment prints its own tombstone; drop the explicit \square the qmd carries for the HTML
tex, n_sq = re.subn(r'\s*\\\(\\square\\\)\s*\n(\\end\{proof\})', r'\n\1', tex)
print('explicit \\square removed before \\end{proof}:', n_sq)

os.makedirs(LTX, exist_ok=True)
open(os.path.join(LTX, 'main.tex'), 'w').write(tex)
os.remove(tex_src)
shutil.copy(os.path.join(HERE, 'references.bib'), os.path.join(LTX, 'references.bib'))
figs = os.path.join(HERE, 'figures')
if os.path.isdir(figs):
    shutil.copytree(figs, os.path.join(LTX, 'figures'), dirs_exist_ok=True)
# citations are already resolved by citeproc, so a plain two-pass run suffices
for _ in range(3):
    r = subprocess.run(['xelatex', '-interaction=nonstopmode', 'main.tex'],
                       cwd=LTX, capture_output=True, text=True)
open(os.path.join(LTX, 'full.log'), 'w').write(r.stdout[-20000:])

log = open(os.path.join(LTX, 'main.log'), errors='ignore').read()
errs  = [l for l in log.splitlines() if l.startswith('!')]
undef = [l for l in log.splitlines() if 'Reference' in l and 'undefined' in l]
print('LaTeX errors:', len(errs), errs[:5])
print('undefined references:', len(undef), undef[:5])
shutil.copy(os.path.join(LTX, 'main.pdf'), os.path.join(HERE, 'main.pdf'))
print('pages:', subprocess.run(['pdfinfo', os.path.join(HERE, 'main.pdf')],
      capture_output=True, text=True).stdout.split('Pages:')[1].split()[0])
