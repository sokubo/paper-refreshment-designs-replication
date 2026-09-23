#!/usr/bin/env python3
"""compare_outputs.py — file-by-file, cell-by-cell comparison of regenerated outputs against shipped ones.
Standard library only.  Used by the release check of the replication archive: every output regenerated
from a clean run is compared with the shipped one.  (The same comparator is used, unchanged in its rules, by
the author's companion replication archive.)

    python3 compare_outputs.py NEW_DIR OLD_DIR [--tol 1e-8] [--verdict FILE]
    python3 compare_outputs.py --selftest DIR            # DIR holds the shipped outputs; corrupted copies are made in a temp dir

Rules (there is no fallback: a file either satisfies them or the comparison FAILS):
  * every *.csv and *.txt present in either directory is compared; a file present in only one directory is a failure;
  * CSV: identical header (column names, in order), identical numbers of rows and columns, and cell by cell in file
    order: a missing cell ("" / NA / NaN) must be missing in both; a cell that parses as a number in both must agree
    to within --tol in absolute value (integers included, so a count of 100 against 101 fails); any other cell must be
    byte-identical (labels, scenario names, verdict strings);
  * TXT (console transcript): identical number of non-blank lines after stripping trailing whitespace and normalising
    line endings; each line is split on whitespace into tokens; the two lines must have the same number of tokens and
    each token is compared by the CSV cell rule (NA / NaN missing in both; numbers within --tol; anything else identical);
  * everything else in the two directories is ignored (the release check copies only outputs into them).
Beyond whitespace (line endings, trailing and inter-token spacing, blank lines) and numeric spellings that agree to
within --tol, the only print-format difference tolerated is therefore an empty cell / NA / NaN for a missing entry,
which R's printers emit inconsistently across versions for a mean over an empty set.  Exit status 0 when every file passes, 1 otherwise.
--verdict writes two lines: the largest absolute numeric difference over all compared cells, and the number of files
that failed (0 means the comparison passed).
"""
import csv, io, math, os, shutil, sys, tempfile

MISSING = {"", "NA", "NaN", "nan", "<NA>"}

def parse_num(tok):
    try:
        v = float(tok)
    except ValueError:
        return None
    if math.isnan(v):
        return None          # "nan" is treated as missing, not as a number
    return v

def cell_diff(a, b, tol):
    """returns (ok, absdiff or None, reason)"""
    ma, mb = a.strip() in MISSING, b.strip() in MISSING
    if ma or mb:
        return (ma and mb, None, None if (ma and mb) else "missing in one file only (%r vs %r)" % (a, b))
    na, nb = parse_num(a), parse_num(b)
    if na is not None and nb is not None:
        if math.isinf(na) or math.isinf(nb):
            ok = (na == nb)
            return (ok, None, None if ok else "infinite value differs (%r vs %r)" % (a, b))
        d = abs(na - nb)
        return (d <= tol, d, None if d <= tol else "numeric difference %.3e > tol (%r vs %r)" % (d, a, b))
    ok = (a == b)
    return (ok, None, None if ok else "text differs (%r vs %r)" % (a, b))

def read_text(path):
    with io.open(path, "r", encoding="utf-8", newline="") as f:
        return f.read().replace("\r\n", "\n").replace("\r", "\n")

def compare_csv(pnew, pold, tol):
    new = list(csv.reader(io.StringIO(read_text(pnew))))
    old = list(csv.reader(io.StringIO(read_text(pold))))
    # drop physical blank lines only (a record with no field or a single empty field); a record that carries
    # delimiters but only empty fields is a data row and must be counted, so that an appended all-empty
    # record changes the row count and fails the comparison
    blank = lambda r: len(r) == 0 or (len(r) == 1 and not r[0].strip())
    new = [r for r in new if not blank(r)]; old = [r for r in old if not blank(r)]
    if not new or not old:
        return False, 0.0, 0, "empty file (%d vs %d rows)" % (len(new), len(old))
    if new[0] != old[0]:
        return False, 0.0, 0, "header differs: %s vs %s" % (new[0], old[0])
    if len(new) != len(old):
        return False, 0.0, 0, "row count differs: %d vs %d data rows" % (len(new) - 1, len(old) - 1)
    md, ncells = 0.0, 0
    for i, (rn, ro) in enumerate(zip(new, old)):
        if len(rn) != len(ro):
            return False, md, ncells, "row %d has %d vs %d cells" % (i, len(rn), len(ro))
        for j, (a, b) in enumerate(zip(rn, ro)):
            ok, d, why = cell_diff(a, b, tol)
            ncells += 1
            if d is not None: md = max(md, d)
            if not ok:
                return False, md, ncells, "row %d, column %r: %s" % (i, new[0][j] if j < len(new[0]) else j, why)
    return True, md, ncells, None

def compare_txt(pnew, pold, tol):
    ln = [l.rstrip() for l in read_text(pnew).split("\n")]; lo = [l.rstrip() for l in read_text(pold).split("\n")]
    ln = [l for l in ln if l.strip()]; lo = [l for l in lo if l.strip()]
    if len(ln) != len(lo):
        return False, 0.0, 0, "non-blank line count differs: %d vs %d" % (len(ln), len(lo))
    md, ntok = 0.0, 0
    for i, (a, b) in enumerate(zip(ln, lo), 1):
        ta, tb = a.split(), b.split()
        if len(ta) != len(tb):
            return False, md, ntok, "line %d has %d vs %d tokens:\n      new: %s\n      old: %s" % (i, len(ta), len(tb), a.strip(), b.strip())
        for j, (x, y) in enumerate(zip(ta, tb), 1):
            ok, d, why = cell_diff(x, y, tol)
            ntok += 1
            if d is not None: md = max(md, d)
            if not ok:
                return False, md, ntok, "line %d, token %d: %s" % (i, j, why)
    return True, md, ntok, None

def compare_dirs(newdir, olddir, tol, out=sys.stdout):
    files = sorted({f for d in (newdir, olddir) for f in os.listdir(d) if f.lower().endswith((".csv", ".txt"))})
    rows, maxd, nfail = [], 0.0, 0
    for f in files:
        pn, po = os.path.join(newdir, f), os.path.join(olddir, f)
        if not os.path.exists(po):
            rows.append((f, "FAIL: regenerated but not shipped", None, 0)); nfail += 1; continue
        if not os.path.exists(pn):
            rows.append((f, "FAIL: shipped but not regenerated", None, 0)); nfail += 1; continue
        if read_text(pn) == read_text(po):
            ok, md, n, why = (True, 0.0, None, None); status = "identical"
        else:
            ok, md, n, why = (compare_csv if f.lower().endswith(".csv") else compare_txt)(pn, po, tol)
            status = "within tolerance" if ok else "FAIL: " + why
        if md: maxd = max(maxd, md)
        if not ok: nfail += 1
        rows.append((f, status, md, n))
    w = max([len(r[0]) for r in rows] + [4])
    out.write("%-*s  %-18s  %-12s  %s\n" % (w, "file", "max |diff|", "cells/tokens", "status"))
    for f, status, md, n in rows:
        out.write("%-*s  %-18s  %-12s  %s\n" % (w, f, ("%.3e" % md) if md is not None else "-", n if n else "-", status))
    out.write("\nfiles compared: %d; identical: %d; within tolerance: %d; FAILED: %d; largest numeric difference: %.3e; tol: %g\n"
              % (len(rows), sum(r[1] == "identical" for r in rows), sum(r[1] == "within tolerance" for r in rows), nfail, maxd, tol))
    out.write("comparison %s\n" % ("passed" if nfail == 0 else "FAILED"))
    return nfail, maxd

# ------------------------------------------------------------------------------------------------------------------
def selftest(srcdir, tol):
    """corrupted copies of the shipped outputs must FAIL; benign variants must PASS"""
    def edit(d, f, fn):
        p = os.path.join(d, f); s = read_text(p); s2 = fn(s)
        if s2 == s: raise RuntimeError("self-test mutation had no effect on %s" % f)
        io.open(p, "w", encoding="utf-8", newline="").write(s2)
    def csv_rows(s): return [r for r in s.split("\n") if r.strip()]
    def csv_cell(s, row, col, fn):
        rows = [list(csv.reader([r]))[0] for r in csv_rows(s)]; rows[row][col] = fn(rows[row][col])
        buf = io.StringIO(); csv.writer(buf, lineterminator="\n").writerows(rows); return buf.getvalue()
    def first_num_col(s, row):
        r = list(csv.reader([csv_rows(s)[row]]))[0]
        return next(j for j, c in enumerate(r) if parse_num(c) is not None and c.strip() not in MISSING and "." in c)
    def first_na(s):
        for i, r in enumerate(csv_rows(s)):
            cells = list(csv.reader([r]))[0]
            for j, c in enumerate(cells):
                if i > 0 and c.strip() in MISSING: return i, j
        raise RuntimeError("no missing cell")
    csvs = sorted(f for f in os.listdir(srcdir) if f.endswith(".csv")); txts = sorted(f for f in os.listdir(srcdir) if f.endswith(".txt"))
    if not csvs or not txts: sys.exit("self-test needs at least one .csv and one .txt in %s" % srcdir)
    # c0: the first CSV (in sorted order) whose first data row holds a text label (a CSV of numbers only cannot
    # carry the label-change case); t0: the first transcript
    def has_label(f):
        r = csv_rows(read_text(os.path.join(srcdir, f)))
        return len(r) > 3 and any(parse_num(c) is None and c.strip() not in MISSING for c in list(csv.reader([r[1]]))[0])
    def has_decimal(f):
        return any(parse_num(tok) is not None and "." in tok for line in read_text(os.path.join(srcdir, f)).split("\n") for tok in line.split())
    c0 = next((f for f in csvs if has_label(f)), None); t0 = next((f for f in txts if has_decimal(f)), None)
    if c0 is None: sys.exit("self-test needs a CSV whose first data row holds a text label in %s" % srcdir)
    if t0 is None: sys.exit("self-test needs a transcript with a decimal number in %s" % srcdir)
    # a CSV with a text label in the first data row and a numeric cell; a CSV with a missing cell (if any)
    s0 = read_text(os.path.join(srcdir, c0)); label_col = next(j for j, c in enumerate(list(csv.reader([csv_rows(s0)[1]]))[0]) if parse_num(c) is None and c.strip() not in MISSING)
    with_na = next((f for f in csvs if any(c.strip() in MISSING for r in csv_rows(read_text(os.path.join(srcdir, f)))[1:] for c in list(csv.reader([r]))[0])), None)
    int_file = next((f for f in csvs if any(c.strip().isdigit() and c.strip() not in ("0", "1") for r in csv_rows(read_text(os.path.join(srcdir, f)))[1:] for c in list(csv.reader([r]))[0])), None)
    tx = read_text(os.path.join(srcdir, t0))
    txt_num = next(tok for line in tx.split("\n") for tok in line.split() if parse_num(tok) is not None and "." in tok)
    txt_word = next((w for w in ("not identified", "identified", "rank") if w in tx), None)
    cases = []
    def case(name, expect_fail, mut): cases.append((name, expect_fail, mut))
    case("unmodified copy", False, lambda d: None)
    case("CSV text label changed in the first data row (%s)" % c0, True, lambda d: edit(d, c0, lambda s: csv_cell(s, 1, label_col, lambda v: "WRONG_LABEL")))
    if int_file:
        si = read_text(os.path.join(srcdir, int_file)); ri = 1
        ci = next(j for j, c in enumerate(list(csv.reader([csv_rows(si)[ri]]))[0]) if c.strip().isdigit() and c.strip() not in ("0", "1"))
        li = next((j for j, c in enumerate(list(csv.reader([csv_rows(si)[ri]]))[0]) if parse_num(c) is None and c.strip() not in MISSING), None)
        case("CSV integer count +1 (%s)" % int_file, True, lambda d: edit(d, int_file, lambda s: csv_cell(s, ri, ci, lambda v: str(int(v) + 1))))
        if li is not None:
            case("CSV integer count +1 together with a label change in the same row (%s)" % int_file, True,
                 lambda d: edit(d, int_file, lambda s: csv_cell(csv_cell(s, ri, ci, lambda v: str(int(v) + 1)), ri, li, lambda v: "WRONG_LABEL")))
    case("CSV numeric cell perturbed by 1e-6 (%s)" % c0, True, lambda d: edit(d, c0, lambda s: csv_cell(s, 1, first_num_col(s, 1), lambda v: repr(float(v) + 1e-6))))
    case("CSV numeric cell perturbed by 1e-12 (%s)" % c0, False, lambda d: edit(d, c0, lambda s: csv_cell(s, 1, first_num_col(s, 1), lambda v: repr(float(v) + 1e-12))))
    case("CSV numeric cell replaced by NA (%s)" % c0, True, lambda d: edit(d, c0, lambda s: csv_cell(s, 1, first_num_col(s, 1), lambda v: "NA")))
    if with_na:
        case("CSV missing cell replaced by 0 (%s)" % with_na, True, lambda d: edit(d, with_na, lambda s: (lambda ij: csv_cell(s, ij[0], ij[1], lambda v: "0"))(first_na(s))))
        case("CSV missing cell written as NaN instead of empty/NA (%s)" % with_na, False, lambda d: edit(d, with_na, lambda s: (lambda ij: csv_cell(s, ij[0], ij[1], lambda v: "NaN"))(first_na(s))))
    case("CSV all-empty record appended (delimiters only) (%s)" % c0, True, lambda d: edit(d, c0, lambda s: s.rstrip("\n") + "\n" + "," * (len(list(csv.reader([csv_rows(s)[0]]))[0]) - 1) + "\n"))
    case("CSV blank physical lines appended (%s)" % c0, False, lambda d: edit(d, c0, lambda s: s.rstrip("\n") + "\n\n\n"))
    case("CSV last data row duplicated (%s)" % c0, True, lambda d: edit(d, c0, lambda s: s.rstrip("\n") + "\n" + csv_rows(s)[-1] + "\n"))
    case("CSV last data row deleted (%s)" % c0, True, lambda d: edit(d, c0, lambda s: "\n".join(csv_rows(s)[:-1]) + "\n"))
    case("CSV two data rows swapped (%s)" % c0, True, lambda d: edit(d, c0, lambda s: (lambda r: "\n".join([r[0], r[2], r[1]] + r[3:]) + "\n")(csv_rows(s))))
    case("CSV last column dropped (%s)" % c0, True, lambda d: edit(d, c0, lambda s: "\n".join(",".join(r.split(",")[:-1]) for r in csv_rows(s)) + "\n"))
    case("CSV column renamed in the header (%s)" % c0, True, lambda d: edit(d, c0, lambda s: (lambda r: "\n".join([r[0] + "_x"] + r[1:]) + "\n")(csv_rows(s))))
    case("CSV file missing from the regenerated set (%s)" % csvs[-1], True, lambda d: os.remove(os.path.join(d, csvs[-1])))
    case("TXT numeric token changed in the last displayed digit (%s)" % t0, True,
         lambda d: edit(d, t0, lambda s: s.replace(txt_num, txt_num[:-1] + ("0" if txt_num[-1] != "0" else "1"), 1)))
    if txt_word:
        case("TXT verdict/word changed (%r) (%s)" % (txt_word, t0), True, lambda d: edit(d, t0, lambda s: s.replace(txt_word, txt_word.upper() + "_X", 1)))
    if " NA " in tx or "\tNA" in tx or tx.rstrip().endswith("NA"):
        case("TXT NA printed as NaN (%s)" % t0, False, lambda d: edit(d, t0, lambda s: s.replace(" NA ", " NaN ").replace(" NA\n", " NaN\n")))
        case("TXT NA replaced by a number (%s)" % t0, True, lambda d: edit(d, t0, lambda s: s.replace(" NA ", " 0.0 ", 1)))
    case("TXT trailing whitespace and CRLF line endings (%s)" % t0, False, lambda d: edit(d, t0, lambda s: s.replace("\n", "   \r\n")))
    # delete the third NON-BLANK line (deleting a blank line is a benign variant, not a corruption)
    case("TXT one non-blank line deleted (%s)" % t0, True, lambda d: edit(d, t0, lambda s: (lambda L: "\n".join(
        l for i, l in enumerate(L) if i != [k for k, x in enumerate(L) if x.strip()][2]))(s.split("\n"))))
    root = tempfile.mkdtemp("compare_selftest"); allok = True
    print("self-test of compare_outputs.py on corrupted copies of %s (%d cases)\n" % (srcdir, len(cases)))
    for name, expect_fail, mut in cases:
        old = os.path.join(root, "shipped"); new = os.path.join(root, "regen")
        for d in (old, new):
            if os.path.isdir(d): shutil.rmtree(d)
            os.mkdir(d)
            for f in csvs + txts: shutil.copy(os.path.join(srcdir, f), d)
        mut(new)
        buf = io.StringIO(); nfail, _ = compare_dirs(new, old, tol, out=buf)
        failed = nfail > 0; asexp = (failed == expect_fail)
        allok &= asexp
        print("   %-84s -> %-6s %s" % (name, "FAIL" if failed else "pass", "as expected" if asexp else "!! NOT AS EXPECTED (must %s)" % ("FAIL" if expect_fail else "pass")))
        if not asexp: print("      " + buf.getvalue().replace("\n", "\n      "))
    shutil.rmtree(root)
    print("\nself-test %s: %d corruptions %s detected, %d benign variants %s accepted" % (
        "passed" if allok else "FAILED", sum(e for _, e, _ in cases), "all" if allok else "not all", sum(not e for _, e, _ in cases), "all" if allok else "not all"))
    return 0 if allok else 1

def main(argv):
    tol = 1e-8; verdict = None; args = []
    i = 0
    while i < len(argv):
        if argv[i] == "--tol": tol = float(argv[i + 1]); i += 2
        elif argv[i] == "--verdict": verdict = argv[i + 1]; i += 2
        else: args.append(argv[i]); i += 1
    if args and args[0] == "--selftest":
        return selftest(args[1] if len(args) > 1 else ".", tol)
    if len(args) != 2:
        sys.exit(__doc__)
    nfail, maxd = compare_dirs(args[0], args[1], tol)
    if verdict:
        io.open(verdict, "w", encoding="utf-8").write("%.6e\n%d\n" % (maxd, nfail))
    return 0 if nfail == 0 else 1

if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))
