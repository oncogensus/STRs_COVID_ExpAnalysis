import json

NB = r"C:\Users\mathe\STRs_COVID_ExpAnalysis\7_variants_analysis\7.4_scRNA-seq_analysis\7.4.3.2_trackviewer_STR_viz.ipynb"

with open(NB, encoding="utf-8") as f:
    nb = json.load(f)

src = "".join(nb["cells"][5]["source"])

# Replace the UCSC dump detection logic
old = """  # Detect UCSC database dumps mislabelled as BED (the downloader saves
  # goldenPath/hg38/database/rmsk.txt.gz as rmsk.hg38.bed.gz). That format has
  # a '#bin swScore ... genoName genoStart genoEnd ... repClass repFamily'
  # header instead of BED columns -> remap to a proper 4-column BED.
  con <- gzfile(path, open = 'rt'); hdr <- readLines(con, n = 1L); close(con)
  fields <- strsplit(sub('^#', '', hdr), '\t', fixed = TRUE)[[1]]
  is_ucsc_dump <- all(c('genoName', 'genoStart', 'genoEnd', 'repClass',
                        'repFamily') %in% fields)
  if (is_ucsc_dump) {
    chr <- as.character(dt[[match('genoName', fields)]])
    s <- suppressWarnings(as.integer(dt[[match('genoStart', fields)]])) + 1L
    e <- suppressWarnings(as.integer(dt[[match('genoEnd', fields)]]))
    nm <- paste(dt[[match('repClass', fields)]],
                dt[[match('repFamily', fields)]], sep = '/')
    ok <- !is.na(s) & !is.na(e) & e >= s & !is.na(chr) &
          grepl('^chr', chr) & !is.na(nm) & nm != 'NA/NA'
    gri <- GRanges(chr[ok], IRanges(s[ok], e[ok]))
    keep <- sort(unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gri, win))))
    res <- data.table(V1 = chr[ok][keep], V2 = s[ok][keep],
                      V3 = e[ok][keep], V4 = nm[ok][keep])
  } else {"""

new = """  # Detect UCSC database dumps mislabelled as BED (the downloader saves
  # goldenPath/hg38/database/rmsk.txt.gz as rmsk.hg38.bed.gz).
  # Two formats exist:
  #  (a) with header:  '#bin swScore ... genoName genoStart genoEnd ... repClass repFamily'
  #  (b) headerless:   data starts immediately, 17 cols, col6=chr..., col12=repClass, col13=repFamily
  # Both -> remap to a proper 4-column BED (chr, start, end, class/family).
  con <- gzfile(path, open = 'rt'); hdr <- readLines(con, n = 1L); close(con)
  fields <- strsplit(sub('^#', '', hdr), '\t', fixed = TRUE)[[1]]
  is_ucsc_dump_header <- all(c('genoName', 'genoStart', 'genoEnd', 'repClass',
                               'repFamily') %in% fields)
  is_ucsc_dump_headerless <- length(fields) >= 13 &&
    grepl('^chr', fields[6]) &&
    grepl('^[A-Za-z]', fields[12]) && grepl('^[A-Za-z]', fields[13])
  is_ucsc_dump <- is_ucsc_dump_header || is_ucsc_dump_headerless
  if (is_ucsc_dump) {
    if (is_ucsc_dump_header) {
      chr <- as.character(dt[[match('genoName', fields)]])
      s <- suppressWarnings(as.integer(dt[[match('genoStart', fields)]])) + 1L
      e <- suppressWarnings(as.integer(dt[[match('genoEnd', fields)]]))
      nm <- paste(dt[[match('repClass', fields)]],
                  dt[[match('repFamily', fields)]], sep = '/')
    } else {
      # headerless: columns are 1-based in R (dt[[1]] = first col)
      # col1=bin, col2=swScore, col3=milliDiv, col4=milliDel, col5=milliIns,
      # col6=genoName, col7=genoStart, col8=genoEnd, col11=repName,
      # col12=repClass, col13=repFamily
      chr <- as.character(dt[[6]])
      s <- suppressWarnings(as.integer(dt[[7]])) + 1L
      e <- suppressWarnings(as.integer(dt[[8]]))
      nm <- paste(dt[[12]], dt[[13]], sep = '/')
    }
    ok <- !is.na(s) & !is.na(e) & e >= s & !is.na(chr) &
          grepl('^chr', chr) & !is.na(nm) & nm != 'NA/NA'
    gri <- GRanges(chr[ok], IRanges(s[ok], e[ok]))
    keep <- sort(unique(S4Vectors::queryHits(
      GenomicRanges::findOverlaps(gri, win))))
    res <- data.table(V1 = chr[ok][keep], V2 = s[ok][keep],
                      V3 = e[ok][keep], V4 = nm[ok][keep])
  } else {"""

if old in src:
    src = src.replace(old, new)
    print("Notebook: UCSC dump detection updated")
else:
    print("Could not find old pattern - trying alternative")
    # Try a different approach - find the key line
    if "is_ucsc_dump <- all(c('genoName'" in src:
        print("Found alternative pattern")

nb["cells"][5]["source"] = src.splitlines(keepends=True)

with open(NB, "w", encoding="utf-8", newline="\n") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("Notebook updated successfully")

# Verify
with open(NB, encoding="utf-8") as f:
    nb = json.load(f)
src = "".join(nb["cells"][5]["source"])
if "is_ucsc_dump_headerless" in src:
    print("Verification: headerless detection added")
if "grepl('^chr', fields[6])" in src:
    print("Verification: headerless detection logic present")