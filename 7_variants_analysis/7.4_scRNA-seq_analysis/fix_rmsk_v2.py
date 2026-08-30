import json

NB = r"C:\Users\mathe\STRs_COVID_ExpAnalysis\7_variants_analysis\7.4_scRNA-seq_analysis\7.4.3.2_trackviewer_STR_viz.ipynb"

with open(NB, encoding="utf-8") as f:
    nb = json.load(f)

src = "".join(nb["cells"][5]["source"])
lines = src.splitlines(keepends=True)

# Find the start and end of the block to replace
# Start at line 73 (0-indexed: 72) - the comment line
# End at line 93 (0-indexed: 92) - the '  } else {' line

new_block = [
    "  # Detect UCSC database dumps mislabelled as BED (the downloader saves\n",
    "  # goldenPath/hg38/database/rmsk.txt.gz as rmsk.hg38.bed.gz).\n",
    "  # Two formats exist:\n",
    "  #  (a) with header:  '#bin swScore ... genoName genoStart genoEnd ... repClass repFamily'\n",
    "  #  (b) headerless:   data starts immediately, 17 cols, col6=chr..., col12=repClass, col13=repFamily\n",
    "  # Both -> remap to a proper 4-column BED (chr, start, end, class/family).\n",
    "  con <- gzfile(path, open = 'rt'); hdr <- readLines(con, n = 1L); close(con)\n",
    "  fields <- strsplit(sub('^#', '', hdr), '\\t', fixed = TRUE)[[1]]\n",
    "  is_ucsc_dump_header <- all(c('genoName', 'genoStart', 'genoEnd', 'repClass',\n",
    "                               'repFamily') %in% fields)\n",
    "  is_ucsc_dump_headerless <- length(fields) >= 13 &&\n",
    "    grepl('^chr', fields[6]) &&\n",
    "    grepl('^[A-Za-z]', fields[12]) && grepl('^[A-Za-z]', fields[13])\n",
    "  is_ucsc_dump <- is_ucsc_dump_header || is_ucsc_dump_headerless\n",
    "  if (is_ucsc_dump) {\n",
    "    if (is_ucsc_dump_header) {\n",
    "      chr <- as.character(dt[[match('genoName', fields)]])\n",
    "      s <- suppressWarnings(as.integer(dt[[match('genoStart', fields)]])) + 1L\n",
    "      e <- suppressWarnings(as.integer(dt[[match('genoEnd', fields)]]))\n",
    "      nm <- paste(dt[[match('repClass', fields)]],\n",
    "                  dt[[match('repFamily', fields)]], sep = '/')\n",
    "    } else {\n",
    "      # headerless: columns are 1-based in R (dt[[1]] = first col)\n",
    "      # col1=bin, col2=swScore, col3=milliDiv, col4=milliDel, col5=milliIns,\n",
    "      # col6=genoName, col7=genoStart, col8=genoEnd, col11=repName,\n",
    "      # col12=repClass, col13=repFamily\n",
    "      chr <- as.character(dt[[6]])\n",
    "      s <- suppressWarnings(as.integer(dt[[7]])) + 1L\n",
    "      e <- suppressWarnings(as.integer(dt[[8]]))\n",
    "      nm <- paste(dt[[12]], dt[[13]], sep = '/')\n",
    "    }\n",
    "    ok <- !is.na(s) & !is.na(e) & e >= s & !is.na(chr) &\n",
    "          grepl('^chr', chr) & !is.na(nm) & nm != 'NA/NA'\n",
    "    gri <- GRanges(chr[ok], IRanges(s[ok], e[ok]))\n",
    "    keep <- sort(unique(S4Vectors::queryHits(\n",
    "      GenomicRanges::findOverlaps(gri, win))))\n",
    "    res <- data.table(V1 = chr[ok][keep], V2 = s[ok][keep],\n",
    "                      V3 = e[ok][keep], V4 = nm[ok][keep])\n",
    "  } else {\n"
]

# Replace lines 72-93 (0-indexed: 72-92 inclusive) with new_block
# lines[72] = "  # goldenPath/hg38/database/rmsk.txt.gz as rmsk.hg38.bed.gz). That format has\n"
# lines[92] = "  } else {"

# Find exact indices
start_idx = None
end_idx = None
for i, line in enumerate(lines):
    if 'goldenPath/hg38/database/rmsk.txt.gz as rmsk.hg38.bed.gz' in line:
        start_idx = i
        break

for i in range(start_idx, len(lines)):
    if lines[i].strip() == '} else {':
        end_idx = i
        break

print(f"start_idx={start_idx}, end_idx={end_idx}")
print(f"Replacing lines {start_idx} to {end_idx}")

# Replace
lines[start_idx:end_idx+1] = [l + '\n' if not l.endswith('\n') else l for l in new_block]

# Reconstruct
src = "".join(lines)

nb["cells"][5]["source"] = src.splitlines(keepends=True)

with open(NB, "w", encoding="utf-8", newline="\n") as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)
    f.write("\n")

print("Done!")
# Verify
with open(NB, encoding="utf-8") as f:
    nb = json.load(f)
src = "".join(nb["cells"][5]["source"])
if "is_ucsc_dump_headerless" in src:
    print("SUCCESS: headerless detection added")
    i = src.index("is_ucsc_dump_headerless")
    print(src[i:i+300])