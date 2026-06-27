import os, glob, re
from pathlib import Path

configfile: "config/config.yaml"

# ── Input paths ─────────────────────────────────────────────────────────────
# Amino-acid marker proteins used to assign contigs → Microcoleus genome IDs.
# Each sequence header must begin with the genome ID (first whitespace-delimited token).
AA_FASTA     = config["aa_fasta"]
ASSEMBLY_DIR = config["assembly_dir"]
OUTDIR       = config["outdir"]
FASTA_EXTS   = config.get("fasta_exts", ["fasta", "fa", "fna"])
ENV          = config.get("conda_env", "microcoleus")

# Wrap every tool call so it runs inside the pre-existing conda environment.
# Use -p (path) when conda_env is an absolute path, -n (name) otherwise.
_ENV_FLAG = "-p" if ENV.startswith("/") else "-n"
DIAMOND  = f"conda run {_ENV_FLAG} {ENV} diamond"
SAMTOOLS = f"conda run {_ENV_FLAG} {ENV} samtools"
BASH     = f"conda run {_ENV_FLAG} {ENV} bash -lc"

def resolve_dir(p):
    try:
        return str(Path(p).resolve())
    except Exception:
        return p

ASSEMBLY_DIR_REAL = resolve_dir(ASSEMBLY_DIR)

def discover_samples():
    """Return sorted list of sample basenames found in the assembly directory."""
    files = []
    for ext in FASTA_EXTS:
        files.extend(glob.glob(os.path.join(ASSEMBLY_DIR_REAL, f"*.{ext}")))
    files = sorted(set(files))
    return sorted({os.path.splitext(os.path.basename(f))[0] for f in files})

SAMPLES = discover_samples()
if not SAMPLES:
    raise ValueError(f"No assembly FASTA files found in {ASSEMBLY_DIR_REAL} with extensions {FASTA_EXTS}")

# ── Output directory layout ──────────────────────────────────────────────────
DBDIR     = os.path.join(OUTDIR, "diamond_db")
TMPDIR    = os.path.join(OUTDIR, "tmp")
PERSAMPLE = os.path.join(OUTDIR, "per_sample")
MATRIXDIR = os.path.join(OUTDIR, "matrices")

GENOME_LIST = os.path.join(TMPDIR, "genome_ids.txt")
DB_PREFIX   = os.path.join(DBDIR, "Dataset2_AA")
DB_DMND     = DB_PREFIX + ".dmnd"
LOGDIR      = os.path.join(OUTDIR, "logs")
SINGLEM_DIR = config.get("singlem_dir", "")
CORRDIR     = os.path.join(OUTDIR, "correlations")

def assembly_fasta_path(sample):
    for ext in FASTA_EXTS:
        p = os.path.join(ASSEMBLY_DIR_REAL, f"{sample}.{ext}")
        if os.path.exists(p):
            return p
    raise ValueError(f"No FASTA found for sample={sample} in {ASSEMBLY_DIR_REAL} ({FASTA_EXTS})")

def bam_path(sample):
    return os.path.join(ASSEMBLY_DIR_REAL, f"{sample}.bam")

# ── Anatoxin detection constants ─────────────────────────────────────────────
# The 10 core genes of the anatoxin-a / dihydroanatoxin-a biosynthetic cluster
# (anaA through anaK, excluding anaH which is variably present and often absent).
# Reference: Méjean et al. 2009, 2016; Kust et al. 2020; Gibis et al. 2026.
ANA_GENES   = ["anaA", "anaB", "anaC", "anaD", "anaE",
               "anaF", "anaG", "anaI", "anaJ", "anaK"]
N_ANA_GENES = len(ANA_GENES)

ANA_DIR     = os.path.join(OUTDIR, "anatoxin")
ANA_REFS    = os.path.join(ANA_DIR, "refs", "ana_refs.faa")
ANA_DB_PFX  = os.path.join(ANA_DIR, "refs", "ana_refs")   # DIAMOND makedb prefix
ANA_DB_DMND = ANA_DB_PFX + ".dmnd"

# ── Reference dataset sources ─────────────────────────────────────────────────
# Alignment datasets from Stanojković et al. (2024) Nat Commun
# "The global speciation continuum of the cyanobacterium Microcoleus"
# DOI: 10.6084/m9.figshare.24710961.v2
FIGSHARE_DATASETS = {
    "Dataset1.AA.alignment.fa":              "https://ndownloader.figshare.com/files/43416798",
    "Dataset2.AA.microcoleus.and.GenBank.fa": "https://ndownloader.figshare.com/files/43416804",
    "Dataset3.AA.microcoleus.without.genbank.fa": "https://ndownloader.figshare.com/files/43416807",
    "Dataset3.all.single.copy.trees.nwk":    "https://ndownloader.figshare.com/files/43416801",
}

# ════════════════════════════════════════════════════════════════════════════════
# Default target: Microcoleus genome abundance across all samples
# ════════════════════════════════════════════════════════════════════════════════

rule all:
    input:
        expand(os.path.join(PERSAMPLE, "{sample}.genome_abundance.tsv"), sample=SAMPLES),
        os.path.join(MATRIXDIR, "genome_Reads_matrix.tsv"),
        os.path.join(MATRIXDIR, "genome_CPM_matrix.tsv"),
        os.path.join(MATRIXDIR, "genome_RPKM_matrix.tsv"),
        os.path.join(MATRIXDIR, "genome_MeanDepth_matrix.tsv"),
        os.path.join(MATRIXDIR, "genome_Breadth_matrix.tsv")

# ════════════════════════════════════════════════════════════════════════════════
# Anatoxin target: run with  snakemake ana_all --cores N
# ════════════════════════════════════════════════════════════════════════════════

rule ana_all:
    # Top-level target for the anatoxin sub-workflow.
    # Runs independently of `rule all`; requires the Microcoleus abundance
    # pipeline to have completed first (it uses contig2genome + coverage outputs).
    input:
        os.path.join(ANA_DIR, "ana_summary.tsv"),
        os.path.join(ANA_DIR, "ana_combined.tsv"),
        expand(os.path.join(ANA_DIR, "per_sample_annotated", "{sample}.ana_abundance.tsv"),
               sample=SAMPLES)

# ════════════════════════════════════════════════════════════════════════════════
# REFERENCE DATASET DOWNLOAD
# ════════════════════════════════════════════════════════════════════════════════

rule datasets_all:
    # Convenience target: download all four reference alignment files.
    # Run with:  snakemake datasets_all --cores 4
    input:
        list(FIGSHARE_DATASETS.keys())

rule download_dataset:
    # Download a single reference alignment file from figshare.
    # Uses wget -c so partial downloads resume rather than restart.
    # Source: Stanojković et al. (2024) Nat Commun, DOI 10.6084/m9.figshare.24710961.v2
    output:
        "{dataset_file}"
    wildcard_constraints:
        # Only match the four known figshare filenames; avoid hijacking other rules.
        dataset_file="|".join(re.escape(f) for f in FIGSHARE_DATASETS)
    log:
        os.path.join(LOGDIR, "download", "{dataset_file}.log")
    params:
        url=lambda wc: FIGSHARE_DATASETS[wc.dataset_file]
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        wget -c -O {output} '{params.url}' > "{log}" 2>&1
        """

# ════════════════════════════════════════════════════════════════════════════════
# MICROCOLEUS ABUNDANCE RULES
# ════════════════════════════════════════════════════════════════════════════════

rule make_dirs:
    # Create the full output directory tree before any other rule runs.
    output:
        ok=os.path.join(OUTDIR, ".dirs.ok")
    log:
        os.path.join(LOGDIR, "make_dirs.log")
    shell:
        r"""
        mkdir -p "{OUTDIR}" "{DBDIR}" "{TMPDIR}" "{PERSAMPLE}" "{MATRIXDIR}" "{LOGDIR}"
        (
        set -euo pipefail
        touch "{output.ok}"
        ) > "{log}" 2>&1
        """

rule genome_list:
    # Extract the unique genome ID (first token of each FASTA header) from the
    # marker protein file.  The resulting list drives the abundance matrix rows.
    input:
        AA_FASTA,
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        GENOME_LIST
    log:
        os.path.join(LOGDIR, "genome_list.log")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        {SAMTOOLS} --version >/dev/null 2>&1 || true
        grep -E '^>' "{input[0]}" | sed 's/^>//' | awk '{{print $1}}' | sort -u > "{output}"

        if [[ ! -s "{output}" ]]; then
          echo "ERROR: genome list is empty. Check AA_FASTA path/content: {input[0]}" >&2
          exit 1
        fi
        ) > "{log}" 2>&1
        """

rule diamond_db:
    # Build the DIAMOND protein database from the Microcoleus marker sequences.
    # This DB is used to assign metagenomic contigs to genome bins.
    input:
        AA_FASTA,
        glist=GENOME_LIST,
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        DB_DMND
    log:
        os.path.join(LOGDIR, "diamond_db.log")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")"
        (
        set -euo pipefail
        {DIAMOND} makedb --in "{input[0]}" -d "{DB_PREFIX}"
        ) > "{log}" 2>&1
        """

rule index_bam:
    # Index each sample's BAM file so samtools coverage can do random access.
    # Handles both .bai and .csi index formats.
    input:
        bam=lambda wc: bam_path(wc.sample),
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        ok=os.path.join(TMPDIR, "bam_indexed", "{sample}.ok")
    log:
        os.path.join(LOGDIR, "index_bam", "{sample}.log")
    threads: 4
    shell:
        r"""
        mkdir -p "$(dirname "{log}")" "{TMPDIR}/bam_indexed"
        (
        set -euo pipefail
        bam="{input.bam}"
        bam_nosuffix="${{bam%.bam}}"

        {SAMTOOLS} index -@ {threads} "$bam"

        if [[ -f "${{bam}}.bai" || -f "${{bam_nosuffix}}.bai" || -f "${{bam}}.csi" ]]; then
          touch "{output.ok}"
        else
          echo "ERROR: BAM index not found after indexing: $bam" >&2
          exit 1
        fi
        ) > "{log}" 2>&1
        """

rule contig2genome:
    # Assign each contig in a sample to a Microcoleus genome bin using DIAMOND
    # blastx (nucleotide contigs vs amino-acid marker proteins).
    # Output: two-column TSV  contig_id <TAB> genome_id  (best hit only).
    input:
        dmnd=DB_DMND,
        fasta=lambda wc: assembly_fasta_path(wc.sample),
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        os.path.join(TMPDIR, "contig2genome", "{sample}.contig2genome.tsv")
    log:
        os.path.join(LOGDIR, "contig2genome", "{sample}.log")
    threads: 8
    params:
        evalue=lambda wc: config["diamond"]["evalue"],
        qcov=lambda wc: config["diamond"]["query_cover"],
        max_t=lambda wc: config["diamond"]["max_target_seqs"]
    shell:
        r"""
        mkdir -p "$(dirname "{log}")" "{TMPDIR}/contig2genome" "{TMPDIR}/diamond_hits"
        (
        set -euo pipefail
        m8="{TMPDIR}/diamond_hits/{wildcards.sample}.m8"

        {DIAMOND} blastx \
          -d "{DB_PREFIX}" \
          -q "{input.fasta}" \
          -o "$m8" \
          --outfmt 6 qseqid sseqid bitscore evalue pident length \
          --max-target-seqs {params.max_t} \
          --evalue {params.evalue} \
          --query-cover {params.qcov} \
          --threads {threads}

        # Keep only the single best-scoring genome hit per contig
        awk 'BEGIN{{OFS="\t"}} {{
              q=$1; g=$2; b=$3;
              if(!(q in best) || b>best[q]) {{ best[q]=b; hit[q]=g; }}
            }}
            END{{ for(q in hit) print q, hit[q]; }}' "$m8" > "{output}"
        ) > "{log}" 2>&1
        """

rule genome_abundance:
    # Compute read-level abundance metrics for each Microcoleus genome bin
    # in each sample.  Joins:
    #   - samtools idxstats  → total mapped reads
    #   - samtools coverage  → per-contig depth and breadth
    #   - contig2genome map  → which genome each contig belongs to
    #
    # Outputs per genome:
    #   Reads, CPM, RPKM, MeanDepth (length-weighted), Breadth (length-weighted)
    #
    # Side-effect: writes {TMPDIR}/coverage/{sample}.coverage.tsv, which is
    # consumed later by rule ana_per_sample (not tracked as a Snakemake output
    # but guaranteed to exist once this rule completes).
    input:
        bam=lambda wc: bam_path(wc.sample),
        bam_ok=os.path.join(TMPDIR, "bam_indexed", "{sample}.ok"),
        c2g=os.path.join(TMPDIR, "contig2genome", "{sample}.contig2genome.tsv"),
        glist=GENOME_LIST,
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        os.path.join(PERSAMPLE, "{sample}.genome_abundance.tsv")
    log:
        os.path.join(LOGDIR, "genome_abundance", "{sample}.log")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")" "{TMPDIR}/coverage" "{TMPDIR}/idxstats"
        (
        set -euo pipefail

        bam="{input.bam}"
        sample="{wildcards.sample}"
        glist="{input.glist}"
        c2g="{input.c2g}"

        final_out="{output}"
        tmpout="${{final_out}}.tmp"
        tmptab="${{final_out}}.tmp.unsorted"

        if [[ ! -s "$glist" ]]; then
            echo "ERROR: genome list is empty: $glist" >&2
            exit 1
        fi

        idx="{TMPDIR}/idxstats/{wildcards.sample}.idxstats.tsv"
        cov="{TMPDIR}/coverage/{wildcards.sample}.coverage.tsv"

        {SAMTOOLS} idxstats "$bam" > "$idx"
        total_mapped=$(awk 'BEGIN{{s=0}} $1!="*" {{s+=$3}} END{{print s+0}}' "$idx")

        {SAMTOOLS} coverage -H "$bam" > "$cov"
        [[ -s "$cov" ]] || {{ echo "ERROR: coverage empty: $cov" >&2; exit 1; }}

        awk -v OFS="\t" -v sample="$sample" -v total="$total_mapped" -v glist="$glist" '
          BEGIN{{
            while((getline g < glist)>0){{ genomes[g]=1 }}
            close(glist)
          }}
          FNR==NR {{ c2g[$1]=$2; next }}
          NR==1 && $1=="#rname" {{ next }}
          {{
            contig=$1
            start=$2; end=$3
            len=end-start+1
            numreads=$4
            covpct=$6
            meandepth=$7

            if(!(contig in c2g)) next
            g=c2g[contig]

            reads[g]+=numreads
            clen[g]+=len
            contigs[g]+=1
            wdepth[g]+=meandepth*len
            wbreadth[g]+=(covpct/100.0)*len
          }}
          END{{
            print "GenomeID","Sample","Reads","CPM","RPKM","MeanDepth","Breadth","ContigsHit","SumContigLen_bp","TotalMappedReads"
            for(g in genomes){{
              r=(g in reads)?reads[g]:0
              L=(g in clen)?clen[g]:0
              n=(g in contigs)?contigs[g]:0
              cpm=(total>0)?(r*1e6/total):0
              rpkm=(L>0 && total>0)?(r*1e9/(L*total)):0
              md=(L>0 && (g in wdepth))?(wdepth[g]/L):0
              br=(L>0 && (g in wbreadth))?(wbreadth[g]/L):0
              print g, sample, r, cpm, rpkm, md, br, n, L, total
            }}
          }}' "$c2g" "$cov" > "$tmptab"

        head -n 1 "$tmptab" > "$tmpout"
        tail -n +2 "$tmptab" | sort -k1,1 >> "$tmpout"

        mv -f "$tmpout" "$final_out"
        rm -f "$tmptab"
        ) > "{log}" 2>&1
        """

rule matrices:
    # Pivot per-sample abundance TSVs into genome × sample matrices.
    # One matrix per metric: Reads, CPM, RPKM, MeanDepth, Breadth.
    input:
        expand(os.path.join(PERSAMPLE, "{sample}.genome_abundance.tsv"), sample=SAMPLES),
        glist=GENOME_LIST,
        dirs=os.path.join(OUTDIR, ".dirs.ok")
    output:
        reads  = os.path.join(MATRIXDIR, "genome_Reads_matrix.tsv"),
        cpm    = os.path.join(MATRIXDIR, "genome_CPM_matrix.tsv"),
        rpkm   = os.path.join(MATRIXDIR, "genome_RPKM_matrix.tsv"),
        depth  = os.path.join(MATRIXDIR, "genome_MeanDepth_matrix.tsv"),
        breadth= os.path.join(MATRIXDIR, "genome_Breadth_matrix.tsv")
    log:
        os.path.join(LOGDIR, "matrices.log")
    run:
        import csv, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        per_files = sorted(glob.glob(os.path.join(PERSAMPLE, "*.genome_abundance.tsv")))
        samples = [os.path.basename(x).replace(".genome_abundance.tsv","") for x in per_files]

        with open(input.glist) as f:
            genomes = [x.strip() for x in f if x.strip()]

        metrics = ["Reads","CPM","RPKM","MeanDepth","Breadth"]
        data = {m:{g:{s:"0" for s in samples} for g in genomes} for m in metrics}

        for fp, s in zip(per_files, samples):
            with open(fp) as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    g = row["GenomeID"]
                    if g not in data["Reads"]:
                        continue
                    for m in metrics:
                        data[m][g][s] = row[m]

        def write(metric, outpath):
            with open(outpath, "w") as out:
                out.write("GenomeID\t" + "\t".join(samples) + "\n")
                for g in genomes:
                    out.write(g)
                    for s in samples:
                        out.write("\t" + data[metric][g][s])
                    out.write("\n")

        write("Reads",   output.reads)
        write("CPM",     output.cpm)
        write("RPKM",    output.rpkm)
        write("MeanDepth", output.depth)
        write("Breadth", output.breadth)
        _lf.close(); _sys.stdout = _sys.__stdout__


# ════════════════════════════════════════════════════════════════════════════════
# ANATOXIN DETECTION RULES
# ════════════════════════════════════════════════════════════════════════════════
#
# Strategy
# ────────
# 1. Download the 10 core ana gene proteins (anaA–anaK) from NCBI, sourced
#    from reference producers: Oscillatoria PCC 6506, Cylindrospermum stagnale
#    PCC 7417, Microcoleus anatoxicus PTRS3, and related strains.
#    Headers are reformatted to  >gene|accession|organism  for easy parsing.
#
# 2. Build a DIAMOND protein database from those references.
#
# 3. For each sample: run DIAMOND blastx (nucleotide contigs vs ana proteins).
#    Each contig is assigned to its single best-matching ana gene.
#
# 4. Join blastx hits with the contig→genome map (from the existing pipeline)
#    and per-contig sequencing depth (from samtools coverage, already computed
#    by rule genome_abundance).
#
# 5. For every genome bin that carries at least one ana gene:
#      NormalizedCopies = length-weighted mean depth of ana-carrying contigs
#                         ─────────────────────────────────────────────────
#                         mean depth of the entire genome bin (all its contigs)
#    A value of ~1.0 means the ana cluster is present in roughly every copy
#    of that genome; <1.0 indicates partial presence (e.g. plasmid copy number
#    below chromosome copy number, or a mixed community).
#
# 6. Aggregate results into genome × sample matrices and fetch NCBI taxonomy
#    for any genome ID that encodes a GCA assembly accession.
#
# ════════════════════════════════════════════════════════════════════════════════

rule fetch_ana_refs:
    # Download anatoxin biosynthetic gene cluster proteins from NCBI.
    #
    # Uses organism-based searches (Kamptonema, Oscillatoria, Microcoleus,
    # Dolichospermum, Anabaena) rather than the [Gene Name] field query,
    # which is unreliable on NCBI's search backend.  Fetches GenBank format
    # records so that /gene= qualifiers (anaA, anaB, …) can be extracted
    # directly.  A keyword-based fallback maps product descriptions to gene
    # names for records lacking explicit /gene= annotations.
    #
    # Header format in output FASTA:  >gene|accession|organism
    # Requires outbound HTTPS to eutils.ncbi.nlm.nih.gov.
    output:
        ANA_REFS
    log:
        os.path.join(LOGDIR, "fetch_ana_refs.log")
    run:
        import urllib.request, urllib.parse, json, time, re, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"
        HDR  = {"User-Agent": "microcoleus-anatoxin-pipeline/1.0"}

        def ncbi_get(endpoint, params, retries=5):
            url = f"{NCBI}/{endpoint}?" + urllib.parse.urlencode(params)
            req = urllib.request.Request(url, headers=HDR)
            for attempt in range(retries):
                try:
                    with urllib.request.urlopen(req, timeout=60) as r:
                        return r.read().decode("utf-8", errors="replace").strip()
                except urllib.error.HTTPError as e:
                    if e.code in (429, 500, 503) and attempt < retries - 1:
                        time.sleep(10 * (attempt + 1))
                        continue
                    if attempt == retries - 1:
                        raise
                    time.sleep(5 * attempt)
                except Exception:
                    if attempt == retries - 1:
                        raise
                    time.sleep(5 * attempt)

        def esearch(query, db="protein", retmax=200):
            # strict=False: NCBI sometimes embeds \t/\n in JSON string values
            for attempt in range(5):
                r = ncbi_get("esearch.fcgi",
                             {"db": db, "term": query, "retmax": retmax, "retmode": "json"})
                if not r.startswith("{"):
                    time.sleep(8 * (attempt + 1))
                    continue
                try:
                    d = json.loads(r, strict=False)
                except json.JSONDecodeError:
                    time.sleep(8 * (attempt + 1))
                    continue
                res = d.get("esearchresult", {})
                if "ERROR" in res:
                    print(f"  NCBI error (attempt {attempt+1}): {res['ERROR'][:60]}", flush=True)
                    time.sleep(10 * (attempt + 1))
                    continue
                return res.get("idlist", [])
            print(f"  WARNING: all attempts failed for: {query}", flush=True)
            return []

        def efetch_gb(ids, db="protein"):
            """Fetch GenBank flat-file records in chunks of 50."""
            chunks = [ids[i:i+50] for i in range(0, len(ids), 50)]
            parts  = []
            for chunk in chunks:
                r = ncbi_get("efetch.fcgi",
                             {"db": db, "id": ",".join(chunk),
                              "rettype": "gb", "retmode": "text"})
                parts.append(r)
                time.sleep(0.5)
            return "\n".join(parts)

        def parse_gb(gb_text):
            """
            Parse GenBank protein flat files.
            Returns list of (accession, gene, product, organism, seq) tuples.
            Skips PDB-format records (chain IDs like '4IRN_A') which lack
            gene annotations and are not useful as alignment references.
            """
            records = []
            for block in gb_text.split("\n//\n"):
                block = block.strip()
                if not block:
                    continue
                m_acc = re.search(r"^ACCESSION\s+(\S+)", block, re.M)
                if not m_acc:
                    continue
                acc = m_acc.group(1)
                # Skip PDB chain records (e.g. 4IRN_A)
                if re.match(r"^\dxx\w{2}", acc, re.I) or "_" in acc:
                    continue
                m_org  = re.search(r"ORGANISM\s+(.+)",  block)
                m_gene = re.search(r'/gene="([^"]+)"',  block)
                m_prod = re.search(r'/product="([^"]+)"', block)
                m_seq  = re.search(r"ORIGIN\s+([\s\S]+)$", block)
                if not m_seq:
                    continue
                seq = re.sub(r"[\d\s/\\]", "", m_seq.group(1)).upper()
                if len(seq) < 50:   # ignore tiny fragments
                    continue
                records.append((
                    acc,
                    m_gene.group(1).lower() if m_gene else "",
                    m_prod.group(1)          if m_prod else "",
                    m_org.group(1).strip().replace(" ", "_") if m_org else "Unknown",
                    seq,
                ))
            return records

        # Keyword patterns for assigning gene names when /gene= is absent.
        # Order matters: more specific patterns first.
        KEYWORD_MAP = [
            (re.compile(r'\bana[Aa]\b|anatoxin.synthetase|non.ribosomal.peptide.synthetase.*atx', re.I), "anaA"),
            (re.compile(r'\bana[Bb]\b|prolyl.ACP.dehydrogenase|proline.oxidase|PLP.depend.*proline|histidine.ammonia.lyase', re.I), "anaB"),
            (re.compile(r'\bana[Cc]\b|adenylation.*(pro|ana)|proline.activating', re.I), "anaC"),
            (re.compile(r'\bana[Dd]\b|enoyl.reduct.*ana',              re.I), "anaD"),
            (re.compile(r'\bana[Ee]\b',                                 re.I), "anaE"),
            (re.compile(r'\bana[Ff]\b|acyl.carrier.*ana',               re.I), "anaF"),
            (re.compile(r'\bana[Gg]\b|FAD.*oxidoreductase.*ana',        re.I), "anaG"),
            (re.compile(r'\bana[Ii]\b|anatoxin.*transporter|transporter.*atx', re.I), "anaI"),
            (re.compile(r'\bana[Jj]\b|acetyltransferase.*ana',          re.I), "anaJ"),
            (re.compile(r'\bana[Kk]\b|F420.*oxidoreductase|dihydroanatoxin.*reduct', re.I), "anaK"),
        ]

        def infer_gene(gene_ann, product):
            combined = f"{gene_ann} {product}"
            # Explicit anaX match (case-insensitive). Uppercase the letter so
            # the result ("anaC") matches the ANA_GENES list keys exactly.
            m = re.search(r'\bana([abcdefgijk])\b', combined, re.I)
            if m:
                return "ana" + m.group(1).upper()
            # Keyword fallback (entries already use correct case)
            for pattern, gene in KEYWORD_MAP:
                if pattern.search(combined):
                    return gene
            return None

        # ── Organism-based searches ──────────────────────────────────────────
        # These bypass the unreliable [Gene Name] backend while still reaching
        # well-characterised anatoxin producers.
        ORG_QUERIES = [
            "Kamptonema[Organism] AND anatoxin",
            "Oscillatoria[Organism] AND anatoxin",
            "Microcoleus anatoxicus[Organism]",
            'Anabaena[Organism] AND "anatoxin synthetase"',
            "Dolichospermum[Organism] AND anatoxin synthetase",
            "Cylindrospermum[Organism] AND dihydroanatoxin",
            "Tychonema[Organism] AND anatoxin",
            "Sphaerospermopsis[Organism] AND anatoxin",
        ]

        all_ids = []
        for q in ORG_QUERIES:
            print(f"  Searching: {q}", flush=True)
            ids = esearch(q, retmax=200)
            all_ids.extend(ids)
            time.sleep(0.5)

        all_ids = list(dict.fromkeys(all_ids))   # deduplicate, preserve order
        print(f"  Unique protein IDs to fetch: {len(all_ids)}", flush=True)

        if not all_ids:
            raise RuntimeError(
                "No protein IDs found. Check network access to eutils.ncbi.nlm.nih.gov"
            )

        # ── Fetch GenBank records ────────────────────────────────────────────
        gb_text = efetch_gb(all_ids, "protein")
        parsed  = parse_gb(gb_text)
        print(f"  Parsed {len(parsed)} valid protein records", flush=True)

        # ── Assign gene names and write FASTA ────────────────────────────────
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        seen_acc    = set()
        gene_buckets = {g: [] for g in ANA_GENES}
        n_unassigned = 0

        for acc, gene_ann, product, organism, seq in parsed:
            if acc in seen_acc:
                continue
            seen_acc.add(acc)
            gene = infer_gene(gene_ann, product)
            if gene and gene in gene_buckets:
                gene_buckets[gene].append((acc, organism, seq))
            else:
                n_unassigned += 1

        print("  Gene assignment summary:", flush=True)
        for g in ANA_GENES:
            print(f"    {g}: {len(gene_buckets[g])} sequences", flush=True)
        print(f"    unassigned (excluded): {n_unassigned}", flush=True)

        # Atomic write
        total = 0
        tmp   = output[0] + ".tmp"
        with open(tmp, "w") as fh:
            for gene in ANA_GENES:
                for acc, org, seq in gene_buckets[gene]:
                    fh.write(f">{gene}|{acc}|{org}\n{seq}\n")
                    total += 1

        if total == 0:
            raise RuntimeError(
                "Gene assignment produced no sequences. "
                "Inspect NCBI responses or add organism queries."
            )

        os.rename(tmp, output[0])
        print(f"  Wrote {total} reference proteins → {output[0]}", flush=True)
        _lf.close(); _sys.stdout = _sys.__stdout__


rule make_ana_db:
    # Build a DIAMOND protein database from the downloaded ana reference proteins.
    input:
        ANA_REFS
    output:
        ANA_DB_DMND
    log:
        os.path.join(LOGDIR, "make_ana_db.log")
    shell:
        r"""
        mkdir -p "$(dirname "{log}")" "$(dirname "{ANA_DB_DMND}")"
        (
        set -euo pipefail
        {DIAMOND} makedb --in "{input}" -d "{ANA_DB_PFX}"
        ) > "{log}" 2>&1
        """


rule ana_blastx:
    # Translate metagenomic contigs in all 6 reading frames (blastx) and align
    # against the ana protein database.  Uses --sensitive mode to capture
    # divergent homologues (down to ~40% amino-acid identity, as seen for anaK
    # across Microcoleus / Cylindrospermum strains).
    #
    # Only the single best protein hit per contig is kept (--max-target-seqs 1).
    # The subject ID encodes the gene name in the first pipe-delimited field,
    # allowing downstream rules to know which ana gene each contig encodes.
    #
    # Output columns: qseqid sseqid pident length qcovhsp evalue bitscore
    input:
        db    = ANA_DB_DMND,
        fasta = lambda wc: assembly_fasta_path(wc.sample)
    output:
        os.path.join(ANA_DIR, "blastx", "{sample}.m8")
    log:
        os.path.join(LOGDIR, "ana_blastx", "{sample}.log")
    threads: 8
    shell:
        r"""
        mkdir -p "$(dirname "{log}")" "{ANA_DIR}/blastx"
        (
        set -euo pipefail
        {DIAMOND} blastx \
            -d "{ANA_DB_PFX}" \
            -q "{input.fasta}" \
            -o "{output}" \
            --outfmt 6 qseqid sseqid pident length qcovhsp evalue bitscore \
            --evalue 1e-5 \
            --query-cover 30 \
            --max-target-seqs 1 \
            --sensitive \
            --threads {threads}
        ) > "{log}" 2>&1
        """


rule ana_per_sample:
    # Core joining rule: for each sample, links blastx hits → genome bins → depth.
    #
    # Inputs reused from the existing pipeline (no recomputation needed):
    #   contig2genome  – maps each contig to a Microcoleus genome bin
    #   coverage.tsv   – per-contig mean sequencing depth (written by genome_abundance)
    #   MeanDepth matrix – genome-level mean depth per sample
    #
    # Output columns:
    #   Sample, GenomeID, GenesFound_N, GeneList, Completeness_pct,
    #   MeanAnaDepth_wtd, GenomeDepth, NormalizedCopies
    #
    # NormalizedCopies = length-weighted mean depth of ana-positive contigs
    #                    / mean depth of the whole genome bin
    # Interpretation: copies of the ana cluster per genome copy in this sample.
    input:
        hits   = os.path.join(ANA_DIR, "blastx", "{sample}.m8"),
        c2g    = os.path.join(TMPDIR, "contig2genome", "{sample}.contig2genome.tsv"),
        # Declare the per-sample abundance TSV as input to guarantee that
        # genome_abundance (which writes coverage.tsv as a side-effect) has run.
        ga     = os.path.join(PERSAMPLE, "{sample}.genome_abundance.tsv"),
        gdepth = os.path.join(MATRIXDIR, "genome_MeanDepth_matrix.tsv")
    output:
        os.path.join(ANA_DIR, "per_sample", "{sample}.ana_abundance.tsv")
    log:
        os.path.join(LOGDIR, "ana_per_sample", "{sample}.log")
    run:
        import csv, collections, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        sample   = wildcards.sample
        cov_path = os.path.join(TMPDIR, "coverage", f"{sample}.coverage.tsv")

        # 1. Parse blastx hits: keep first (best-scoring) hit per contig.
        #    Extract the gene name from the first pipe-delimited field of sseqid.
        contig_gene = {}
        with open(input.hits) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                contig  = parts[0]
                gene    = parts[1].split("|")[0]   # e.g. "anaA"
                pident  = float(parts[2])
                qcov    = float(parts[4])
                if contig not in contig_gene:
                    contig_gene[contig] = (gene, pident, qcov)

        # 2. Parse contig → genome bin mapping
        c2g = {}
        with open(input.c2g) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t")
                if len(parts) >= 2:
                    c2g[parts[0]] = parts[1]

        # 3. Parse per-contig coverage (samtools coverage -H output).
        #    Columns: #rname startpos endpos numreads covbases coverage meandepth ...
        cov = {}   # contig → (meandepth, length_bp)
        with open(cov_path) as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 7:
                    continue
                contig    = parts[0]
                length    = int(parts[2]) - int(parts[1]) + 1
                meandepth = float(parts[6])
                cov[contig] = (meandepth, length)

        # 4. Parse genome mean depth for this sample from the existing matrix
        genome_depth = {}
        with open(input.gdepth) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                try:
                    val = float(row.get(sample, 0) or 0)
                except ValueError:
                    val = 0.0
                genome_depth[row["GenomeID"]] = val

        # 5. Group ana hits by genome bin.
        #    genome → gene → list of (meandepth, length, contig, pident, qcov)
        genome_data = collections.defaultdict(
            lambda: collections.defaultdict(list)
        )
        for contig, (gene, pident, qcov) in contig_gene.items():
            genome = c2g.get(contig)
            if genome is None:
                continue   # contig not assigned to any Microcoleus genome
            depth, length = cov.get(contig, (0.0, 0))
            genome_data[genome][gene].append((depth, length, contig, pident, qcov))

        # 6. Compute per-genome summary and write output
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        cols = ["Sample", "GenomeID", "GenesFound_N", "GeneList",
                "Completeness_pct", "MeanAnaDepth_wtd", "GenomeDepth", "NormalizedCopies"]
        with open(output[0], "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for genome in sorted(genome_data):
                gene_dict   = genome_data[genome]
                genes_found = sorted(gene_dict)

                # Length-weighted mean depth across all ana-carrying contigs
                total_wtd  = sum(d * l
                                 for hits in gene_dict.values()
                                 for d, l, *_ in hits)
                total_len  = sum(l
                                 for hits in gene_dict.values()
                                 for d, l, *_ in hits)
                mean_depth = total_wtd / total_len if total_len > 0 else 0.0

                gdepth       = genome_depth.get(genome, 0.0)
                norm_copies  = mean_depth / gdepth if gdepth > 0 else 0.0
                completeness = 100.0 * len(genes_found) / N_ANA_GENES

                fh.write("\t".join([
                    sample, genome,
                    str(len(genes_found)),
                    ",".join(genes_found),
                    f"{completeness:.1f}",
                    f"{mean_depth:.4f}",
                    f"{gdepth:.4f}",
                    f"{norm_copies:.4f}",
                ]) + "\n")
        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_matrix:
    # Pivot per-sample anatoxin results into three genome × sample matrices:
    #
    #   ana_normalized_copies_matrix.tsv  – NormalizedCopies  (primary abundance metric)
    #   ana_completeness_matrix.tsv       – % of 10 ana genes detected
    #   ana_genes_present_matrix.tsv      – comma-separated list of genes detected
    #
    # Only genome bins with at least one ana hit in at least one sample are included.
    input:
        per_sample = expand(
            os.path.join(ANA_DIR, "per_sample", "{sample}.ana_abundance.tsv"),
            sample=SAMPLES
        )
    output:
        copies       = os.path.join(ANA_DIR, "matrices", "ana_normalized_copies_matrix.tsv"),
        completeness = os.path.join(ANA_DIR, "matrices", "ana_completeness_matrix.tsv"),
        genes        = os.path.join(ANA_DIR, "matrices", "ana_genes_present_matrix.tsv")
    log:
        os.path.join(LOGDIR, "ana_matrix.log")
    run:
        import csv, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        # Collect all rows from per-sample files
        all_rows = []
        for fp in input.per_sample:
            with open(fp) as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    all_rows.append(row)

        # Determine which genomes had any ana hit
        hit_genomes = sorted({r["GenomeID"] for r in all_rows})

        # Build per-genome per-sample value dictionaries
        copies_d  = {g: {s: "0"  for s in SAMPLES} for g in hit_genomes}
        comp_d    = {g: {s: "0"  for s in SAMPLES} for g in hit_genomes}
        genes_d   = {g: {s: ""   for s in SAMPLES} for g in hit_genomes}

        for row in all_rows:
            g = row["GenomeID"]
            s = row["Sample"]
            if g in copies_d and s in copies_d[g]:
                copies_d[g][s] = row["NormalizedCopies"]
                comp_d[g][s]   = row["Completeness_pct"]
                genes_d[g][s]  = row["GeneList"]

        os.makedirs(os.path.dirname(output.copies), exist_ok=True)

        def write_matrix(path, data):
            with open(path, "w") as fh:
                fh.write("GenomeID\t" + "\t".join(SAMPLES) + "\n")
                for g in hit_genomes:
                    vals = [data[g][s] for s in SAMPLES]
                    fh.write(g + "\t" + "\t".join(vals) + "\n")

        write_matrix(output.copies,       copies_d)
        write_matrix(output.completeness, comp_d)
        write_matrix(output.genes,        genes_d)
        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_taxonomy:
    # For each genome bin that tested positive for ana genes, retrieve the
    # organism name from NCBI.
    #
    # Genome IDs that embed a GCA assembly accession (e.g. "PH2017_39_GCA_020737835")
    # are queried against the NCBI Assembly database → organism name returned.
    # Genome IDs without a GCA accession (e.g. MAGs "POL7_B2") are reported
    # as-is with no species annotation.
    input:
        os.path.join(ANA_DIR, "matrices", "ana_normalized_copies_matrix.tsv")
    output:
        os.path.join(ANA_DIR, "taxonomy", "genome_taxonomy.tsv")
    log:
        os.path.join(LOGDIR, "ana_taxonomy.log")
    run:
        import csv, re, urllib.request, urllib.parse, json, time, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        NCBI = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

        def ncbi_get(endpoint, params, retries=4):
            url = f"{NCBI}/{endpoint}?" + urllib.parse.urlencode(params)
            req = urllib.request.Request(
                url,
                headers={"User-Agent": "microcoleus-anatoxin-pipeline/1.0 (contact: research)"}
            )
            for attempt in range(retries):
                try:
                    with urllib.request.urlopen(req, timeout=60) as r:
                        return r.read().decode("utf-8", errors="replace").strip()
                except urllib.error.HTTPError as e:
                    if e.code == 429:
                        time.sleep(10 * (attempt + 1))
                        continue
                    if attempt == retries - 1:
                        raise
                    time.sleep(3 ** attempt)
                except Exception:
                    if attempt == retries - 1:
                        raise
                    time.sleep(3 ** attempt)

        # Read genome IDs from the matrix
        genomes = []
        with open(input[0]) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                genomes.append(row["GenomeID"])

        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        with open(output[0], "w") as fh:
            fh.write("GenomeID\tSpecies\tGCA_Accession\n")
            for gid in genomes:
                m = re.search(r"GCA_(\d+)", gid)
                if not m:
                    # MAG or local bin: no GCA accession available
                    fh.write(f"{gid}\t\t\n")
                    continue
                gca = f"GCA_{m.group(1)}"
                try:
                    r = ncbi_get("esearch.fcgi",
                                 {"db": "assembly",
                                  "term": f"{gca}[Assembly Accession]",
                                  "retmax": 1, "retmode": "json"})
                    uids = json.loads(r, strict=False).get("esearchresult", {}).get("idlist", [])
                    time.sleep(0.4)
                    if not uids:
                        fh.write(f"{gid}\tNot found in NCBI\t{gca}\n")
                        continue
                    summ = ncbi_get("esummary.fcgi",
                                    {"db": "assembly", "id": uids[0], "retmode": "json"})
                    d   = json.loads(summ, strict=False)["result"][uids[0]]
                    org = d.get("organism", "Unknown")
                    time.sleep(0.4)
                    fh.write(f"{gid}\t{org}\t{gca}\n")
                except Exception as e:
                    fh.write(f"{gid}\tError:{e}\t{gca}\n")
        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_annotate_samples:
    # Join taxonomy names onto each per-sample anatoxin file.
    # Adds a 'Species' column immediately after 'GenomeID'.
    # Genome IDs embedding a GCA accession get the NCBI organism name;
    # custom isolate codes (e.g. Aus8_D4) are labelled 'Microcoleus sp.'
    # because they predate public deposition.
    # Output goes to per_sample_annotated/ to keep the raw files intact.
    input:
        tax = os.path.join(ANA_DIR, "taxonomy", "genome_taxonomy.tsv"),
        tsvs = expand(os.path.join(ANA_DIR, "per_sample", "{sample}.ana_abundance.tsv"),
                      sample=SAMPLES)
    output:
        expand(os.path.join(ANA_DIR, "per_sample_annotated", "{sample}.ana_abundance.tsv"),
               sample=SAMPLES)
    log:
        os.path.join(LOGDIR, "ana_annotate_samples.log")
    run:
        import csv, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        # Build genome → species mapping from the taxonomy TSV.
        # GCA-bearing IDs get the NCBI organism name; others get a fallback.
        tax = {}
        with open(input.tax) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                sp = row.get("Species", "").strip()
                tax[row["GenomeID"]] = sp if sp else "Microcoleus sp."

        out_dir = os.path.join(ANA_DIR, "per_sample_annotated")
        os.makedirs(out_dir, exist_ok=True)

        for tsv_path in input.tsvs:
            sample = os.path.basename(tsv_path).replace(".ana_abundance.tsv", "")
            out_path = os.path.join(out_dir, f"{sample}.ana_abundance.tsv")
            with open(tsv_path) as fin, open(out_path, "w") as fout:
                reader = csv.DictReader(fin, delimiter="\t")
                # Insert Species right after GenomeID
                orig_fields = list(reader.fieldnames)
                idx = orig_fields.index("GenomeID") + 1
                new_fields = orig_fields[:idx] + ["Species"] + orig_fields[idx:]
                writer = csv.DictWriter(fout, fieldnames=new_fields,
                                        delimiter="\t", extrasaction="ignore")
                writer.writeheader()
                for row in reader:
                    gid = row["GenomeID"]
                    row["Species"] = tax.get(gid, "Microcoleus sp.")
                    writer.writerow(row)
            print(f"  Annotated {sample}", flush=True)

        print(f"  Done — wrote {len(input.tsvs)} annotated per-sample files", flush=True)
        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_combined:
    # Concatenate all annotated per-sample anatoxin files into one flat table.
    # The header is written once; each sample's rows follow in sample name order.
    # This is the easiest file to load into R/Python for downstream analysis.
    input:
        expand(os.path.join(ANA_DIR, "per_sample_annotated", "{sample}.ana_abundance.tsv"),
               sample=SAMPLES)
    output:
        os.path.join(ANA_DIR, "ana_combined.tsv")
    log:
        os.path.join(LOGDIR, "ana_combined.log")
    run:
        import sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        header_written = False
        total_rows = 0
        with open(output[0], "w") as fout:
            for path in sorted(input):
                with open(path) as fin:
                    lines = fin.readlines()
                if not lines:
                    continue
                if not header_written:
                    fout.write(lines[0])   # header
                    header_written = True
                for line in lines[1:]:     # data rows
                    fout.write(line)
                    total_rows += 1

        print(f"  Combined {len(input)} sample files → {total_rows} rows", flush=True)
        print(f"  Output: {output[0]}", flush=True)
        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_summary:
    # Final output: one row per genome bin that carries the ana cluster.
    # Columns:
    #   GenomeID, Species, GCA_Accession,
    #   N_Samples_Detected,          – number of samples where NormalizedCopies > 0
    #   Max_Completeness_pct,        – highest fraction of the 10 genes detected
    #   Max_NormalizedCopies,        – peak copies-per-genome across all samples
    #   Mean_NormalizedCopies_det,   – mean over samples where the genome was detected
    #   NormCopies_{sample}…         – per-sample normalized copies (0 = absent)
    input:
        copies       = os.path.join(ANA_DIR, "matrices", "ana_normalized_copies_matrix.tsv"),
        completeness = os.path.join(ANA_DIR, "matrices", "ana_completeness_matrix.tsv"),
        genes        = os.path.join(ANA_DIR, "matrices", "ana_genes_present_matrix.tsv"),
        taxonomy     = os.path.join(ANA_DIR, "taxonomy", "genome_taxonomy.tsv")
    output:
        os.path.join(ANA_DIR, "ana_summary.tsv")
    log:
        os.path.join(LOGDIR, "ana_summary.log")
    run:
        import csv, sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        # Load taxonomy lookup
        tax = {}
        with open(input.taxonomy) as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                tax[row["GenomeID"]] = (row.get("Species", ""), row.get("GCA_Accession", ""))

        # Load matrices
        def load_matrix(path):
            rows = {}
            with open(path) as fh:
                reader = csv.DictReader(fh, delimiter="\t")
                sample_cols = [c for c in reader.fieldnames if c != "GenomeID"]
                for row in reader:
                    rows[row["GenomeID"]] = row
            return rows, sample_cols

        copies_rows, sample_cols = load_matrix(input.copies)
        comp_rows, _             = load_matrix(input.completeness)
        genes_rows, _            = load_matrix(input.genes)

        with open(output[0], "w") as fh:
            header = (["GenomeID", "Species", "GCA_Accession",
                       "N_Samples_Detected", "Max_Completeness_pct",
                       "Max_NormalizedCopies", "Mean_NormalizedCopies_det",
                       "GenesDetected_anysample"] +
                      [f"NormCopies_{s}" for s in sample_cols])
            fh.write("\t".join(header) + "\n")

            for gid in sorted(copies_rows):
                row_c    = copies_rows[gid]
                row_comp = comp_rows.get(gid, {})
                row_g    = genes_rows.get(gid, {})
                species, gca = tax.get(gid, ("", ""))

                vals      = [float(row_c.get(s, 0) or 0) for s in sample_cols]
                comp_vals = [float(row_comp.get(s, 0) or 0) for s in sample_cols]

                detected     = [v for v in vals if v > 0]
                n_detected   = len(detected)
                max_copies   = max(vals) if vals else 0
                mean_det     = sum(detected) / len(detected) if detected else 0
                max_comp     = max(comp_vals) if comp_vals else 0

                # Union of genes detected across all samples
                all_genes = sorted({
                    g for s in sample_cols
                    for g in (row_g.get(s, "") or "").split(",")
                    if g
                })

                fh.write("\t".join([
                    gid, species, gca,
                    str(n_detected),
                    f"{max_comp:.1f}",
                    f"{max_copies:.4f}",
                    f"{mean_det:.4f}",
                    ",".join(all_genes),
                ] + [f"{v:.4f}" for v in vals]) + "\n")
        _lf.close(); _sys.stdout = _sys.__stdout__


# ════════════════════════════════════════════════════════════════════════════════
# COMMUNITY CORRELATION RULES
# ════════════════════════════════════════════════════════════════════════════════
#
# Uses SingleM contig-level taxonomic profiles (one TSV per sample) to test
# whether anatoxin abundance correlates with community composition.
#
# singlem_matrix      – pivots 450 profile TSVs into a sample × phylum matrix
#                       (raw coverage and fraction-of-bacteria) plus a
#                       dedicated Cyanobacteriota summary file.
#
# ana_singlem_correlation – joins the anatoxin NormalizedCopies matrix with
#                           the SingleM phylum matrix and computes Spearman ρ
#                           for every (genome, phylum) pair, plus an aggregate
#                           "any anatoxin in sample" series.
#
# Run with:  snakemake correlation_all --cores 4
# ════════════════════════════════════════════════════════════════════════════════

rule correlation_all:
    input:
        os.path.join(CORRDIR, "ana_singlem_correlations.tsv"),
        os.path.join(CORRDIR, "ana_singlem_merged.tsv"),
        os.path.join(CORRDIR, "singlem_cyanobacteria.tsv")


rule singlem_matrix:
    # Pivot all per-sample SingleM profile TSVs into two sample × phylum matrices:
    #   singlem_phylum_coverage.tsv  – raw SingleM coverage values
    #   singlem_phylum_fraction.tsv  – fraction of total bacterial coverage
    #   singlem_cyanobacteria.tsv    – Cyanobacteriota coverage + fraction only
    #
    # Extracts phylum-level rows (taxonomy depth = 3: Root; d__X; p__Y).
    # Strips the _noOrganellar suffix SingleM appends to sample names.
    input:
        profiles = [os.path.join(SINGLEM_DIR, f"{s}_singlem_contigs_profile.tsv")
                    for s in SAMPLES]
    output:
        cov   = os.path.join(CORRDIR, "singlem_phylum_coverage.tsv"),
        frac  = os.path.join(CORRDIR, "singlem_phylum_fraction.tsv"),
        cyano = os.path.join(CORRDIR, "singlem_cyanobacteria.tsv")
    log:
        os.path.join(LOGDIR, "singlem_matrix.log")
    run:
        import sys as _sys
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        cov_data = {}   # sample → {phylum: coverage}
        bact_tot = {}   # sample → total bacterial coverage

        for path in input.profiles:
            fname  = os.path.basename(path)
            sample = fname.replace("_singlem_contigs_profile.tsv", "") \
                          .replace("_noOrganellar", "")
            cov_data[sample] = {}
            bact_tot[sample] = 0.0

            with open(path) as fh:
                next(fh)   # skip header
                for line in fh:
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 3:
                        continue
                    cov = float(parts[1])
                    tax = parts[2]
                    if tax == "Root; d__Bacteria":
                        bact_tot[sample] = cov
                    if tax.count(";") == 2:   # phylum-level row exactly
                        phylum = tax.split("; ")[-1]
                        cov_data[sample][phylum] = cov

        all_phyla   = sorted({p for d in cov_data.values() for p in d})
        all_samples = sorted(cov_data.keys())
        print(f"  Samples: {len(all_samples)}, Phyla: {len(all_phyla)}", flush=True)

        def write_matrix(path, getter):
            os.makedirs(os.path.dirname(path), exist_ok=True)
            with open(path, "w") as fh:
                fh.write("Sample\t" + "\t".join(all_phyla) + "\n")
                for s in all_samples:
                    vals = [f"{getter(s, p):.6f}" for p in all_phyla]
                    fh.write(s + "\t" + "\t".join(vals) + "\n")

        write_matrix(output.cov,  lambda s, p: cov_data[s].get(p, 0.0))
        write_matrix(output.frac, lambda s, p: (cov_data[s].get(p, 0.0) / bact_tot[s]
                                                if bact_tot[s] > 0 else 0.0))

        cyano_key = "p__Cyanobacteriota"
        with open(output.cyano, "w") as fh:
            fh.write("Sample\tCyanobacteriota_cov\tTotal_bacteria_cov\tCyanobacteriota_frac\n")
            for s in all_samples:
                cyano = cov_data[s].get(cyano_key, 0.0)
                tot   = bact_tot[s]
                frac  = cyano / tot if tot > 0 else 0.0
                fh.write(f"{s}\t{cyano:.4f}\t{tot:.4f}\t{frac:.6f}\n")

        _lf.close(); _sys.stdout = _sys.__stdout__


rule ana_singlem_correlation:
    # Join anatoxin NormalizedCopies with SingleM phylum coverage and compute
    # Spearman rank correlations.
    #
    # Two correlation series:
    #   genome  – per ana-positive genome vs every phylum
    #   aggregate – max NormalizedCopies across all genomes per sample
    #               (i.e. "any anatoxin detected") vs every phylum
    #
    # Results are sorted by |ρ| descending so the strongest associations
    # appear first regardless of direction.
    #
    # Also writes ana_singlem_merged.tsv – a flat table with anatoxin copies
    # AND community composition side-by-side, ready for R/Python.
    input:
        copies = os.path.join(ANA_DIR, "matrices", "ana_normalized_copies_matrix.tsv"),
        cov    = os.path.join(CORRDIR, "singlem_phylum_coverage.tsv"),
        frac   = os.path.join(CORRDIR, "singlem_phylum_fraction.tsv")
    output:
        merged = os.path.join(CORRDIR, "ana_singlem_merged.tsv"),
        corr   = os.path.join(CORRDIR, "ana_singlem_correlations.tsv")
    log:
        os.path.join(LOGDIR, "ana_singlem_correlation.log")
    run:
        import csv, sys as _sys
        from scipy.stats import spearmanr
        os.makedirs(os.path.dirname(log[0]), exist_ok=True)
        _lf = open(log[0], "w"); _sys.stdout = _lf

        # Anatoxin copies matrix: genome → {sample: NormalizedCopies}
        copies, samples_ana = {}, []
        with open(input.copies) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            samples_ana = [c for c in reader.fieldnames if c != "GenomeID"]
            for row in reader:
                copies[row["GenomeID"]] = {s: float(row[s]) for s in samples_ana}

        # SingleM phylum coverage matrix: sample → {phylum: coverage}
        singlem_cov, singlem_frac, phyla = {}, {}, []
        with open(input.cov) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            phyla = [c for c in reader.fieldnames if c != "Sample"]
            for row in reader:
                singlem_cov[row["Sample"]] = {p: float(row[p]) for p in phyla}
        with open(input.frac) as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                singlem_frac[row["Sample"]] = {p: float(row[p]) for p in phyla}

        shared = sorted(set(samples_ana) & set(singlem_cov))
        print(f"  Shared samples: {len(shared)}, ana genomes: {len(copies)}, phyla: {len(phyla)}",
              flush=True)

        # ── Merged flat table (one row per sample × genome) ───────────────────
        cyano_col = "p__Cyanobacteriota"
        os.makedirs(os.path.dirname(output.merged), exist_ok=True)
        with open(output.merged, "w") as fh:
            header = (["Sample", "GenomeID", "NormalizedCopies", "Ana_detected",
                       "Cyano_cov", "Cyano_frac"] +
                      [f"cov_{p}" for p in phyla])
            fh.write("\t".join(header) + "\n")
            for genome in sorted(copies):
                for s in shared:
                    nc   = copies[genome].get(s, 0.0)
                    cyano = singlem_cov[s].get(cyano_col, 0.0)
                    cyano_fr = singlem_frac[s].get(cyano_col, 0.0)
                    cov_vals = [f"{singlem_cov[s].get(p, 0.0):.4f}" for p in phyla]
                    fh.write("\t".join([s, genome, f"{nc:.4f}",
                                        "1" if nc > 0 else "0",
                                        f"{cyano:.4f}", f"{cyano_fr:.6f}"]
                                       + cov_vals) + "\n")

        # ── Spearman correlations ─────────────────────────────────────────────
        corr_rows = []

        # Per-genome: only for genomes with at least one detection
        for genome in sorted(copies):
            nc_vec = [copies[genome].get(s, 0.0) for s in shared]
            n_det  = sum(1 for x in nc_vec if x > 0)
            if n_det == 0:
                continue
            for p in phyla:
                cov_vec = [singlem_cov[s].get(p, 0.0) for s in shared]
                rho, pval = spearmanr(nc_vec, cov_vec)
                corr_rows.append((genome, p, rho, pval, len(shared), n_det, "genome"))

        # Aggregate: max NormalizedCopies across all ana genomes per sample
        max_nc   = {s: max((copies[g].get(s, 0.0) for g in copies), default=0.0)
                    for s in shared}
        n_det_any = sum(1 for v in max_nc.values() if v > 0)
        for p in phyla:
            cov_vec = [singlem_cov[s].get(p, 0.0) for s in shared]
            nc_vec  = [max_nc[s] for s in shared]
            rho, pval = spearmanr(nc_vec, cov_vec)
            corr_rows.append(("ANY_genome", p, rho, pval, len(shared), n_det_any, "aggregate"))

        corr_rows.sort(key=lambda r: abs(r[2]), reverse=True)

        with open(output.corr, "w") as fh:
            fh.write("GenomeID\tPhylum\tSpearman_rho\tP_value\t"
                     "N_samples\tN_detected\tLevel\n")
            for genome, p, rho, pval, n, nd, lvl in corr_rows:
                fh.write(f"{genome}\t{p}\t{rho:.4f}\t{pval:.4e}\t{n}\t{nd}\t{lvl}\n")

        print(f"  Wrote {len(corr_rows)} correlation rows to {output.corr}", flush=True)
        _lf.close(); _sys.stdout = _sys.__stdout__
