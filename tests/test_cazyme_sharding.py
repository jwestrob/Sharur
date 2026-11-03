"""Sharded dbCAN consensus equals the single-process run (scripts/classify_cazymes.py).

Stand-in ``diamond`` and ``aksha`` executables keep the property that matters: the HMM
i-evalue is P x Z, with Z the number of target sequences in the search, so splitting
the proteins of an HMM search would change calls near the 1e-15 cut.
"""

import importlib.util
import json
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT / "scripts"))  # merge workers import the module by name
spec = importlib.util.spec_from_file_location("classify_cazymes", REPO_ROOT / "scripts/classify_cazymes.py")
cc = importlib.util.module_from_spec(spec)
sys.modules["classify_cazymes"] = cc
spec.loader.exec_module(cc)

# PL1 has dbCAN.hmm profiles only: its value reaches the table from one tool, so a tie
# between its two profiles (in different parts) must resolve as in the single table
MOTIFS = {"GH5": "WWGH", "GT2": "CCGT", "CE4": "KKCE", "PL1": "PPLY"}

FAKE_DIAMOND = r'''
import sys
a = sys.argv[1:]
cut = float(a[a.index("--evalue") + 1])
MOTIFS = {MOTIFS}
def emit(name, seq):
    for fam, motif in MOTIFS.items():
        k = seq.count(motif) if fam != "PL1" else 0
        ev = float("%.2e" % (10.0 ** -(60 + 30 * k)))
        if k and ev <= cut:
            print(f"{{name}}\tREF{{k}}|{{fam}}_1\t90.0\t{{ev:.2e}}\t{{300 + k}}")
name, seq = None, []
for line in open(a[a.index("--query") + 1]):
    if line.startswith(">"):
        if name:
            emit(name, "".join(seq))
        name, seq = line[1:].split()[0], []
    else:
        seq.append(line.strip())
if name:
    emit(name, "".join(seq))
'''

FAKE_AKSHA = r'''
import re, sys
from pathlib import Path
a = sys.argv[1:]
assert a[0] == "search"
prot, out = Path(a[a.index("--prot_in") + 1]), Path(a[a.index("--outdir") + 1])
profiles = re.findall(r"^NAME\s+(\S+)\s*\nLENG\s+(\d+)", open(a[a.index("--hmm_in") + 1]).read(), re.M)
seqs = {{}}
for faa in sorted(prot.glob("*.faa")):
    name = None
    for line in open(faa):
        if line.startswith(">"):
            name = line[1:].split()[0]; seqs[name] = ""
        else:
            seqs[name] += line.strip().rstrip("*")
Z = len(seqs)  # pyhmmer's default: the number of target sequences searched
MOTIFS = {MOTIFS}
with open(out / "user_hmms_hits_df.tsv", "w") as fh:
    fh.write("sequence_id\thmm_name\tbitscore\tevalue\tc_evalue\ti_evalue\tenv_from\tenv_to\t"
             "dom_bitscore\tali_from\tali_to\thmm_from\thmm_to\n")
    for prof, leng in profiles:
        fam = re.match(r"(GH|GT|CE|PL)\d+", prof).group(0)
        rows = []
        for sid, seq in seqs.items():
            k = seq.count(MOTIFS[fam]) if fam in MOTIFS else 0
            if not k:
                continue
            p = 10.0 ** -(16 + k) * (1 + len(seq) % 7)
            # same P for every profile of a family (i-evalue ties); the score is the profile's own
            score = 50 + 10 * k + sum(map(ord, prof)) % 7
            hmm_to = int(leng) // 2 if len(seq) % 5 else int(leng) // 5  # some fail coverage
            rows.append((-score, sid, p, score, hmm_to))
        for _, sid, p, score, hmm_to in sorted(rows):
            fh.write(f"{{sid}}\t{{prof}}\t{{score + 1:.2f}}\t{{p * Z:.2e}}\t{{p * 3:.2e}}\t{{p * Z:.2e}}\t"
                     f"1\t{{len(seqs[sid])}}\t{{score:.2f}}\t1\t{{len(seqs[sid])}}\t1\t{{hmm_to}}\n")
'''


def _hmm_file(path: Path, profiles: list[tuple[str, int]]) -> None:
    path.write_text("".join(f"HMMER3/f [fake]\nNAME  {n}\nLENG  {l}\nHMM   A C D\n//\n" for n, l in profiles))


@pytest.fixture
def toy(tmp_path, monkeypatch):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    for name, body in (("diamond", FAKE_DIAMOND), ("aksha", FAKE_AKSHA)):
        exe = bin_dir / name
        exe.write_text(f"#!{sys.executable}\n" + body.format(MOTIFS=repr(MOTIFS)))
        exe.chmod(0o755)
    monkeypatch.setenv("PATH", f"{bin_dir}:{Path(sys.executable).parent}:/usr/bin:/bin")

    dbcan = tmp_path / "dbcan"
    dbcan.mkdir()
    (dbcan / "CAZy.dmnd").write_text("")
    # profiles of one family sit in different parts, so cross-part ties are exercised
    # with 3 parts, PL1.hmm sits below the GH5.hmm rows of part 0 while PL1_2.hmm heads part 2
    # (GT99.hmm has no hits): their tie must resolve by file order, not by row within a part
    _hmm_file(dbcan / "dbCAN.hmm", [("GH5.hmm", 100), ("GT2.hmm", 110), ("GT99.hmm", 70), ("PL1.hmm", 60),
                                    ("GH5_2.hmm", 100), ("PL1_2.hmm", 60), ("GT2_1.hmm", 100),
                                    ("CE4.hmm", 90), ("CE4_3.hmm", 95)])
    _hmm_file(dbcan / "dbCAN_sub.hmm", [("GH5_e1.hmm|GH5:3", 80), ("GT2_e7.hmm|GT2:4", 120),
                                        ("GH5_e2.hmm|GH5:1", 90), ("CE4_e1.hmm|CE4:2", 70),
                                        ("GT2_e9.hmm|GT2:2", 100)])

    # 90 proteins carrying zero to two copies of each motif, varied lengths
    state, records = 7, []
    for i in range(90):
        parts = []
        for motif in MOTIFS.values():
            state = (state * 1103515245 + 12345) % 2**31
            parts += [motif] * (state % 3 if state % 4 else 0)
        state = (state * 1103515245 + 12345) % 2**31
        filler = "AL" * (20 + state % 23)
        records.append(f">p{i:03d}\n{filler}{''.join(parts)}{filler[:state % 9]}*\n")
    faa = tmp_path / "proteins.faa"
    faa.write_text("".join(records))
    return faa, dbcan, tmp_path


def test_sharded_run_equals_single_process(toy):
    faa, dbcan, tmp = toy
    single = tmp / "single"
    single.mkdir()
    ref = cc.classify_cazymes(None, threads=1, verbose=False, dbcan_dir=dbcan, write_db=False,
                              results_out=single / "cazy.tsv", proteins=str(faa))
    assert set(ref["n_tools"]) == {2, 3}

    for n in (1, 3, 4):
        work, out = tmp / f"work{n}", tmp / f"sharded{n}"
        out.mkdir()
        for stage in ("search", "sub"):
            for i in range(n):
                cc.run_stage(stage, None, str(faa), work, n, shard_index=i, threads=1,
                             dbcan_dir=dbcan, verbose=False)
        # the merge parses shard tables in worker processes when given threads
        cc.run_stage("merge", None, str(faa), work, n, threads=1 if n == 1 else 2, dbcan_dir=dbcan,
                     results_out=out / "cazy.tsv", verbose=False)
        # side files carry the --out stem; the reusable ones carry provenance stamps
        names = ["cazy.tsv"] + [f"cazy.{t}" for t in (cc.DIAMOND_CACHE, cc.HMM_TABLE, cc.SUB_TABLE)]
        names += [f"cazy.{t}.stamp.json" for t in (cc.DIAMOND_CACHE, cc.HMM_TABLE)]
        for name in names:
            assert (out / name).read_text() == (single / name).read_text(), (n, name)
    assert "family_tools" in ref.columns


def _run_shards(faa, dbcan, work, n, stages=("search", "sub")):
    for stage in stages:
        for i in range(n):
            cc.run_stage(stage, None, str(faa), work, n, shard_index=i, threads=1, dbcan_dir=dbcan,
                         verbose=False)


def _install(bin_dir: Path, name: str, body: str) -> None:
    exe = bin_dir / name
    exe.write_text(f"#!{sys.executable}\n" + body.format(MOTIFS=repr(MOTIFS)))
    exe.chmod(0o755)


def test_merge_raises_when_a_reported_sub_table_is_missing(toy):
    """Finding 1: a sub shard that reports a search but lost its table stops the merge."""
    faa, dbcan, tmp = toy
    work = tmp / "work"
    _run_shards(faa, dbcan, work, 3)
    (work / "sub_0001" / "dbcan_sub.tsv").unlink()
    with pytest.raises(RuntimeError, match="dbCAN-sub shard table"):
        cc.run_stage("merge", None, str(faa), work, 3, threads=1, dbcan_dir=dbcan,
                     results_out=tmp / "cazy.tsv", verbose=False)


def test_merge_raises_on_sub_shards_from_other_search_results(toy):
    """Finding 1: sub shards built on an earlier search (here a sparser DIAMOND) stop the
    merge once the search stage is rerun, instead of losing DIAMOND + dbCAN-sub proteins."""
    faa, dbcan, tmp = toy
    bin_dir, work = tmp / "bin", tmp / "work"
    _install(bin_dir, "diamond", FAKE_DIAMOND.replace("if k and ev <= cut:", "if k >= 3 and ev <= cut:"))
    _run_shards(faa, dbcan, work, 3)
    _install(bin_dir, "diamond", FAKE_DIAMOND)
    _run_shards(faa, dbcan, work, 3, stages=("search",))
    with pytest.raises(RuntimeError, match="rerun the sub stage"):
        cc.run_stage("merge", None, str(faa), work, 3, threads=1, dbcan_dir=dbcan,
                     results_out=tmp / "cazy.tsv", verbose=False)
    _run_shards(faa, dbcan, work, 3, stages=("sub",))  # the rerun makes the merge whole again
    merged = cc.run_stage("merge", None, str(faa), work, 3, threads=1, dbcan_dir=dbcan,
                          results_out=tmp / "cazy.tsv", verbose=False)
    (tmp / "one").mkdir()
    single = cc.classify_cazymes(None, threads=1, verbose=False, dbcan_dir=dbcan, write_db=False,
                                 results_out=tmp / "one" / "cazy.tsv", proteins=str(faa))
    assert merged.equals(single)


def test_merge_raises_when_sub_shards_disagree_on_candidates(toy):
    faa, dbcan, tmp = toy
    work = tmp / "work"
    _run_shards(faa, dbcan, work, 2)
    done = work / "sub_0001" / "DONE.json"
    info = json.loads(done.read_text())
    info["candidates_hash"] = "0" * 32
    done.write_text(json.dumps(info))
    with pytest.raises(RuntimeError, match="disagree on candidates_hash"):
        cc.run_stage("merge", None, str(faa), work, 2, threads=1, dbcan_dir=dbcan,
                     results_out=tmp / "cazy.tsv", verbose=False)


def test_missing_diamond_raises_single_and_sharded(toy, monkeypatch):
    """Finding 3: no diamond binary stops the run instead of yielding an HMM-only consensus."""
    faa, dbcan, tmp = toy
    (tmp / "bin" / "diamond").unlink()
    monkeypatch.setenv("PATH", f"{tmp / 'bin'}:/usr/bin:/bin")
    (tmp / "single").mkdir()
    with pytest.raises(RuntimeError, match="diamond"):
        cc.classify_cazymes(None, threads=1, verbose=False, dbcan_dir=dbcan, write_db=False,
                            results_out=tmp / "single" / "cazy.tsv", proteins=str(faa))
    work = tmp / "work"
    with pytest.raises(RuntimeError, match="diamond"):
        cc.run_stage("search", None, str(faa), work, 2, shard_index=0, threads=1, dbcan_dir=dbcan,
                     verbose=False)
    assert not (work / "search_0000" / "DONE.json").exists()


def _second_proteome(tmp: Path) -> Path:
    """The toy proteome with renamed proteins and one motif copy fewer in every fifth protein."""
    lines = (tmp / "proteins.faa").read_text().splitlines(keepends=True)
    out = []
    for k in range(0, len(lines), 2):
        seq = lines[k + 1]
        if k % 10 == 0:
            for motif in MOTIFS.values():
                if motif in seq:
                    seq = seq.replace(motif, "", 1)
                    break
        out += [lines[k].replace(">p", ">q"), seq]
    faa = tmp / "other.faa"
    faa.write_text("".join(out))
    return faa


def test_runs_sharing_an_output_directory_keep_their_own_evidence(toy):
    """Finding 2: two inputs written into one directory get stem-prefixed side files, and
    each consensus equals the one it gets alone."""
    faa, dbcan, tmp = toy
    other = _second_proteome(tmp)
    shared, alone = tmp / "shared", tmp / "alone"
    shared.mkdir()
    alone.mkdir()
    kw = dict(threads=1, verbose=False, dbcan_dir=dbcan, write_db=False)
    cc.classify_cazymes(None, results_out=shared / "a.tsv", proteins=str(faa), **kw)
    b = cc.classify_cazymes(None, results_out=shared / "b.tsv", proteins=str(other), **kw)
    b_alone = cc.classify_cazymes(None, results_out=alone / "b.tsv", proteins=str(other), **kw)
    assert b.equals(b_alone) and "diamond" in ",".join(b["tools"])
    assert (shared / f"a.{cc.DIAMOND_CACHE}").exists() and (shared / f"b.{cc.DIAMOND_CACHE}").exists()


def test_caches_are_reused_only_with_a_matching_stamp(toy, capsys):
    """Finding 2: a DIAMOND cache or dbCAN.hmm table stamped for other proteins, another
    CAZy.dmnd, or carrying no stamp is ignored with a log line; a matching one is reused."""
    faa, dbcan, tmp = toy
    out = tmp / "out"
    out.mkdir()
    kw = dict(threads=1, verbose=False, dbcan_dir=dbcan, write_db=False, results_out=out / "cazy.tsv")
    ref = cc.classify_cazymes(None, proteins=str(faa), **kw)
    cache, table = out / f"cazy.{cc.DIAMOND_CACHE}", out / f"cazy.{cc.HMM_TABLE}"
    capsys.readouterr()

    # a matching stamp: both are reused, and the result is unchanged
    assert cc.classify_cazymes(None, proteins=str(faa), **kw).equals(ref)
    log = capsys.readouterr().out
    assert "Ignoring" not in log and "Reusing dbCAN.hmm table" in log

    # another proteome: both caches are ignored and the result equals a clean run
    other = _second_proteome(tmp)
    got = cc.classify_cazymes(None, proteins=str(other), **kw)
    log = capsys.readouterr().out
    assert "Ignoring DIAMOND cache" in log and "stamp differs in proteins" in log
    (tmp / "clean").mkdir()
    clean = cc.classify_cazymes(None, proteins=str(other), **{**kw, "results_out": tmp / "clean" / "cazy.tsv"})
    assert got.equals(clean)

    # another CAZy.dmnd: the DIAMOND cache is ignored
    (dbcan / "CAZy.dmnd").write_text("another release")
    cc.classify_cazymes(None, proteins=str(other), **kw)
    assert "stamp differs in reference" in capsys.readouterr().out

    # an unstamped cache (one written before stamps existed) is ignored
    (out / f"cazy.{cc.DIAMOND_CACHE}.stamp.json").unlink()
    cc.classify_cazymes(None, proteins=str(other), **kw)
    assert "carries no provenance stamp" in capsys.readouterr().out
    assert cache.exists() and table.exists()


def test_both_dbcan_sub_file_names_stop_the_run(tmp_path):
    """Finding 7: dbCAN-sub.hmm beside dbCAN_sub.hmm is ambiguous between releases."""
    for name in ("dbCAN.hmm", "dbCAN-sub.hmm", "dbCAN_sub.hmm"):
        (tmp_path / name).write_text("x")
    with pytest.raises(ValueError, match="both dbCAN-sub.hmm and dbCAN_sub.hmm"):
        cc.dbcan_hmm_files(tmp_path)
    (tmp_path / "dbCAN-sub.hmm").unlink()
    assert cc.dbcan_hmm_files(tmp_path)[1] == tmp_path / "dbCAN_sub.hmm"


def test_protein_split_hmm_search_changes_calls(toy):
    """The control: the same HMM search split by protein reports smaller i-evalues."""
    faa, dbcan, tmp = toy
    prot_dir, _, _ = cc.resolve_proteins(None, str(faa), verbose=False)
    whole = cc.run_astra_hmm(prot_dir, "whole", hmm_path=dbcan / "dbCAN.hmm", threads=1, verbose=False)
    lines = faa.read_text().splitlines(keepends=True)
    split: dict = {}
    for i in range(3):
        d = tmp / f"chunk{i}"
        d.mkdir()
        (d / "p.faa").write_text("".join(lines[2 * j] + lines[2 * j + 1] for j in range(i, len(lines) // 2, 3)))
        split.update(cc.run_astra_hmm(str(d), "chunk", hmm_path=dbcan / "dbCAN.hmm", threads=1, verbose=False))
    assert set(split) > set(whole)


def test_profile_parts_cover_every_profile_once(tmp_path):
    hmm = tmp_path / "x.hmm"
    names = [f"GH{i}.hmm" for i in range(10)]
    _hmm_file(hmm, [(name, 50 + i) for i, name in enumerate(names)])
    record = {name: f"HMMER3/f [fake]\nNAME  {name}\nLENG  {50 + i}\nHMM   A C D\n//\n" for i, name in enumerate(names)}
    assert cc.hmm_profile_order(hmm) == {name: i for i, name in enumerate(names)}
    for n in (1, 3, 4, 10, 12):
        for p in range(n):
            out = tmp_path / f"part{n}_{p}.hmm"
            assert cc.write_hmm_part(hmm, p, n, out) == len(names[p::n])
            assert out.read_text() == "".join(record[name] for name in names[p::n])


def test_run_dbcan_z_makes_protein_split_hmm_search_agree(toy):
    """On run_dbcan's Z (the profile count) an HMM search gives the same calls however
    the proteins are split, unlike Aksha's raw Z (the control above)."""
    faa, dbcan, tmp = toy
    hmm = dbcan / "dbCAN.hmm"
    prot_dir, _, _ = cc.resolve_proteins(None, str(faa), verbose=False)
    whole = cc.run_astra_hmm(prot_dir, "whole", hmm_path=hmm, threads=1, verbose=False,
                             ievalue_scale=cc.dbcan_ievalue_scale(hmm, cc.count_fasta_records(prot_dir)))
    lines = faa.read_text().splitlines(keepends=True)
    split: dict = {}
    for i in range(3):
        d = tmp / f"zchunk{i}"
        d.mkdir()
        (d / "p.faa").write_text("".join(lines[2 * j] + lines[2 * j + 1] for j in range(i, len(lines) // 2, 3)))
        part = cc.run_astra_hmm(str(d), "chunk", hmm_path=hmm, threads=1, verbose=False,
                                ievalue_scale=cc.dbcan_ievalue_scale(hmm, cc.count_fasta_records(str(d))))
        split.update(part)
    assert set(split) == set(whole)
    for pid, fams in whole.items():
        assert set(split[pid]) == set(fams)
        for fam, (ev, sc) in fams.items():
            assert split[pid][fam][0] == pytest.approx(ev, rel=0.02) and split[pid][fam][1] == sc


def test_dbcan_ievalue_scale_is_profiles_over_targets(tmp_path):
    hmm = tmp_path / "x.hmm"
    _hmm_file(hmm, [(f"GH{i}.hmm", 50) for i in range(8)])
    assert cc.dbcan_ievalue_scale(hmm, 2000) == pytest.approx(8 / 2000)
    assert cc.dbcan_ievalue_scale(hmm, 4, n_profiles=8) == 2.0
    with pytest.raises(ValueError):
        cc.dbcan_ievalue_scale(hmm, 0)
