"""
Jobs und Worker (Plan P0.8c): JSON-Lines-Protokoll, Staging und Commit,
Sammel-CSV nach dem Commit, Abbruch, gekillter Worker, Fehlerklassen,
Worker = API.  Alle Laeufe auf dem Demo-Fall.
"""

from __future__ import annotations

import io
import json
import os
import queue
import subprocess
import sys
import threading
import time
from pathlib import Path

import pydicom
import pytest

from dicom_file_modifier import dose_indices as di
from dicom_file_modifier.api import dose, jobs, selection, structures, transform, worker
from dicom_file_modifier.api.outputs import OutputSpec
from dicom_file_modifier.demo import DemoSpec, make_demo_case

EVENT_TYPES = {"hello", "stage", "progress", "log", "issue", "artifact", "heartbeat", "result"}


@pytest.fixture(scope="module")
def demo(tmp_path_factory):
    return make_demo_case(tmp_path_factory.mktemp("demo") / "case", DemoSpec())


def _job(tmp_path: Path, workflow: str, settings, sel, **out) -> tuple:
    job = jobs.new_job(workflow, settings, sel, OutputSpec(str(tmp_path / "results"), **out))
    return job, job.save(tmp_path / "jobs" / f"{job.job_id}.json")


def _worker_cmd(job_path) -> list:
    return [sys.executable, "-m", "dicom_file_modifier", "worker", "--job", str(job_path)]


def _run_worker(job_path) -> tuple:
    proc = subprocess.run(_worker_cmd(job_path), capture_output=True, text=True, encoding="utf-8",
                          timeout=600, stdin=subprocess.DEVNULL)
    return proc.returncode, [json.loads(line) for line in proc.stdout.splitlines()], proc.stderr


class _Live:
    """
    Worker als Popen; Ereignisse ueber einen Lese-Thread (kein Haengen bei
    readline).  stdin bleibt eine offene Pipe wie bei QProcess: frueher hing der
    Worker damit schon beim numpy-Import (Lese-Thread auf stdin, siehe ``worker``).
    """

    def __init__(self, job_path):
        self.proc = subprocess.Popen(_worker_cmd(job_path), stdin=subprocess.PIPE,
                                     stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                                     text=True, encoding="utf-8")
        self.events: queue.Queue = queue.Queue()
        threading.Thread(target=self._read, daemon=True).start()
        threading.Thread(target=self.proc.stderr.read, daemon=True).start()

    def _read(self):
        for line in self.proc.stdout:
            self.events.put(json.loads(line))
        self.events.put(None)

    def wait_for(self, pred, timeout=300):
        end = time.monotonic() + timeout
        while time.monotonic() < end:
            ev = self.events.get(timeout=end - time.monotonic())
            if ev is None:
                raise AssertionError("Worker beendet, ohne dass das Ereignis kam")
            if pred(ev):
                return ev
        raise AssertionError("Zeitueberschreitung")


def _slow_dose_job(tmp_path, demo):
    """Rechnet einige Sekunden (0.1-mm-Raster), damit Abbruch und Kill mitten hinein treffen."""
    return _job(tmp_path, "dose", dose.Settings(grid_mm=0.1, viz=False, write_rs=False),
                selection.for_dose(str(demo.root)))


def _report(path: Path) -> dict:
    r = json.loads(path.read_text(encoding="utf-8"))
    r["meta"].pop("timestamp")
    r.pop("outputs")
    return r


# -- Protokoll, Staging, Commit ----------------------------------------------------

def test_worker_protocol_commit_and_collection_csv(demo, tmp_path):
    csv_path = tmp_path / "sammel.csv"
    job, path = _job(tmp_path, "dose", dose.Settings(viz=False, append_csv=str(csv_path)),
                     selection.for_dose(str(demo.root)))
    code, events, stderr = _run_worker(path)
    assert code == 0, stderr
    assert {e["type"] for e in events} <= EVENT_TYPES
    hello, result = events[0], events[-1]
    assert hello["type"] == "hello" and hello["v"] == worker.PROTOCOL_VERSION and hello["pid"] > 0
    assert hello["versions"]["dfm"] and all(hello["threads"].values())
    assert result["type"] == "result" and result["status"] == "ok"      # Eclipse-Abgleich ohne Fehlalarm
    stages = [(e["key"], e["i"], e["n"]) for e in events if e["type"] == "stage"]
    assert stages == [("prepare", 1, 4), ("compute", 2, 4), ("rs_export", 3, 4), ("reports", 4, 4)]
    assert any(e["type"] == "log" and "Dosisindex-Berechnung" in e["msg"] for e in events)

    res = result["result"]
    final = Path(res["output_dir"])
    assert final == tmp_path / "results" / f"{demo.root.name}_IDX"
    assert not job.staging_dir.exists() and sorted(p.name for p in final.parent.iterdir()) == [final.name]
    artifacts = {e["role"]: Path(e["path"]) for e in events if e["type"] == "artifact"}
    assert set(artifacts) == {"rs", "json", "txt", "csv", "append_csv"} and all(p.exists() for p in artifacts.values())
    assert {m["path"] for m in res["manifest"]} == {p.name for p in final.iterdir()} - {"run.json", "run.log"}

    run = json.loads((final / "run.json").read_text(encoding="utf-8"))
    assert run["job_id"] == job.job_id and run["result"]["status"] == "ok"
    assert run["command"] == res["command"] and "--append-csv" in run["command_text"]
    assert {i["role"] for i in run["inputs"]} == {"rs", "rd", "rp", "ct"}
    assert all(i.get("sop_instance_uid") for i in run["inputs"] if i["role"] != "ct")
    assert "Dosisindex-Berechnung" in (final / "run.log").read_text(encoding="utf-8")

    # die Sammel-CSV entsteht erst nach dem Commit, mit denselben Zeilen wie indices.csv
    assert csv_path.read_text(encoding="utf-8") == (final / "indices.csv").read_text(encoding="utf-8")

    # Worker = API (gleiche Einstellungen ohne Sammel-CSV)
    api = dose.run(dose.Settings(viz=False), selection.for_dose(str(demo.root)), tmp_path / "api")
    assert api.ok
    cid = demo.root.name
    assert _report(final / f"{cid}_indices.json") == _report(tmp_path / "api" / f"{cid}_indices.json")


def test_second_run_gets_a_new_folder_or_replaces_it(demo, tmp_path):
    sel = selection.from_ct_dir(str(demo.ct_dir))
    s = transform.Settings(tx=1.0, method="metadata", viz=False)
    first = jobs.run_job_inprocess(_job(tmp_path, "transform", s, sel)[0])
    second = jobs.run_job_inprocess(_job(tmp_path, "transform", s, sel)[0])
    assert Path(second.output_dir).name == Path(first.output_dir).name + "_2"
    marker = Path(first.output_dir) / "alt.txt"
    marker.write_text("alt", encoding="utf-8")
    third = jobs.run_job_inprocess(_job(tmp_path, "transform", s, sel, policy="overwrite")[0])
    assert third.ok and third.output_dir == first.output_dir and not marker.exists()
    assert sorted(p.name for p in (tmp_path / "results").iterdir()) == [
        Path(first.output_dir).name, Path(second.output_dir).name]


def test_inprocess_jobs_equal_the_api(demo, tmp_path, capsys):
    s = structures.Settings()
    job, _ = _job(tmp_path, "structures", s, selection.from_rtstruct(str(demo.rs)))
    res = jobs.run_job_inprocess(job)
    assert res.ok and capsys.readouterr().out == ""               # print landet in run.log
    assert structures.run(s, selection.from_rtstruct(str(demo.rs)), tmp_path / "api").ok
    capsys.readouterr()
    final = Path(res.output_dir)
    got = {p.name: p.read_bytes() for p in final.iterdir() if p.name not in ("run.json", "run.log")}
    assert got == {p.name: p.read_bytes() for p in (tmp_path / "api").iterdir()}
    assert "RTSTRUCT Analyse" in (final / "run.log").read_text(encoding="utf-8")


# -- Abbruch, Kill, Aufraeumen ---------------------------------------------------------

def test_cancel_file_leaves_nothing(demo, tmp_path):
    job, path = _slow_dose_job(tmp_path, demo)
    live = _Live(path)
    live.wait_for(lambda e: e["type"] == "stage" and e["key"] == "compute")
    assert job.staging_dir.exists()
    worker.cancel_file(path).touch()
    result = live.wait_for(lambda e: e["type"] == "result")
    assert live.proc.wait(timeout=60) == 3
    assert result["status"] == "cancelled" and result["result"]["issues"][-1]["code"] == "JOB.CANCELLED"
    assert not job.staging_dir.exists() and not any((tmp_path / "results").iterdir())


def test_killed_worker_leaves_nothing_after_cleanup(demo, tmp_path):
    job, path = _slow_dose_job(tmp_path, demo)
    live = _Live(path)
    live.wait_for(lambda e: e["type"] == "stage" and e["key"] == "compute")
    live.proc.kill()
    live.proc.wait(timeout=60)
    assert job.staging_dir.exists()                  # hart beendet: der Worker raeumt nicht mehr auf
    assert jobs.remove_tree(job.staging_dir)         # das macht die App nach dem Kill
    assert not any((tmp_path / "results").iterdir())


def test_cleanup_staging_removes_only_job_folders(tmp_path):
    for name in (".stgabc123", ".old12345678_case_IDX", "case_IDX", ".hidden"):
        (tmp_path / name).mkdir()
    (tmp_path / ".stgabc123" / "x.dcm").write_bytes(b"1")
    os.chmod(tmp_path / ".stgabc123" / "x.dcm", 0o444)             # schreibgeschuetzt
    removed = jobs.cleanup_staging(tmp_path)
    assert sorted(p.name for p in removed) == [".old12345678_case_IDX", ".stgabc123"]
    assert sorted(p.name for p in tmp_path.iterdir()) == [".hidden", "case_IDX"]


# -- Fehlerklassen ---------------------------------------------------------------------

def test_input_error_is_exit_2_without_writing(demo, tmp_path):
    sel = selection.for_dose(str(demo.root))
    sel.rs = str(tmp_path / "fehlt.dcm")
    _, path = _job(tmp_path, "dose", dose.Settings(viz=False), sel)
    code, events, _ = _run_worker(path)
    assert code == 2 and events[-1]["status"] == "failed"
    assert events[-1]["result"]["issues"][-1]["code"] == "FILE.NOT_FOUND"
    assert not (tmp_path / "results").exists() or not any((tmp_path / "results").iterdir())


def test_bad_collection_csv_blocks_before_the_run(demo, tmp_path):
    csv_path = tmp_path / "alt.csv"
    csv_path.write_text("case_id,target\nx,y\n", encoding="utf-8")
    job, _ = _job(tmp_path, "dose", dose.Settings(viz=False, append_csv=str(csv_path)),
                  selection.for_dose(str(demo.root)))
    res = jobs.run_job_inprocess(job)
    assert (res.status, res.exit_code) == ("failed", 2)
    assert [(i.code, i.field) for i in res.issues if i.level == "error"] == [("INPUT.INVALID", "append_csv")]
    assert not job.staging_dir.exists() and csv_path.read_text(encoding="utf-8") == "case_id,target\nx,y\n"


def test_missing_or_broken_job_file_is_exit_2(tmp_path):
    code, events, _ = _run_worker(tmp_path / "gibtsnicht.json")
    assert code == 2 and [e["type"] for e in events][::len(events) - 1] == ["hello", "result"]
    bad = tmp_path / "bad.json"
    bad.write_text(json.dumps({"v": 99}), encoding="utf-8")
    code, events, _ = _run_worker(bad)
    assert code == 2 and "Job-Version" in events[-1]["result"]["issues"][-1]["message_de"]


def test_internal_errors_are_exit_1(demo, tmp_path, monkeypatch):
    job, path = _job(tmp_path, "dose", dose.Settings(viz=False, write_rs=False),
                     selection.for_dose(str(demo.root)))
    monkeypatch.setattr(di, "compute_dose_indices", lambda *a, **k: 1 / 0)
    proto = io.StringIO()
    assert worker.run_worker(str(path), proto) == 1
    last = json.loads(proto.getvalue().splitlines()[-1])
    assert last["status"] == "failed" and last["result"]["issues"][-1]["code"] == "INTERNAL"
    assert "ZeroDivisionError" in last["result"]["issues"][-1]["detail"]
    assert not job.staging_dir.exists()

    monkeypatch.setattr(jobs, "execute_job", lambda *a, **k: 1 / 0)    # Fehler im Job-Rahmen
    proto = io.StringIO()
    assert worker.run_worker(str(path), proto) == 1
    assert json.loads(proto.getvalue().splitlines()[-1])["result"]["issues"][-1]["code"] == "INTERNAL"


def test_protocol_survives_writes_to_fd1(tmp_path):
    code = ("import os, sys\n"
            "from dicom_file_modifier.api import worker\n"
            "proto = worker.protect_stdout()\n"
            "os.write(1, b'c-bibliothek\\n')\n"
            "print('print-ausgabe', flush=True)\n"
            "proto.write('{\"type\": \"hello\"}\\n'); proto.flush()\n")
    out = subprocess.run([sys.executable, "-c", code], capture_output=True, text=True, timeout=120)
    assert out.returncode == 0, out.stderr
    assert out.stdout == '{"type": "hello"}\n'
    assert "c-bibliothek" in out.stderr and "print-ausgabe" in out.stderr


def test_transform_job_equals_the_api(demo, tmp_path, capsys):
    s = transform.Settings(tx=3.0, rz=5.0, method="metadata", viz=False)
    job, _ = _job(tmp_path, "transform", s, selection.for_transform(str(demo.root)))
    res = jobs.run_job_inprocess(job)
    api = transform.run(s, selection.for_transform(str(demo.root)), tmp_path / "api")
    capsys.readouterr()
    assert res.ok and api.ok and res.outputs == api.outputs

    def geometry(d: Path) -> list:
        return [(p.name, list(pydicom.dcmread(str(p), stop_before_pixels=True).ImagePositionPatient))
                for p in sorted((d / "CT").glob("*.dcm"))]

    assert geometry(Path(res.output_dir)) == geometry(tmp_path / "api")
