import importlib.util
import multiprocessing as mp
from pathlib import Path
from types import SimpleNamespace

import pytest


@pytest.fixture
def web_app(monkeypatch):
    pytest.importorskip("flask")
    path = Path(__file__).resolve().parents[2] / "web" / "app.py"
    spec = importlib.util.spec_from_file_location("schubmult_web_test", path)
    module = importlib.util.module_from_spec(spec)
    with monkeypatch.context() as setup:
        setup.setenv("SCHUBMULT_ACCESS_LOG", "-")
        setup.delenv("SCHUBMULT_MAX_DOWNLOAD_BYTES", raising=False)
        setup.setattr(mp, "current_process", lambda: SimpleNamespace(name="TestProcess"))
        spec.loader.exec_module(module)
    return module


def test_default_download_limit(web_app):
    assert web_app.MAX_DOWNLOAD_BYTES == 100_000_000
    buffer = web_app._CappedBuffer(web_app.MAX_DOWNLOAD_BYTES, stop_on_cap=True, count_bytes=True)
    chunk = "x" * 1_000_000
    for _ in range(100):
        assert buffer.write(chunk) == len(chunk)
    assert buffer.tell() == 100_000_000
    assert not buffer._truncated
    with pytest.raises(BrokenPipeError):
        buffer.write("x")
    assert buffer._truncated


def test_utf8_limit_preserves_complete_characters(web_app):
    buffer = web_app._CappedBuffer(5, count_bytes=True)
    buffer.write("\u03b1")
    buffer.write("\u03b2\u03b3")
    content = buffer.getvalue()
    assert content.startswith("\u03b1\u03b2\n")
    assert "truncated at 5 bytes" in content
    assert "\ufffd" not in content
    buffer.write("ignored")
    assert buffer.getvalue() == content


@pytest.mark.parametrize("download,limit", [(False, 10), (True, 30)])
def test_inline_selects_limit_and_keeps_stderr_small(web_app, monkeypatch, download, limit):
    monkeypatch.setattr(web_app, "MAX_OUTPUT_BYTES", 10)
    monkeypatch.setattr(web_app, "MAX_DOWNLOAD_BYTES", 30)

    def main(argv):
        import sys

        sys.stderr.write("e" * 11)
        sys.stdout.write("x" * 31)

    monkeypatch.setattr(web_app, "_script_module", lambda flavor: SimpleNamespace(main=main))
    out, err, timed_out, _ = web_app._run_inline("py", ["schubmult_py"], download=download)
    assert out.startswith("x" * limit + "\n")
    assert f"truncated at {limit}" in out
    assert err.startswith("e" * 10 + "\n")
    assert "truncated at 10 characters" in err
    assert not timed_out


@pytest.mark.parametrize("download", [False, True])
@pytest.mark.parametrize("inline", [False, True])
def test_runner_passes_download_choice(web_app, monkeypatch, download, inline):
    result = ("output", "", False, 0.1)
    if inline:
        monkeypatch.setenv("SCHUBMULT_DISABLE_SUBPROCESS", "1")

        def run(flavor, argv, *, download: bool):
            assert download is expected
            return result

        expected = download
        monkeypatch.setattr(web_app, "_run_inline", run)
    else:
        monkeypatch.delenv("SCHUBMULT_DISABLE_SUBPROCESS", raising=False)
        messages = []
        worker = SimpleNamespace(conn=SimpleNamespace(
            send=messages.append, poll=lambda timeout: True,
            recv=lambda: ("output", "", 0.1),
        ))
        monkeypatch.setattr(web_app, "_acquire_worker", lambda: worker)
        monkeypatch.setattr(web_app, "_release_worker", lambda worker: None)
    assert web_app._run_script("py", ["schubmult_py"], download=download) == result
    if not inline:
        assert messages == [("py", ["schubmult_py"], download)]


@pytest.mark.parametrize("download", [None, False, True])
def test_api_passes_download_choice(web_app, monkeypatch, download):
    calls = []

    def run(flavor, argv, *, download):
        calls.append(download)
        return ("1  (4, 1, 2, 3)\n", "", False, 0.1)

    monkeypatch.setattr(web_app, "_run_script", run)
    payload = {"flavor": "py", "perms": "3 1 2 - 2 1 3"}
    if download is not None:
        payload["download"] = download
    response = web_app.app.test_client().post("/api/compute", json=payload)
    assert response.status_code == 200
    assert response.json["stdout"] == "1  (4, 1, 2, 3)\n"
    assert calls == [download is True]


def test_api_rejects_non_boolean_download(web_app, monkeypatch):
    def run(*args, **kwargs):
        pytest.fail("Invalid request must not compute")

    monkeypatch.setattr(web_app, "_run_script", run)
    response = web_app.app.test_client().post("/api/compute", json={"download": "false"})
    assert response.status_code == 400
    assert response.json == {"ok": False, "error": "download must be a boolean"}


def test_worker_passes_download_choice(web_app, monkeypatch):
    messages = iter([("py", ["schubmult_py"], True)])
    replies = []

    def recv():
        try:
            return next(messages)
        except StopIteration:
            raise EOFError from None

    def run(flavor, argv, *, download):
        assert download is True
        return ("output", "", False, 0.1)

    monkeypatch.setattr(web_app, "_run_inline", run)
    web_app._serve(SimpleNamespace(recv=recv, send=replies.append))
    assert replies == [("output", "", 0.1)]


def test_real_computation_matches_display_and_download(web_app, monkeypatch):
    monkeypatch.setattr(web_app, "COMPUTE_TIMEOUT", 30)
    monkeypatch.setattr(web_app, "MP_START_METHOD", "fork")
    monkeypatch.delenv("SCHUBMULT_DISABLE_SUBPROCESS", raising=False)
    argv = ["schubmult_py", "3", "1", "2", "-", "2", "1", "3"]
    display = web_app._run_inline("py", argv)
    download = web_app._run_inline("py", argv, download=True)
    assert display[0]
    assert display[:3] == download[:3]
    assert display[1:3] == ("", False)
    if "fork" not in mp.get_all_start_methods():
        return
    try:
        worker_download = web_app._run_script("py", argv, download=True)
        assert worker_download[:3] == display[:3]
    finally:
        for worker in web_app._idle_workers:
            worker.kill()
        web_app._idle_workers.clear()
