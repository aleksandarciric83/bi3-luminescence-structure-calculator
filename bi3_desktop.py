"""Bi3+ Luminescence-Structure Calculator as a desktop program (own window, no browser).

    python bi3_desktop.py            open the calculator in its own window
    Bi3_Calculator.exe               the same, after build_exe.bat has packaged it
    Bi3_Calculator.exe --browser     use the web browser instead of the own window (this happens
                                     automatically when the Microsoft Edge WebView2 component is missing,
                                     and on Linux; `--window` tries an own window on Linux)
    On macOS the packaged program is Bi3_Calculator.app (build_app.py), on Linux the file Bi3_Calculator.

The calculation and the user interface are the unchanged bi3_app.py; this file only shows the page
in a native window (pywebview) instead of the web browser. The local server listens on 127.0.0.1 only.
"""
import os
import socket
import subprocess
import sys
import time
from pathlib import Path

# Listed so that PyInstaller bundles them; bi3_app.py imports them at run time.
import bi3_engine, bi3_export, bi3_format, bi3_plots, bi3_validate  # noqa: F401,E401

BASE = Path(getattr(sys, "_MEIPASS", Path(__file__).resolve().parent))
APP = BASE / "bi3_app.py"
TITLE = "Bi3+ Luminescence-Structure Calculator"
ICON = BASE / "bi3_icon.ico"
WEBVIEW2_ID = "{F3017226-FE2A-4295-8BDF-00C3A9A7E4C5}"


def free_port(start=8601):
    for port in range(start, start + 100):
        with socket.socket() as s:
            s.settimeout(0.3)
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise RuntimeError("no free port")


def serve(port):
    """Run the Streamlit server in this process (called in the child process)."""
    for name in ("stdout", "stderr"):  # no console in the packaged program; keep Unicode output safe
        stream = getattr(sys, name)
        if stream is None:
            setattr(sys, name, open(os.devnull, "w", encoding="utf-8"))
        else:
            stream.reconfigure(encoding="utf-8", errors="replace")
    sys.path.insert(0, str(BASE))
    os.chdir(str(BASE))
    from streamlit.web import bootstrap
    flags = {"server.port": port, "server.address": "127.0.0.1", "server.headless": True,
             "browser.gatherUsageStats": False, "global.developmentMode": False,
             "server.fileWatcherType": "none", "client.toolbarMode": "minimal"}
    bootstrap.load_config_options(flags)
    bootstrap.run(str(APP), False, [], flags)


def wait_for_port(port, timeout=90.0):
    t_end = time.time() + timeout
    while time.time() < t_end:
        with socket.socket() as s:
            s.settimeout(1.0)
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return True
        time.sleep(0.4)
    return False


def webview2_available():
    """True when the Microsoft Edge WebView2 runtime (part of Windows 11 and of updated Windows 10) is installed."""
    if sys.platform != "win32":
        return True
    import winreg
    clients = "Microsoft\\EdgeUpdate\\Clients\\" + WEBVIEW2_ID
    for hive, path in (
        (winreg.HKEY_LOCAL_MACHINE, "SOFTWARE\\WOW6432Node\\" + clients),
        (winreg.HKEY_LOCAL_MACHINE, "SOFTWARE\\" + clients),
        (winreg.HKEY_CURRENT_USER, "Software\\" + clients),
    ):
        try:
            with winreg.OpenKey(hive, path) as key:
                version = winreg.QueryValueEx(key, "pv")[0]
            if version and version != "0.0.0.0":
                return True
        except OSError:
            continue
    return False


def use_browser_mode(argv=None, platform=None):
    """True when the calculator is shown in the web browser instead of an own window.

    Windows: own window (WebView2) when the component exists. macOS: own window (WKWebView is part of macOS).
    Linux: web browser, because there is no web-view component that every distribution has;
    `--window` tries an own window anyway (needs pywebview with GTK or Qt).
    """
    argv = sys.argv if argv is None else argv
    platform = sys.platform if platform is None else platform
    if "--browser" in argv:
        return True
    if "--window" in argv:
        return False
    if platform == "win32":
        return not webview2_available()
    return platform != "darwin"


def run_in_browser(port, server=None):
    """Open the default web browser and keep a small window with a Stop button (or the terminal, without a display)."""
    import webbrowser

    url = "http://127.0.0.1:%d" % port
    try:
        import tkinter as tk
        root = tk.Tk()
    except Exception:  # no tkinter or no display (for example a server or a minimal Linux)
        print("Calculator: %s   (press Ctrl+C to stop)" % url, flush=True)
        if wait_for_port(port):
            webbrowser.open(url)
        try:
            while server is None or server.poll() is None:
                time.sleep(0.5)
        except KeyboardInterrupt:
            pass
        return
    root.title(TITLE)
    try:
        root.iconbitmap(str(ICON))
    except Exception:
        pass
    root.geometry("460x150")
    msg = tk.Label(root, text="Starting the calculator ...", font=("Segoe UI", 11), wraplength=420, justify="center")
    msg.pack(pady=(24, 10))
    tk.Button(root, text="Stop the calculator", command=root.destroy, width=22).pack()
    state = {"opened": False, "t0": time.time()}

    def poll():
        with socket.socket() as s:
            s.settimeout(0.3)
            up = s.connect_ex(("127.0.0.1", port)) == 0
        if up:
            webbrowser.open(url)
            msg.config(text="The calculator is open in your web browser:\n%s\nClose this window to stop it." % url)
        elif time.time() - state["t0"] > 90:
            msg.config(text="The calculator did not start.")
        else:
            root.after(400, poll)

    root.after(300, poll)
    root.mainloop()


def selftest(outfile=None):
    """Calculate one built-in structure and write every export file (PNG, SVG, PDF, ZIP); proves that a packaged
    program has all libraries it needs. Result goes to `outfile` (a packaged Windows program has no console)."""
    try:
        c = bi3_engine.get_boutinaud_2023_benchmarks()["CaYAlO4 (ICDD 01-081-0742, I4/mmm)"]
        res = bi3_engine.run_boutinaud_2023_calculation(c, custom_exp_ex_nm=280.0, custom_exp_em_nm=398.0)
        sp = res["summary_df"].drop(columns=["Raw Site"])
        files = bi3_export.figure_files(sp, c.formula, 280.0, 398.0)
        assert files["energy_levels.png"][:4] == b"\x89PNG" and files["EA_vs_he.pdf"][:4] == b"%PDF"
        assert b"<svg" in files["energy_levels.svg"][:2000] and all(len(v) > 1000 for v in files.values())
        zipped = bi3_export.build_zip(sp, res["matrix_I_df"], res["subshells_df"], c.formula, c.space_group, "V = 1", "",
                                      280.0, 398.0, "selftest", "selftest")
        assert len(zipped) > 10000
        import streamlit  # noqa: F401
        message, code = "Self-test passed: %d figure files, ZIP %d bytes, streamlit %s" % (
            len(files), len(zipped), streamlit.__version__), 0
    except Exception as e:  # report every kind of failure
        message, code = "Self-test FAILED: %r" % (e,), 1
    if outfile:
        Path(outfile).write_text(message + "\n", encoding="utf-8")
    else:
        print(message)
    return code


def main():
    if len(sys.argv) == 3 and sys.argv[1] == "--serve":
        serve(int(sys.argv[2]))
        return 0
    if "--selftest" in sys.argv:
        i = sys.argv.index("--selftest")
        return selftest(sys.argv[i + 1] if i + 1 < len(sys.argv) else None)

    if sys.platform != "win32":  # stop the server too when this program is stopped with `kill`
        import signal
        signal.signal(signal.SIGTERM, lambda *_: sys.exit(0))
    use_browser = use_browser_mode()
    port = free_port()
    exe = [sys.executable] if getattr(sys, "frozen", False) else [sys.executable, str(Path(__file__).resolve())]
    flags = subprocess.CREATE_NO_WINDOW if sys.platform == "win32" else 0
    server = subprocess.Popen(exe + ["--serve", str(port)], creationflags=flags)
    url = "http://127.0.0.1:%d" % port

    def load(window):
        if wait_for_port(port):
            window.load_url(url)
        else:
            window.load_html("<h3>The calculator did not start.</h3>")

    try:
        if not use_browser:
            try:
                import webview
                webview.settings["ALLOW_DOWNLOADS"] = True
                window = webview.create_window(
                    TITLE, html="<body style='font-family:sans-serif;margin:3em'><h2>Starting the calculator ...</h2></body>",
                    width=1500, height=950, min_size=(1000, 650))
                webview.start(load, window)
            except Exception as e:  # no usable web-view component: use the browser instead
                print("Own window not available (%s); using the web browser." % e, flush=True)
                use_browser = True
        if use_browser:
            run_in_browser(port, server)
    finally:
        if sys.platform == "win32":
            subprocess.run(["taskkill", "/PID", str(server.pid), "/T", "/F"], capture_output=True)
        else:
            server.terminate()
    return 0


if __name__ == "__main__":
    sys.exit(main())
