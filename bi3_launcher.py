"""Start the Bi3+ Luminescence-Structure Calculator in the default web browser.

Usage:
    python bi3_launcher.py              start the calculator and open the browser
    python bi3_launcher.py --check      run a short self-test and exit
    python bi3_launcher.py --port 8600  use a fixed port
    python bi3_launcher.py --no-browser do not open the browser automatically

After "pip install ." the same program is available as the command bi3-calculator.
"""
import argparse
import os
import socket
import subprocess
import sys
import time
import webbrowser
from pathlib import Path

HERE = Path(__file__).resolve().parent
APP = HERE / "bi3_app.py"


def find_free_port(start=8501, tries=60):
    """Return the first port at or above `start` on which nothing is listening."""
    for port in range(start, start + tries):
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(0.5)
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise RuntimeError("No free port found between %d and %d" % (start, start + tries - 1))


def wait_for_server(port, timeout=90.0):
    """True once the local server accepts connections."""
    t_end = time.time() + timeout
    while time.time() < t_end:
        with socket.socket(socket.AF_INET, socket.SOCK_STREAM) as s:
            s.settimeout(1.0)
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return True
        time.sleep(0.5)
    return False


def self_check():
    """Import all libraries and recompute one published benchmark (YPO4, Boutinaud Tables 5 and 6)."""
    import numpy, pandas, matplotlib, streamlit  # noqa: F401
    sys.path.insert(0, str(HERE))
    import bi3_engine as eng

    bench = eng.get_boutinaud_2023_benchmarks()
    key = [k for k in bench if k.startswith("YPO4")][0]
    res = eng.run_boutinaud_2023_calculation(bench[key])
    row = res["summary_df"].iloc[0]
    he_m, e_a = float(row["he(M) [Calc]"]), float(row["E_A (eV) [Eq. 18]"])
    ok = abs(he_m - 1.12) < 0.01 and abs(e_a - 3.93) < 0.01
    print("Python %s | streamlit %s | numpy %s | pandas %s | matplotlib %s" % (
        sys.version.split()[0], streamlit.__version__, numpy.__version__,
        pandas.__version__, matplotlib.__version__))
    print("YPO4 Y[8]: he(M) = %.3f (published 1.12), E_A = %.2f eV (published 3.93)" % (he_m, e_a))
    print("Self-test passed." if ok else "Self-test FAILED: results differ from the published values.")
    return 0 if ok else 1


def main(argv=None):
    ap = argparse.ArgumentParser(description="Bi3+ Luminescence-Structure Calculator")
    ap.add_argument("--port", type=int, default=None, help="fixed port (default: first free port from 8501)")
    ap.add_argument("--no-browser", action="store_true", help="do not open the web browser")
    ap.add_argument("--check", action="store_true", help="run the self-test and exit")
    args = ap.parse_args(argv)

    if args.check:
        return self_check()

    port = args.port or find_free_port()
    env = dict(os.environ)
    env["STREAMLIT_BROWSER_GATHER_USAGE_STATS"] = "false"
    cmd = [sys.executable, "-m", "streamlit", "run", str(APP),
           "--server.port", str(port), "--server.headless", "true",
           "--server.address", "localhost",
           "--browser.gatherUsageStats", "false"]
    url = "http://localhost:%d" % port
    print("Starting the calculator at %s  (close this window or press Ctrl+C to stop)" % url)
    proc = subprocess.Popen(cmd, cwd=str(HERE), env=env)
    try:
        if wait_for_server(port):
            if not args.no_browser:
                webbrowser.open(url)
        else:
            print("The server did not start within 90 s. See the messages above.")
        return proc.wait()
    except KeyboardInterrupt:
        proc.terminate()
        return 0


if __name__ == "__main__":
    sys.exit(main())
