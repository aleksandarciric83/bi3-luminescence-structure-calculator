"""Take the screenshots for the paper from the running calculator.

Run this on your own computer (a live browser session is needed):

    python -m pip install playwright
    python make_screenshots.py

It starts the calculator on a free local port, opens it with the Google Chrome or Microsoft Edge that is
already installed (no extra download), saves one full-page PNG per tab into the folder "screenshots" and stops
the program again. Options:  --out FOLDER  --width 1500  --height 1000  --scale 2
"""
import argparse
import re
import socket
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent

TABS = [
    ("01_optical_transitions", "Optical Transitions"),
    ("02_matrices_bond_valence", "Matrices I"),
    ("03_dielectric_nephelauxetic", "Dielectric"),
    ("04_energy_diagrams_export", "Energy Diagrams"),
    ("05_benchmark_suite", "Benchmark Suite"),
    ("06_help_and_citation", "Help"),
]


def free_port(start=8765):
    for port in range(start, start + 50):
        with socket.socket() as s:
            s.settimeout(0.5)
            if s.connect_ex(("127.0.0.1", port)) != 0:
                return port
    raise RuntimeError("no free port")


def wait_for_port(port, timeout=90.0):
    t_end = time.time() + timeout
    while time.time() < t_end:
        with socket.socket() as s:
            s.settimeout(1.0)
            if s.connect_ex(("127.0.0.1", port)) == 0:
                return True
        time.sleep(0.5)
    return False


def settle(page, extra_ms=2500):
    """Wait until Streamlit has finished drawing."""
    try:
        page.wait_for_load_state("networkidle", timeout=8000)
    except Exception:
        pass
    try:  # Streamlit shows its "Running" indicator while a rerun is in progress
        page.wait_for_selector('[data-testid="stStatusWidget"]', state="detached", timeout=20000)
    except Exception:
        pass
    page.wait_for_timeout(extra_ms)


def shoot(page, path, min_height):
    """Screenshot of the whole page. Streamlit scrolls an inner container, so full_page alone cuts the page off:
    the viewport is made as tall as the content first."""
    width = page.viewport_size["width"]
    height = page.evaluate("""() => {
        const m = document.querySelector('[data-testid="stMain"]') || document.querySelector('section.main');
        return m ? m.scrollHeight : document.documentElement.scrollHeight;
    }""")
    page.set_viewport_size({"width": width, "height": max(min_height, int(height) + 10)})
    page.wait_for_timeout(1200)
    page.screenshot(path=str(path), full_page=True)
    page.set_viewport_size({"width": width, "height": min_height})


def open_browser(p):
    for channel in ("chrome", "msedge"):
        try:
            return p.chromium.launch(channel=channel, headless=True)
        except Exception:
            continue
    return p.chromium.launch(headless=True)  # works after:  python -m playwright install chromium


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--out", default="screenshots")
    ap.add_argument("--width", type=int, default=1500)
    ap.add_argument("--height", type=int, default=1000)
    ap.add_argument("--scale", type=float, default=2.0, help="device pixel ratio (2 = sharp, for print)")
    args = ap.parse_args()

    try:
        from playwright.sync_api import sync_playwright
    except ImportError:
        sys.exit("Playwright is not installed. Run:  python -m pip install playwright")

    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    port = free_port()
    proc = subprocess.Popen([sys.executable, str(HERE / "bi3_launcher.py"), "--no-browser", "--port", str(port)], cwd=str(HERE))
    try:
        if not wait_for_port(port):
            sys.exit("The calculator did not start. Run Start_Calculator.bat once and try again.")
        with sync_playwright() as p:
            browser = open_browser(p)
            ctx = browser.new_context(viewport={"width": args.width, "height": args.height}, device_scale_factor=args.scale)
            page = ctx.new_page()
            page.goto("http://localhost:%d" % port)
            settle(page, 5000)
            shoot(page, str(out / "00_start_page.png"), args.height)
            for name, label in TABS:
                page.get_by_role("tab", name=re.compile(label)).click()
                settle(page, 4000 if "Diagrams" in label else 2500)
                shoot(page, str(out / (name + ".png")), args.height)
                print("saved", name + ".png")
            page.get_by_role("tab", name=re.compile(TABS[0][1])).click()  # back to the first tab
            page.get_by_text("Upload .json File (MP / Custom)").first.click()
            settle(page, 7000)
            shoot(page, str(out / "07_json_input_cell_parameters.png"), args.height)
            print("saved 07_json_input_cell_parameters.png")
            browser.close()
    finally:
        if sys.platform == "win32":  # also stop the Streamlit child, else it stays running
            subprocess.run(["taskkill", "/PID", str(proc.pid), "/T", "/F"], capture_output=True)
        else:
            proc.terminate()
    print("Done. Screenshots are in", out.resolve())


if __name__ == "__main__":
    main()
