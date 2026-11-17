LBNH Process Twin - AI-built Python simulator (Ripa Sutar, INEN 6301, Lamar University)
=====================================================================================
LNG plant with boil-off-gas helium recovery and liquefaction (my Aspen Plus design,
Sutar, Modi & Xu, I&EC Research), rebuilt as a Python simulator by Claude (Anthropic).

SETUP (once)
  1. Install Python 3.10 or newer (python.org). On Windows, tick "Add Python to PATH".
  2. Double-click 1_INSTALL.bat   (or: python -m pip install -r requirements.txt)
     Needs internet, about 2 minutes. Packages: numpy, scipy, CoolProp, thermo, chemicals.

RUN (Windows: double-click; or type the command in a terminal opened in this folder)
  5_LIVE_SIMULATOR.bat   python live_app.py
        Opens the simulator in the browser (http://localhost:8765). Keep the black window open.
        Tab 1  Live simulator: 3D plant; change feed flow / helium / nitrogen and press Solve -
               Python re-solves the whole plant (30-40 s). Prices are editable.
        Tab 2  My Aspen model: Aspen flowsheet next to the flowsheet generated from Python,
               stream-by-stream comparison, Aspen vs AI results, and how the simulator was
               built (with the checks you can run).
        Tab 3  Python code, all files.
  3_SAME_DESIGN.bat      python same_design.py
        Every Aspen compressor and pump recomputed by Python from the same inlet state,
        outlet pressure and efficiencies (aspen_settings.json): 331.4 vs 332.5 MW.
  4_CHANGE_FEED.bat      python change_feed.py 660 0.137 7.63
        Solve the whole plant for another feed (flow t/h, helium wt %, nitrogen wt %).
  2_RUN_DEMO.bat         python demo.py
        The AI simulator's OWN design of the plant (its own pressures and efficiencies,
        from the blind build) printed next to my Aspen results.
  python compare.py                    full AI vs Aspen comparison table
  python flowsheet/stream_report.py    re-solve my design and record every stream (~30 s)
  python flowsheet/draw_pfd.py         redraw the flowsheet PDF
  python flowsheet/make_excel.py       rebuild the Excel stream table

HOW IT WAS BUILT (see also Tab 2)
  1. Blind build: only the paper's Table 1 feed, block diagram, technologies and property
     method were given. Every other number is in design_basis.py, labelled GIVEN / CHOSEN /
     TYPICAL. Predictions were sealed before any Aspen result was seen:
     predictions_sealed.json, SHA-256 in predictions_sealed.sha256.
     Check:  certutil -hashfile predictions_sealed.json SHA256   (Windows)
  2. Same-design check: Aspen settings given as INPUTS only (aspen_settings.json);
     Python computes the power itself. Change a value there and re-run same_design.py:
     the power changes.
  3. What-if studies: full re-solve in Python. Tab 1 shows the result on my Aspen basis
     (Aspen design-feed numbers x the % change computed by Python).
  Aspen results appear only in live_app.py (ASPEN block), demo.py/compare.py (comparison
  columns) and the settings files (aspen_settings.json, section_models/aspen/streams.json) -
  not in the thermodynamics, flash, column, exchanger, refrigeration-cycle or liquefier code.

MAIN FILES
  eos.py, props.py      Peng-Robinson thermodynamics (helium fitted to GERG-2008)
  flash.py              phase equilibrium (Rachford-Rice, stability test)
  column.py, nru.py     rigorous N2 rejection column (Naphtali-Sandholm)
  lngp.py               C3MR liquefaction, main exchanger composite curves, compressors
  coldbox.py            N2 and He refrigeration cycles of the NRHR cold box
  hls.py                Collins helium liquefier
  plant.py              the whole complex, boil-off recycle included (Wegstein)
  section_models/       LNG section on my Aspen settings (econ_opt.py) and helium-part models
                        (he_opt.py), used by the flowsheet stream report. These files also
                        contain exploratory optimization routines that are NOT used or reported:
                        a trustworthy optimization needs costs from my Aspen economic analysis.
  flowsheet/            stream report, flowsheet PDF, Excel stream table, Aspen comparison
  live_app.py, live_page.html, plant3d.js   the browser simulator (works offline)

Limits: steady state; no equipment rating (feed changes assume design margin); economics
taken from my paper; storage and helium parts laid out differently from Aspen.
