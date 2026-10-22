"""Desktop UI; run python app.py or python app.py --demo."""
from __future__ import annotations

import argparse
import json
import os
from pathlib import Path
import tkinter as tk
from tkinter import filedialog, messagebox, ttk

import numpy as np
import matplotlib
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure

from core import (Settings, analyze, demo, geometry, image_digest, inspect_czi, load_czi, overlay, source_annotation,
                  preprocess, save_result, suggest_threshold)

matplotlib.rcParams["font.sans-serif"] = ["Microsoft YaHei", "SimHei", "DejaVu Sans"]
matplotlib.rcParams["axes.unicode_minus"] = False

PARAMETERS = [
    ("length_um", "Segment length (um)"), ("start_offset_um", "Start offset (um)"),
    ("half_width_um", "ROI half-width (um)"), ("threshold", "Positive threshold"),
    ("shaft_threshold", "Shaft reference threshold (0=linked)"),
    ("gaussian_um", "Gaussian sigma (um; 0=off)"),
    ("background_radius_um", "Background radius (um; 0=off)"),
    ("shaft_halfwidth_um", "Fixed shaft half-width (um; 0=auto)"),
    ("shaft_window_um", "Shaft median window (um)"),
    ("shaft_min_halfwidth_um", "Min. shaft half-width (um)"),
    ("shaft_max_halfwidth_um", "Max. shaft half-width (um)"),
    ("shaft_margin_um", "Shaft margin (um; signed)"),
    ("spine_reach_um", "Max. spine reach (um)"),
    ("min_area_um2", "Min. candidate area (um2)"),
    ("max_area_um2", "Max. candidate area (um2; 0=none)"),
]


def restore_session(path):
    saved = json.loads(Path(path).read_text(encoding="utf-8"))
    p = Settings(**saved["settings"])
    if saved["source"] == "SYNTHETIC_DEMO":
        data, _, _ = demo()
    else:
        sel = saved["source_selection"]
        data = load_czi(saved["source"], scene=sel["scene"], channel=sel["channel"],
                        time=sel["time"], z_start=sel["z_start"],
                        z_end=sel["z_end_inclusive"], projection=sel["projection"])
    if list(data.image.shape) != saved["source_shape_yx"] or not np.allclose(
            [data.sx, data.sy], [saved["sx_um_per_pixel"], saved["sy_um_per_pixel"]], rtol=1e-9):
        raise ValueError("Image shape or calibration differs from the saved session.")
    if saved.get("projection_sha256") != image_digest(data.image):
        raise ValueError("Projection pixels differ from the saved session; use the same CZI.")
    return data, p, saved


class App:
    def __init__(self, root):
        self.root = root
        root.title("CZI spine-positive area / dendrite length")
        root.geometry("1520x960")
        root.minsize(1120, 740)
        self.data = None
        self.result = None
        self.points = []
        self.exclusions = []
        self.brush_strokes = []
        self.polygon = []
        self.upper_points, self.lower_points = [], []
        self.cached_preprocess = None
        self.cached_geometry = None
        self.last_geometry_key = None
        self.mode = tk.StringVar(value="trace")
        self.vars = {}
        self.pending = None
        sidebar = ttk.Frame(root)
        sidebar.pack(side="left", fill="y")
        scroller = tk.Canvas(sidebar, width=365, highlightthickness=0)
        scroll = ttk.Scrollbar(sidebar, orient="vertical", command=scroller.yview)
        scroll.pack(side="right", fill="y")
        scroller.pack(side="left", fill="y", expand=True)
        scroller.configure(yscrollcommand=scroll.set)
        controls = ttk.Frame(scroller, padding=10)
        control_window = scroller.create_window((0, 0), window=controls, anchor="nw")
        controls.bind("<Configure>", lambda event: scroller.configure(scrollregion=scroller.bbox("all")))
        scroller.bind("<Configure>", lambda event: scroller.itemconfigure(control_window, width=event.width))
        main = ttk.Frame(root)
        main.pack(side="left", fill="both", expand=True)
        ttk.Button(controls, text="1  Open CZI", command=self.open_file).pack(fill="x", pady=3)
        tools = ttk.Frame(controls)
        tools.pack(fill="x")
        ttk.Button(tools, text="Synthetic demo", command=self.load_demo).pack(side="left", expand=True, fill="x")
        ttk.Button(tools, text="Restore session", command=self.load_session).pack(side="left", expand=True, fill="x")
        self.info = ttk.Label(controls, text="Open an original CZI.\nAll indices are zero-based.", wraplength=270)
        self.info.pack(fill="x", pady=7)
        ttk.Separator(controls).pack(fill="x", pady=4)
        defaults = Settings()
        for key, title in PARAMETERS:
            row = ttk.Frame(controls)
            row.pack(fill="x", pady=2)
            ttk.Label(row, text=title, width=29).pack(side="left")
            var = tk.StringVar(value=f"{getattr(defaults, key):g}")
            entry = ttk.Entry(row, textvariable=var, width=9)
            entry.pack(side="right")
            entry.bind("<Return>", lambda event: self.run_analysis())
            self.vars[key] = var
        self.threshold_value = tk.DoubleVar(value=defaults.threshold)
        self.threshold_scale = tk.Scale(controls, from_=1, to=65535, orient="horizontal",
                                        variable=self.threshold_value,
                                        showvalue=False, resolution=-1, digits=9,
                                        command=self.threshold_changed)
        self.threshold_scale.pack(fill="x")
        ttk.Button(controls, text="Suggest Otsu threshold (ROI)", command=self.otsu).pack(fill="x", pady=2)
        ttk.Button(controls, text="2  Update analysis", command=self.run_analysis).pack(fill="x", pady=4)
        ttk.Label(controls, text="Mode:", font=("Microsoft YaHei", 10, "bold")).pack(anchor="w")
        for title, value in (("Source: trace centerline (proximal to distal)", "trace"),
                             ("Straight: upper edge (negative side)", "upper"),
                             ("Straight: lower edge (positive side)", "lower"),
                             ("Source: exclusion polygon", "exclude")):
            ttk.Radiobutton(controls, text=title, variable=self.mode, value=value).pack(anchor="w")
        row = ttk.Frame(controls)
        row.pack(fill="x", pady=3)
        ttk.Button(row, text="Undo point", command=self.undo_point).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Finish exclusion", command=self.finish_polygon).pack(side="left", fill="x", expand=True)
        row = ttk.Frame(controls)
        row.pack(fill="x", pady=3)
        ttk.Button(row, text="Clear trace", command=self.clear_trace).pack(side="left", fill="x", expand=True)
        ttk.Button(row, text="Auto shaft edges", command=self.clear_edges).pack(side="left", fill="x", expand=True)
        ttk.Button(controls, text="Undo last exclusion", command=self.undo_exclusion).pack(fill="x", pady=2)
        ttk.Button(controls, text="3  Save results + source annotation", command=self.export).pack(fill="x", pady=5)
        ttk.Button(controls, text="Open guide", command=self.open_guide).pack(fill="x")
        ttk.Label(controls, text="Green: shaft; magenta: spine candidates\nAmber: other positive; blue: excluded\nDisable zoom/pan before clicking points.",
                  wraplength=280).pack(anchor="w", pady=7)
        self.figure = Figure(figsize=(12, 8), layout="constrained", facecolor="#f6f7fa")
        gs = self.figure.add_gridspec(2, 2, width_ratios=[1, 1.2], height_ratios=[1, .65])
        self.source_ax = self.figure.add_subplot(gs[:, 0])
        self.crop_ax = self.figure.add_subplot(gs[0, 1])
        self.straight_ax = self.figure.add_subplot(gs[1, 1])
        self.canvas = FigureCanvasTkAgg(self.figure, master=main)
        self.canvas.get_tk_widget().pack(fill="both", expand=True)
        self.toolbar = NavigationToolbar2Tk(self.canvas, main, pack_toolbar=False)
        self.toolbar.pack(fill="x")
        self.canvas.mpl_connect("button_press_event", self.clicked)
        self.status = tk.StringVar(value="Trace one shaft near a branch; cover at least offset + segment length.")
        ttk.Label(main, textvariable=self.status, wraplength=1080,
                  font=("Microsoft YaHei", 11)).pack(fill="x", padx=8, pady=6)
        self.warning = tk.StringVar(value="2D positive candidate area; review shaft edges and overlapping structures.")
        ttk.Label(main, textvariable=self.warning, wraplength=1080,
                  foreground="#9b5900").pack(fill="x", padx=8, pady=5)
        self.draw()

    def get_settings(self):
        try:
            p = Settings(**{key: float(v.get()) for key, v in self.vars.items()})
        except ValueError:
            raise ValueError("All parameters must be numeric.") from None
        p.validate()
        return p

    def show_error(self, exc):
        messagebox.showerror("Cannot complete", str(exc), parent=self.root)

    def choose_selection(self, scenes):
        dlg = tk.Toplevel(self.root)
        dlg.title("Select CZI scene / channel / time / Z")
        dlg.transient(self.root)
        dlg.grab_set()
        frame = ttk.Frame(dlg, padding=16)
        frame.pack(fill="both", expand=True)
        scenevar = tk.StringVar(value=str(scenes[0]["scene"]))
        ttk.Label(frame, text="Scene (original CZI index)").grid(row=0, column=0, sticky="w", pady=4)
        combo = ttk.Combobox(frame, textvariable=scenevar,
                             values=[str(s["scene"]) for s in scenes], state="readonly")
        combo.grid(row=0, column=1)
        values, spins = {}, {}
        for row, (key, title) in enumerate((("channel", "Channel index"), ("time", "Time index"),
                                           ("z_start", "Z start index"), ("z_end", "Z end index (inclusive)")), 1):
            ttk.Label(frame, text=title).grid(row=row, column=0, sticky="w", pady=4)
            values[key] = tk.StringVar(value="0")
            spins[key] = ttk.Spinbox(frame, from_=0, to=0, textvariable=values[key], width=20)
            spins[key].grid(row=row, column=1)
        proj = tk.StringVar(value="max")
        ttk.Label(frame, text="max: maximum projection; single: Z start only").grid(row=5, column=0, sticky="w", pady=4)
        ttk.Combobox(frame, textvariable=proj, values=["max", "single"], state="readonly").grid(row=5, column=1)
        explanation = ttk.Label(frame, wraplength=540)
        explanation.grid(row=6, column=0, columnspan=2, pady=10)
        def changed(event=None):
            s = next(s for s in scenes if s["scene"] == int(scenevar.get()))
            for key, axis in (("channel", "C"), ("time", "T"), ("z_start", "Z"), ("z_end", "Z")):
                last = s["sizes"].get(axis, 1) - 1
                spins[key].configure(to=last)
                values[key].set(str(last if key == "z_end" else 0))
            explanation.configure(text=f"Dimensions: {s['sizes']}\nChannels: {s['channels']}\n"
                                  f"X/Y calibration (um/pixel): {s['mpp']}\nZero-based indices; keep channel and projection consistent.")
        combo.bind("<<ComboboxSelected>>", changed)
        changed()
        selected = []
        def accept():
            try:
                selected.append({"scene": int(scenevar.get()), **{k: int(v.get()) for k, v in values.items()},
                                 "projection": proj.get()})
                dlg.destroy()
            except ValueError:
                messagebox.showerror("Invalid index", "Indices must be integers.", parent=dlg)
        ttk.Button(frame, text="Load selection", command=accept).grid(row=7, column=0, columnspan=2, sticky="ew")
        self.root.wait_window(dlg)
        return selected[0] if selected else None

    def set_data(self, data, p=None):
        self.data = data
        self.selection_box = None
        self.points, self.exclusions, self.polygon = [], [], []
        self.brush_strokes = []
        self.upper_points, self.lower_points = [], []
        self.result = None
        self.cached_preprocess = self.cached_geometry = self.last_geometry_key = None
        if p is not None:
            for key, value in vars(p).items():
                self.vars[key].set(f"{value:g}")
        self.threshold_scale.configure(to=max(1, float(np.max(data.image))))
        self.threshold_value.set(float(self.vars["threshold"].get()))
        label = Path(data.source).name
        sel = data.selection
        self.info.configure(text=f"{label}\n{data.image.shape[1]} × {data.image.shape[0]} px\n"
                            f"X={data.sx:.7g}，Y={data.sy:.7g} µm/px\n"
                            f"C={sel.get('channel', '-')}，T={sel.get('time', '-')}，"
                            f"Z={sel.get('z_start', '-')}..{sel.get('z_end_inclusive', '-')} {sel.get('projection', '')}")
        self.mode.set("trace")
        self.status.set("Trace near a branch, avoiding the enlarged junction; right-click to undo.")
        self.draw()

    def open_file(self):
        path = filedialog.askopenfilename(title="Open original CZI", filetypes=[("Zeiss CZI", "*.czi")])
        if not path:
            return
        try:
            self.status.set("Reading CZI metadata...")
            self.root.update_idletasks()
            selection = self.choose_selection(inspect_czi(path))
            if selection is None:
                return
            self.status.set("Loading selected scene/channel/Z...")
            self.root.update_idletasks()
            data = load_czi(path, **selection)
            p = self.get_settings()
            p.threshold = suggest_threshold(preprocess(data, p))
            self.set_data(data, p)
        except Exception as exc:
            self.show_error(exc)

    def load_demo(self):
        data, points, p = demo()
        self.set_data(data, p)
        self.points = points
        self.run_analysis()
        self.warning.set("Synthetic demonstration; biological segmentation accuracy is unvalidated.")

    def load_session(self):
        path = filedialog.askopenfilename(title="Restore analysis session", filetypes=[("Session JSON", "session.json *.json")])
        if not path:
            return
        try:
            data, p, saved = restore_session(path)
            self.set_data(data, p)
            self.points = saved["trace_xy_source_px"]
            self.selection_box = saved.get("selection_box_source_px_xyxy")
            self.exclusions = saved["exclusion_polygons_xy_source_px"]
            self.upper_points = saved["upper_boundary_points_s_t_um"]
            self.lower_points = saved["lower_boundary_points_s_t_um"]
            self.brush_strokes = saved.get("brush_strokes", [])
            self.run_analysis()
        except Exception as exc:
            self.show_error(exc)

    def corrected(self, p):
        key = (id(self.data), p.gaussian_um, p.background_radius_um)
        if self.cached_preprocess is None or self.cached_preprocess[0] != key:
            self.cached_preprocess = (key, preprocess(self.data, p))
        return self.cached_preprocess[1]

    def selected_geometry(self, p):
        box = getattr(self, "selection_box", None)
        key = (tuple(map(tuple, self.points)), p.length_um, p.start_offset_um, p.half_width_um,
               tuple(box) if box else None)
        if self.cached_geometry is None or self.cached_geometry[0] != key:
            if self.last_geometry_key is not None and self.last_geometry_key[1:] != key[1:]:
                self.upper_points, self.lower_points = [], []
            self.cached_geometry = (key, geometry(self.data, self.points, p, box))
            self.last_geometry_key = key
        return self.cached_geometry[1]

    def run_analysis(self, quiet=False):
        if self.data is None:
            if not quiet:
                self.show_error("Open a CZI or synthetic demo first.")
            return False
        try:
            p = self.get_settings()
            g = self.selected_geometry(p)
            self.result = analyze(self.data, self.corrected(p), g, p, self.exclusions,
                                  self.upper_points, self.lower_points, self.brush_strokes)
            self.status.set(f"Spine candidate area {self.result.summary['spine_area_um2']:.4f} µm² ÷ "
                            f"Centerline length {g.length:.4f} µm = "
                            f"{self.result.summary['area_per_length_um2_per_um']:.5f} µm²/µm")
            self.warning.set("; ".join(self.result.warnings) or "Review magenta candidates and green shaft edges before export.")
            self.draw()
            return True
        except Exception as exc:
            self.result = None
            self.status.set(str(exc))
            if not quiet:
                self.show_error(exc)
            self.draw()
            return False

    def threshold_changed(self, value):
        self.vars["threshold"].set(f"{float(value):g}")
        if self.pending:
            self.root.after_cancel(self.pending)
        self.pending = self.root.after(250, self.threshold_apply)

    def threshold_apply(self):
        self.pending = None
        if len(self.points) >= 2:
            self.run_analysis(quiet=True)

    def otsu(self):
        if self.data is None:
            return
        try:
            p = self.get_settings()
            g = self.selected_geometry(p) if len(self.points) >= 2 else None
            value = suggest_threshold(self.corrected(p), g)
            self.threshold_value.set(value)
            self.vars["threshold"].set(f"{value:.7g}")
            if g is not None:
                self.run_analysis()
        except Exception as exc:
            self.show_error(exc)

    def clicked(self, event):
        if self.data is None or event.inaxes is None or self.toolbar.mode:
            return
        if event.button == 3:
            self.undo_point()
            return
        if event.button != 1:
            return
        point = [float(event.xdata), float(event.ydata)]
        mode = self.mode.get()
        if mode in ("trace", "exclude") and event.inaxes in (self.source_ax, self.crop_ax):
            if event.inaxes is self.crop_ax:
                if self.result is None:
                    return
                point[0] += self.result.geometry.bbox[0]
                point[1] += self.result.geometry.bbox[1]
            h, w = self.data.image.shape
            if not (0 <= point[0] < w and 0 <= point[1] < h):
                return
            if mode == "trace":
                self.points.append(point)
                self.upper_points, self.lower_points = [], []
                self.result = None
                self.status.set(f"Selected {len(self.points)} trace points; click Update analysis when the trace is long enough.")
            else:
                self.polygon.append(point)
                self.status.set(f"Exclusion polygon: {len(self.polygon)} points; click Finish exclusion after at least three.")
            self.draw(preserve=True)
        elif mode in ("upper", "lower") and event.inaxes is self.straight_ax and self.result is not None:
            if (mode == "upper" and point[1] >= 0) or (mode == "lower" and point[1] <= 0):
                self.status.set("Upper edge: negative side; lower edge: positive side.")
                return
            target = self.upper_points if mode == "upper" else self.lower_points
            target.append(point)
            self.run_analysis(quiet=True)

    def undo_point(self):
        mode = self.mode.get()
        target = {"trace": self.points, "exclude": self.polygon,
                  "upper": self.upper_points, "lower": self.lower_points}[mode]
        if target:
            target.pop()
        if mode == "trace":
            self.result = None
            self.upper_points, self.lower_points = [], []
            self.draw(preserve=True)
        elif mode == "exclude":
            self.draw(preserve=True)
        else:
            self.run_analysis(quiet=True)

    def finish_polygon(self):
        if len(self.polygon) < 3:
            self.show_error("Select at least three polygon vertices in the source or crop.")
            return
        self.exclusions.append(self.polygon.copy())
        self.polygon = []
        self.run_analysis(quiet=True)

    def undo_exclusion(self):
        if self.exclusions:
            self.exclusions.pop()
        self.run_analysis(quiet=True)

    def clear_trace(self):
        self.points, self.polygon, self.upper_points, self.lower_points = [], [], [], []
        self.brush_strokes = []
        self.result = None
        self.cached_geometry = self.last_geometry_key = None
        self.status.set("Trace cleared; select points from proximal to distal.")
        self.draw()

    def clear_edges(self):
        self.upper_points, self.lower_points = [], []
        self.run_analysis(quiet=True)

    def draw(self, preserve=False):
        limits = (self.source_ax.get_xlim(), self.source_ax.get_ylim()) if preserve else None
        for ax in (self.source_ax, self.crop_ax, self.straight_ax):
            ax.clear()
        self.source_ax.set_title("Source: shaft green / spines magenta")
        self.crop_ax.set_title("Source-pixel crop / segmentation")
        self.straight_ax.set_title("Straight view: green shaft edges, magenta candidates")
        if self.data is None:
            self.source_ax.text(.5, .5, "Open CZI or select Synthetic demo", transform=self.source_ax.transAxes,
                                ha="center", va="center")
            self.canvas.draw_idle()
            return
        a = self.data.image
        high = max(1., float(np.percentile(a, 99.7)))
        if self.result is None:
            self.source_ax.imshow(a, cmap="gray", vmin=0, vmax=high)
        else:
            self.source_ax.imshow(source_annotation(self.data, self.result))
        if self.points:
            xy = np.asarray(self.points)
            self.source_ax.plot(xy[:, 0], xy[:, 1], "o-", color="cyan", ms=3, lw=1)
        for poly in self.exclusions + ([self.polygon] if self.polygon else []):
            xy = np.asarray(poly + [poly[0]])
            self.source_ax.plot(xy[:, 0], xy[:, 1], color="#48a6ff", lw=1)
        if self.result is not None:
            r = self.result
            g = r.geometry
            self.source_ax.plot(g.path[:, 0], g.path[:, 1], color="lime", lw=2)
            self.source_ax.scatter(*g.path[0], color="lime", s=45, marker="s")
            self.source_ax.scatter(*g.path[-1], color="magenta", s=45, marker="s")
            self.crop_ax.imshow(r.raw, cmap="gray", vmin=0, vmax=high)
            self.crop_ax.imshow(overlay(r))
            self.crop_ax.plot(g.path[:, 0] - g.bbox[0], g.path[:, 1] - g.bbox[1], color="cyan", lw=.7)
            ext = (0, g.length, g.t_axis[-1], g.t_axis[0])
            self.straight_ax.imshow(g.straighten(a.astype(np.float32)), cmap="gray", vmin=0,
                                    vmax=high, extent=ext, aspect="equal")
            rgba = overlay(r)
            sr = np.stack([g.straighten(rgba[:, :, i], order=0, crop=True) for i in range(4)], axis=-1)
            self.straight_ax.imshow(sr, extent=ext, aspect="equal")
            self.straight_ax.plot(g.s_axis, -r.upper, color="lime", lw=1)
            self.straight_ax.plot(g.s_axis, r.lower, color="lime", lw=1)
            self.straight_ax.axhline(0, color="cyan", ls="--", lw=.6)
            for pts in (self.upper_points, self.lower_points):
                if pts:
                    pts = np.asarray(pts)
                    self.straight_ax.scatter(pts[:, 0], pts[:, 1], s=25, c="cyan")
            self.straight_ax.set(xlabel="Centerline distance (um)", ylabel="Normal offset (um)")
        if limits:
            self.source_ax.set_xlim(limits[0])
            self.source_ax.set_ylim(limits[1])
        self.canvas.draw_idle()

    def export(self):
        if self.polygon:
            self.show_error("Finish the pending exclusion polygon or undo its vertices before export.")
            return
        if not self.run_analysis():
            return
        parent = filedialog.askdirectory(title="Choose parent folder for a new result directory")
        if not parent:
            return
        from datetime import datetime
        folder = Path(parent) / (Path(self.data.source).stem + "_" + datetime.now().strftime("%Y%m%d_%H%M%S_%f"))
        try:
            save_result(folder, self.data, self.points, self.get_settings(), self.result,
                        self.exclusions, self.upper_points, self.lower_points)
            self.status.set(f"Exported: {folder}")
            messagebox.showinfo("Exported", f"{folder}\n\nresults.csv: area / length\nsession.json: replayable settings\n"
                                "source_annotated.png/.tif: full-size annotation\nsource_spine_mask.tif: calibrated spine mask", parent=self.root)
        except Exception as exc:
            self.show_error(exc)

    def open_guide(self):
        os.startfile(Path(__file__).parent / "IMAGEJ_GUIDE.md")


def main():
    parser = argparse.ArgumentParser(description="Calibrated CZI spine-positive area / dendrite length")
    parser.add_argument("--demo", action="store_true")
    parser.add_argument("--session", type=Path, help="Recompute a saved session without the GUI")
    parser.add_argument("--out", type=Path, help="New output directory (must not exist)")
    parser.add_argument("--gui-smoke", action="store_true", help=argparse.SUPPRESS)
    args = parser.parse_args()
    if args.session:
        if args.out is None:
            parser.error("--session requires --out.")
        data, p, saved = restore_session(args.session)
        points = saved["trace_xy_source_px"]
        excluded = saved["exclusion_polygons_xy_source_px"]
        up, down = saved["upper_boundary_points_s_t_um"], saved["lower_boundary_points_s_t_um"]
        result = analyze(data, preprocess(data, p), geometry(data, points, p, saved.get("selection_box_source_px_xyxy")), p,
                         excluded, up, down, saved.get("brush_strokes", []))
        save_result(args.out, data, points, p, result, excluded, up, down)
        print(json.dumps(result.summary, ensure_ascii=False, indent=2))
        return
    root = tk.Tk()
    app = App(root)
    if args.demo or args.gui_smoke:
        app.load_demo()
    if args.gui_smoke:
        root.update()
        assert app.result is not None and app.result.summary["spine_area_um2"] > 0
        app.figure.savefig(Path(__file__).parent / "gui_preview.png", dpi=120)
        root.destroy()
        print("GUI smoke check passed")
    else:
        root.mainloop()


if __name__ == "__main__":
    main()
