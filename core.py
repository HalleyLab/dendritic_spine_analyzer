"""Calibrated 2-D dendritic protrusion area measurement on source pixels."""
from __future__ import annotations

import csv
import hashlib
import json
import math
import platform
import sys
from dataclasses import asdict, dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from xml.etree import ElementTree as ET

import numpy as np
from scipy import ndimage as ndi
from skimage.filters import threshold_otsu
import tifffile

VENDOR = Path(__file__).parent / "vendor"
if VENDOR.exists():
    sys.path.insert(0, str(VENDOR))
VERSION = "1.7.0"
STRUCTURE = np.ones((3, 3), bool)


@dataclass
class ImageData:
    image: np.ndarray
    sx: float
    sy: float
    source: str
    selection: dict = field(default_factory=dict)
    metadata_xml: str = ""

    def __post_init__(self):
        if self.image.ndim != 2 or min(self.image.shape) < 2:
            raise ValueError("A single-channel 2D image is required.")
        if not np.isfinite([self.sx, self.sy]).all() or min(self.sx, self.sy) <= 0:
            raise ValueError("X/Y pixel sizes must be positive (um/pixel).")
        if not np.isfinite(self.image).all():
            raise ValueError("The image contains NaN or Inf.")


@dataclass
class Settings:
    length_um: float = 20.0
    start_offset_um: float = 2.0
    half_width_um: float = 6.0
    threshold: float = 1000.0  # Fluorescence-positive pixels, including dim spine necks.
    shaft_threshold: float = 0.0  # 0 follows threshold; a positive value fixes the shaft reference.
    gaussian_um: float = 0.0
    background_radius_um: float = 0.0
    shaft_halfwidth_um: float = 0.0
    shaft_window_um: float = 2.0
    shaft_min_halfwidth_um: float = 0.15
    shaft_max_halfwidth_um: float = 2.5
    shaft_margin_um: float = 0.0
    spine_reach_um: float = 3.0
    min_area_um2: float = 0.03
    max_area_um2: float = 4.0

    def validate(self):
        vals = asdict(self)
        if not np.isfinite(list(vals.values())).all():
            raise ValueError("Parameters must be finite numbers.")
        if min(self.length_um, self.half_width_um, self.threshold,
               self.shaft_window_um, self.shaft_min_halfwidth_um,
               self.shaft_max_halfwidth_um, self.spine_reach_um) <= 0:
            raise ValueError("Length, ROI half-width, threshold, shaft width/window and reach must be positive.")
        nonnegative = (self.start_offset_um, self.gaussian_um,
                       self.background_radius_um, self.shaft_halfwidth_um,
                       self.min_area_um2, self.max_area_um2, self.shaft_threshold)
        if min(nonnegative) < 0:
            raise ValueError("Offset, smoothing, background radius, fixed half-width and area limits cannot be negative.")
        if self.shaft_min_halfwidth_um > self.shaft_max_halfwidth_um:
            raise ValueError("Minimum shaft half-width cannot exceed the maximum.")
        if self.shaft_max_halfwidth_um >= self.half_width_um:
            raise ValueError("ROI half-width must exceed maximum shaft half-width.")
        if self.shaft_halfwidth_um >= self.half_width_um:
            raise ValueError("Fixed shaft half-width must be smaller than ROI half-width.")
        if self.max_area_um2 and self.max_area_um2 < self.min_area_um2:
            raise ValueError("Maximum candidate area cannot be smaller than the minimum.")


def scaling_from_xml(xml: str) -> tuple[float, float]:
    """CZI Scaling/Distance/Value stores metres; display-unit fields are not used."""
    root = ET.fromstring(xml)
    for node in root.iter():
        node.tag = node.tag.rsplit("}", 1)[-1]
    values = {}
    for item in root.findall(".//Scaling/Items/Distance"):
        if item.get("Id") in ("X", "Y"):
            values[item.get("Id")] = float(item.findtext("Value", "nan")) * 1e6
    xy = (values.get("X", float("nan")), values.get("Y", float("nan")))
    if not np.isfinite(xy).all() or min(xy) <= 0:
        raise ValueError("CZI has no valid X/Y calibration; verify it in ZEN or Bio-Formats.")
    return xy


def inspect_czi(path) -> list[dict]:
    from czifile import CziFile
    with CziFile(path) as czi:
        result = []
        for key, view in czi.scenes.items():
            result.append({"scene": int(key), "sizes": dict(view.sizes),
                           "starts": dict(zip(view.dims, view.start)),
                           "mpp": view.mpp, "channels": list(view.channels)})
        return result


def load_czi(path, scene=0, channel=0, time=0, z_start=0, z_end=None,
             projection="max") -> ImageData:
    """Channel/time/Z arguments are zero-based positions within the chosen scene."""
    from czifile import CziFile
    if projection not in ("max", "single"):
        raise ValueError("Projection must be max or single.")
    with CziFile(path) as czi:
        view = czi.scenes[int(scene)]
        original_sizes = dict(view.sizes)
        starts = dict(zip(view.dims, view.start))
        choice = {}
        for axis, index in (("C", channel), ("T", time)):
            size = original_sizes.get(axis, 1)
            if not isinstance(index, int) or not 0 <= index < size:
                raise ValueError(f"{axis} index is outside 0..{size - 1}.")
            if axis in original_sizes:
                choice[axis] = starts[axis] + index
        nz = original_sizes.get("Z", 1)
        z_end = nz - 1 if z_end is None else z_end
        if not 0 <= z_start <= z_end < nz:
            raise ValueError(f"Z range is outside 0..{nz - 1}; the end index is inclusive.")
        if "Z" in original_sizes:
            a = starts["Z"] + z_start
            choice["Z"] = a if projection == "single" else slice(a, starts["Z"] + z_end + 1)
        view = view(**choice)
        for axis, size in view.sizes.items():
            if axis not in ("Z", "Y", "X") and size > 1:
                raise ValueError(f"Unselected {axis} dimension ({size}); RGB, illumination and phase are not merged automatically.")
        xml = czi.metadata()
        sx, sy = view.mpp or scaling_from_xml(xml)
        arr = view.asarray()
        dims = list(view.dims)
        for i in reversed(range(len(dims))):
            if dims[i] not in ("Z", "Y", "X"):
                arr = np.take(arr, 0, axis=i)
                dims.pop(i)
        if "Z" in dims:
            i = dims.index("Z")
            arr = np.max(arr, axis=i) if projection == "max" else np.take(arr, 0, axis=i)
            dims.pop(i)
        if set(dims) != {"Y", "X"}:
            raise ValueError(f"Unsupported spatial dimensions: {dims}")
        arr = np.transpose(arr, (dims.index("Y"), dims.index("X")))
        selection = {"scene": scene, "channel": channel, "time": time,
                     "z_start": z_start, "z_end_inclusive": z_end,
                     "projection": projection, "source_sizes": original_sizes,
                     "source_bbox_global_px": list(view.bbox),
                     "calibration_source": "CZI full-resolution view.mpp (µm/pixel)"}
        return ImageData(arr, sx, sy, str(Path(path).resolve()), selection, xml)


def preprocess(data: ImageData, p: Settings) -> np.ndarray:
    p.validate()
    a = data.image.astype(np.float32)
    if p.gaussian_um:
        a = ndi.gaussian_filter(a, (p.gaussian_um / data.sy, p.gaussian_um / data.sx),
                                mode="reflect")
    if p.background_radius_um:
        ry = math.ceil(p.background_radius_um / data.sy)
        rx = math.ceil(p.background_radius_um / data.sx)

        background = ndi.grey_opening(a, size=(2 * ry + 1, 2 * rx + 1), mode="reflect")
        a = np.maximum(a - background, 0)
    return a


def trim_path(points, sx, sy, length, offset=0.0):
    a = np.asarray(points, float)
    if a.ndim != 2 or a.shape[1] != 2 or len(a) < 2 or not np.isfinite(a).all():
        raise ValueError("Select at least two valid centerline points along one shaft.")
    if min(sx, sy, length) <= 0 or offset < 0:
        raise ValueError("Pixel sizes and length must be positive; offset must be nonnegative.")
    physical = a * [sx, sy]
    keep = np.r_[True, np.linalg.norm(np.diff(physical, axis=0), axis=1) > 1e-9]
    a, physical = a[keep], physical[keep]
    cum = np.r_[0, np.cumsum(np.linalg.norm(np.diff(physical, axis=0), axis=1))]
    end = offset + length
    if len(a) < 2 or cum[-1] < end - 1e-8:
        raise ValueError(f"Centerline length is {cum[-1]:.3f} um; at least {end:.3f} um is required. Add distal points.")
    s = np.r_[offset, cum[(cum > offset + 1e-9) & (cum < end - 1e-9)], end]
    path = np.column_stack([np.interp(s, cum, physical[:, i]) for i in (0, 1)])

    for i in range(len(path) - 1):
        for j in range(i + 2, len(path) - 1):
            u, v = path[i + 1] - path[i], path[j + 1] - path[j]
            cross = u[0] * v[1] - u[1] * v[0]
            if abs(cross) < 1e-10:
                continue
            d = path[j] - path[i]
            t = (d[0] * v[1] - d[1] * v[0]) / cross
            q = (d[0] * u[1] - d[1] * u[0]) / cross
            if 0 <= t <= 1 and 0 <= q <= 1:
                raise ValueError("The centerline crosses itself; trace one shaft in order.")
    return path / [sx, sy]


@dataclass
class Geometry:
    path: np.ndarray
    bbox: tuple[int, int, int, int]
    s: np.ndarray
    t: np.ndarray
    corridor: np.ndarray
    grid_x: np.ndarray
    grid_y: np.ndarray
    s_axis: np.ndarray
    t_axis: np.ndarray
    length: float
    clipped: bool
    selection_box: tuple[int, int, int, int] | None = None

    def straighten(self, image, order=1, crop=False):
        x0, y0, _, _ = self.bbox
        return ndi.map_coordinates(image, [self.grid_y - (y0 if crop else 0),
                                           self.grid_x - (x0 if crop else 0)],
                                   order=order, mode="constant", cval=0,
                                   prefilter=False)


def geometry(data: ImageData, points, p: Settings, selection_box=None) -> Geometry:
    p.validate()
    path = trim_path(points, data.sx, data.sy, p.length_um, p.start_offset_um)
    h, w = data.image.shape
    if (path < 0).any() or (path[:, 0] > w - 1).any() or (path[:, 1] > h - 1).any():
        raise ValueError("The centerline must stay inside the image.")
    phys = path * [data.sx, data.sy]
    v = np.diff(phys, axis=0)
    lengths = np.linalg.norm(v, axis=1)
    cum = np.r_[0, np.cumsum(lengths)]
    x0, y0 = np.floor(path.min(axis=0) - p.half_width_um / np.array([data.sx, data.sy]) - 1).astype(int)
    x1, y1 = np.ceil(path.max(axis=0) + p.half_width_um / np.array([data.sx, data.sy]) + 2).astype(int)
    clipped = x0 < 0 or y0 < 0 or x1 > w or y1 > h
    x0, y0, x1, y1 = max(0, x0), max(0, y0), min(w, x1), min(h, y1)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    xx, yy = xx * data.sx, yy * data.sy
    best = np.full(xx.shape, np.inf)
    longitudinal = np.zeros(xx.shape)
    transverse = np.zeros(xx.shape)
    endpoint_ok = np.zeros(xx.shape, bool)
    for i, (start, delta, seglen) in enumerate(zip(phys[:-1], v, lengths)):
        dx, dy = xx - start[0], yy - start[1]
        u = (dx * delta[0] + dy * delta[1]) / seglen ** 2
        clipped_u = np.clip(u, 0, 1)
        ex, ey = dx - clipped_u * delta[0], dy - clipped_u * delta[1]
        d2 = ex * ex + ey * ey
        use = d2 < best
        best[use] = d2[use]
        longitudinal[use] = cum[i] + clipped_u[use] * seglen
        sign = np.where(delta[0] * dy - delta[1] * dx < 0, -1, 1)
        transverse[use] = (np.sqrt(d2) * sign)[use]
        ok = ((u >= 0) if i == 0 else np.ones_like(u, bool))
        if i == len(lengths) - 1:
            ok &= u < 1
        endpoint_ok[use] = ok[use]
    corridor = endpoint_ok & (best <= p.half_width_um ** 2)
    if selection_box is not None:
        box = np.asarray(selection_box, float)
        if box.shape != (4,) or not np.isfinite(box).all():
            raise ValueError("Select a valid rectangle on the source image.")
        bx0, by0, bx1, by1 = box
        if not (0 <= bx0 < bx1 <= w and 0 <= by0 < by1 <= h):
            raise ValueError("The selected rectangle must stay inside the source image.")
        if np.any((path[:, 0] < bx0) | (path[:, 0] >= bx1) | (path[:, 1] < by0) | (path[:, 1] >= by1)):
            raise ValueError("The measured centerline must stay inside the selected rectangle.")
        corridor &= (xx / data.sx >= bx0) & (xx / data.sx < bx1) & (yy / data.sy >= by0) & (yy / data.sy < by1)
        selection_box = tuple(map(float, box))
    step = min(data.sx, data.sy)
    n = max(2, math.ceil(p.length_um / step))
    sa = (np.arange(n) + 0.5) * p.length_um / n
    nr = max(1, math.ceil(p.half_width_um / step))
    ta = np.linspace(-p.half_width_um, p.half_width_um, 2 * nr + 1)
    seg = np.minimum(np.searchsorted(cum, sa, side="right") - 1, len(v) - 1)
    frac = (sa - cum[seg]) / lengths[seg]
    centers = phys[seg] + frac[:, None] * v[seg]
    normals = np.column_stack([-v[seg, 1], v[seg, 0]]) / lengths[seg, None]
    coords = centers[None, :, :] + ta[:, None, None] * normals[None, :, :]
    return Geometry(path, tuple(map(int, (x0, y0, x1, y1))), longitudinal, transverse,
                    corridor, coords[:, :, 0] / data.sx, coords[:, :, 1] / data.sy,
                    sa, ta, float(cum[-1]), clipped, selection_box)


def _manual_edge(points, sa, default, half_width):
    if not points:
        return default
    a = np.asarray(points, float)
    if not np.isfinite(a).all() or (a[:, 0] < 0).any() or (a[:, 0] > sa[-1] + sa[0]).any():
        raise ValueError("Manual edge points exceed the centerline length.")

    unique = dict((float(s), abs(float(t))) for s, t in a)
    x = np.array(sorted(unique))
    y = np.array([unique[s] for s in x])
    if (y <= 0).any() or (y >= half_width).any():
        raise ValueError("Manual edges must be on the correct side and inside the ROI.")
    return np.interp(sa, x, y)


@dataclass
class Result:
    geometry: Geometry
    raw: np.ndarray
    corrected: np.ndarray
    positive: np.ndarray
    attached: np.ndarray
    shaft_domain: np.ndarray
    shaft: np.ndarray
    candidates: np.ndarray
    spines: np.ndarray
    excluded: np.ndarray
    upper: np.ndarray
    lower: np.ndarray
    components: list[dict]
    warnings: list[str]
    summary: dict
    automatic_spines: np.ndarray | None = None
    manual_added: np.ndarray | None = None
    manual_removed: np.ndarray | None = None
    brush_strokes: list = field(default_factory=list)


def brush_overrides(data, g, strokes):
    """Paint source-pixel centers using a calibrated round brush; later strokes win."""
    out = np.zeros(g.corridor.shape, np.int8)
    gx0, gy0, gx1, gy1 = g.bbox
    h, w = data.image.shape
    scale = np.array([data.sx, data.sy])
    for stroke in strokes:
        mode = stroke["mode"]
        radius = float(stroke["radius_um"])
        points = np.asarray(stroke["points_xy_source_px"], float)
        if mode not in ("add", "remove") or not np.isfinite(radius) or radius <= 0:
            raise ValueError("Brush mode must be add/remove and radius must be positive.")
        if points.ndim != 2 or points.shape[1] != 2 or len(points) == 0 or not np.isfinite(points).all():
            raise ValueError("Brush strokes require finite source X/Y points.")
        if (points < 0).any() or (points[:, 0] > w-1).any() or (points[:, 1] > h-1).any():
            raise ValueError("Brush points must stay inside the source image.")
        ends = points[1:] if len(points) > 1 else points
        for start, end in zip(points, ends):
            low = np.maximum([gx0, gy0], np.floor(np.minimum(start, end) - radius/scale)).astype(int)
            high = np.minimum([gx1, gy1], np.ceil(np.maximum(start, end) + radius/scale) + 1).astype(int)
            x0, y0 = low; x1, y1 = high
            if x0 >= x1 or y0 >= y1:
                continue
            yy, xx = np.ogrid[y0:y1, x0:x1]
            dx, dy = (xx-start[0])*data.sx, (yy-start[1])*data.sy
            vx, vy = (end-start)*scale
            length2 = vx*vx + vy*vy
            fraction = np.clip((dx*vx + dy*vy) / (length2 if length2 else 1), 0, 1)
            hit = (dx-fraction*vx)**2 + (dy-fraction*vy)**2 <= radius**2 + 1e-12
            patch = out[y0-gy0:y1-gy0, x0-gx0:x1-gx0]
            patch[hit] = 1 if mode == "add" else -1
    out[~g.corridor] = 0
    return out


def analyze(data: ImageData, corrected, g: Geometry, p: Settings,
            exclusions=None, upper_points=None, lower_points=None, brush_strokes=None) -> Result:
    p.validate()
    from matplotlib.path import Path as Polygon
    x0, y0, x1, y1 = g.bbox
    raw = data.image[y0:y1, x0:x1]
    a = corrected[y0:y1, x0:x1]
    exclusion = np.zeros(raw.shape, bool)
    yy, xx = np.mgrid[y0:y1, x0:x1]
    for vertices in exclusions or []:
        v = np.asarray(vertices, float)
        if v.ndim != 2 or v.shape[1] != 2 or len(v) < 3 or not np.isfinite(v).all():
            raise ValueError("An exclusion polygon requires at least three valid points.")
        exclusion |= Polygon(v).contains_points(np.column_stack([xx.ravel(), yy.ravel()]),
                                                 radius=1e-7).reshape(raw.shape)
    # Threshold on source pixels; no hole filling or intensity rescaling for measurement.
    positive = (a >= p.threshold) & g.corridor & ~exclusion
    display = g.straighten(corrected)
    # Estimate the shaft independently so lowering the spine threshold does not widen it.
    above = display >= (p.shaft_threshold or p.threshold)
    mid = len(g.t_axis) // 2
    dt = g.t_axis[1] - g.t_axis[0]
    warnings = []
    if p.shaft_halfwidth_um > 0:
        up = np.full(len(g.s_axis), p.shaft_halfwidth_um)
        down = up.copy()
    else:

        up = np.full(len(g.s_axis), np.nan)
        down = up.copy()
        for j in range(len(g.s_axis)):
            if not above[mid, j]:
                continue
            u, d = mid, mid
            while u > 0 and above[u - 1, j]:
                u -= 1
            while d + 1 < len(g.t_axis) and above[d + 1, j]:
                d += 1
            up[j] = (mid - u + 0.5) * dt
            down[j] = (d - mid + 0.5) * dt
        valid = np.isfinite(up)
        if valid.sum() < 2:
            if not (upper_points and lower_points):
                raise ValueError("Insufficient positive centerline signal; lower the threshold, adjust the trace or set a fixed shaft half-width.")
            up[:] = down[:] = p.shaft_min_halfwidth_um
        else:
            for widths in (up, down):
                widths[:] = np.interp(g.s_axis, g.s_axis[valid], widths[valid])
                ds = g.length / len(g.s_axis)
                win = max(1, round(p.shaft_window_um / ds)) | 1
                widths[:] = ndi.median_filter(widths, size=win, mode="nearest")
        if valid.mean() < 0.9:
            warnings.append(f"Only {valid.mean():.0%} of centerline samples are positive; check threshold and trace.")
        if np.any(up > p.shaft_max_halfwidth_um) or np.any(down > p.shaft_max_halfwidth_um):
            warnings.append("Automatic shaft half-width reached its limit; inspect branches, thick shaft or merged spines.")
        up = np.clip(up, p.shaft_min_halfwidth_um, p.shaft_max_halfwidth_um)
        down = np.clip(down, p.shaft_min_halfwidth_um, p.shaft_max_halfwidth_um)
    up = np.clip(up + p.shaft_margin_um, p.shaft_min_halfwidth_um, p.half_width_um - dt / 2)
    down = np.clip(down + p.shaft_margin_um, p.shaft_min_halfwidth_um, p.half_width_um - dt / 2)
    up = _manual_edge(upper_points, g.s_axis, up, p.half_width_um)
    down = _manual_edge(lower_points, g.s_axis, down, p.half_width_um)
    mapped_up = np.interp(g.s, g.s_axis, up)
    mapped_down = np.interp(g.s, g.s_axis, down)
    band = (g.t >= -mapped_up) & (g.t <= mapped_down) & g.corridor & ~exclusion
    labels, _ = ndi.label(positive, STRUCTURE)
    # Only positive components connected to the shaft can become spine candidates.
    seeds = np.unique(labels[band & positive])
    seeds = seeds[seeds > 0]
    attached = np.isin(labels, seeds) & positive
    shaft = attached & band
    reach = (g.t >= -mapped_up - p.spine_reach_um) & (g.t <= mapped_down + p.spine_reach_um)

    candidate = attached & ~band
    lab, n = ndi.label(candidate, STRUCTURE)
    sizes = np.bincount(lab.ravel(), minlength=n + 1)
    keep = np.zeros(n + 1, bool)
    components = []
    for i, location in enumerate(ndi.find_objects(lab), 1):
        if location is None:
            continue
        loc = lab[location] == i
        area = float(sizes[i] * data.sx * data.sy)
        reason = "accepted"
        # Reject the entire component, rather than clipping a long branch into a spine.
        if np.any(loc & ~reach[location]):
            reason = "exceeds_reach"
        elif area < p.min_area_um2:
            reason = "below_min_area"
        elif p.max_area_um2 and area > p.max_area_um2:
            reason = "above_max_area"
        keep[i] = reason == "accepted"
        by, bx = np.nonzero(loc)
        components.append({"component_id": i, "area_um2": area,
                           "centroid_x_source_px": float(bx.mean() + location[1].start + x0),
                           "centroid_y_source_px": float(by.mean() + location[0].start + y0),
                           "status": reason})
    spine = keep[lab]
    automatic_spines = spine.copy() if brush_strokes else None
    manual_added = manual_removed = None
    if brush_strokes:
        edits = brush_overrides(data, g, brush_strokes)
        # Adding overrides automatic filters, but still requires positive signal outside the shaft.
        spine = (spine | ((edits == 1) & positive & ~band)) & (edits != -1)
        manual_added = spine & ~automatic_spines
        manual_removed = automatic_spines & ~spine
    rejected = sum(row["status"] != "accepted" for row in components)
    if rejected:
        warnings.append(f"Automatic filters rejected {rejected} protrusion components; see components.csv.")
    if g.clipped:
        warnings.append("The ROI is near an image edge; inspect truncated spines.")
    edge = g.corridor & ~ndi.binary_erosion(g.corridor, STRUCTURE)
    if np.any(spine & edge):
        warnings.append("Candidates touch the ROI boundary; outside area is excluded. Consider a wider ROI or another segment.")
    if np.any(positive & ~attached & ~spine):
        warnings.append("Detached positive signal was excluded; inspect dim necks and nearby fibers.")
    if np.issubdtype(data.image.dtype, np.integer):
        frac = float(np.mean(raw[g.corridor] == np.iinfo(data.image.dtype).max))
        if frac > 0.001:
            warnings.append(f"ROI: {frac:.1%} of pixels reach the storage maximum; inspect acquisition saturation.")
    area = float(spine.sum() * data.sx * data.sy)
    summary = {"metric": "2D fluorescence-positive protrusion area per centerline length",
               "spine_area_um2": area, "dendrite_length_um": g.length,
               "area_per_length_um2_per_um": area / g.length,
               "spine_pixels": int(spine.sum()), "pixel_area_um2": data.sx * data.sy,
               "shaft_positive_area_um2": float(shaft.sum() * data.sx * data.sy),
               "accepted_components_not_spine_count": int(ndi.label(spine, STRUCTURE)[1]) if brush_strokes else int(keep.sum()),
               "sx_um_per_pixel": data.sx, "sy_um_per_pixel": data.sy,
               "signal_threshold": p.threshold, "shaft_threshold": p.shaft_threshold or p.threshold}
    if brush_strokes:
        summary.update(automatic_spine_area_um2=float(automatic_spines.sum()*data.sx*data.sy),
                       manual_added_area_um2=float(manual_added.sum()*data.sx*data.sy),
                       manual_removed_area_um2=float(manual_removed.sum()*data.sx*data.sy),
                       manual_brush_strokes=len(brush_strokes))
    recorded_strokes = [{"mode": s["mode"], "radius_um": float(s["radius_um"]),
                         "points_xy_source_px": np.asarray(s["points_xy_source_px"], float).tolist()}
                        for s in brush_strokes or []]
    return Result(g, raw, a, positive, attached, band, shaft, candidate, spine,
                  exclusion, up, down, components, warnings, summary,
                  automatic_spines, manual_added, manual_removed, recorded_strokes)


def suggest_threshold(corrected, g=None):
    a = corrected if g is None else corrected[g.bbox[1]:g.bbox[3], g.bbox[0]:g.bbox[2]][g.corridor]
    if np.ptp(a) == 0:
        raise ValueError("Otsu requires nonconstant image/ROI intensities.")
    return max(float(threshold_otsu(a)), float(np.finfo(np.float32).eps))


def write_tiff(path, array, sx, sy):
    tifffile.imwrite(path, np.asarray(array), imagej=True, resolution=(1 / sx, 1 / sy),
                     metadata={"axes": "YX", "unit": "um"})


def image_digest(image):
    return hashlib.sha256(np.ascontiguousarray(image).tobytes()).hexdigest()


def source_labels(data, result):
    labels = np.zeros(data.image.shape, np.uint8)
    x0, y0, x1, y1 = result.geometry.bbox
    crop = labels[y0:y1, x0:x1]
    crop[result.shaft] = 1
    crop[result.spines] = 2
    return labels


def source_annotation(data, result):
    from PIL import Image, ImageDraw
    high = max(1., float(np.percentile(data.image, 99.7)))
    gray = np.clip(data.image.astype(np.float32) * (255 / high), 0, 255).astype(np.uint8)
    rgb = np.zeros((*gray.shape, 3), np.uint8)
    rgb[:, :, 0] = gray
    x0, y0, x1, y1 = result.geometry.bbox
    crop = rgb[y0:y1, x0:x1]
    for mask, color in ((result.shaft, (26, 255, 51)), (result.spines, (255, 0, 255))):
        edge = mask & ~ndi.binary_erosion(mask, structure=STRUCTURE, border_value=0)
        crop[edge] = color
    image = Image.fromarray(rgb)
    draw = ImageDraw.Draw(image)
    bx0, by0, bx1, by1 = result.geometry.selection_box or result.geometry.bbox
    draw.rectangle((bx0, by0, bx1 - 1, by1 - 1), outline=(255, 255, 0), width=1)
    components, _ = ndi.label(result.spines, STRUCTURE)
    for location in ndi.find_objects(components):
        if location is not None:
            rows, cols = location
            draw.rectangle((max(x0, x0 + cols.start - 2), max(y0, y0 + rows.start - 2),
                            min(x1 - 1, x0 + cols.stop + 1), min(y1 - 1, y0 + rows.stop + 1)),
                           outline=(255, 0, 255), width=1)
    path = [tuple(point) for point in result.geometry.path]
    draw.line(path, fill=(0, 255, 255), width=1)
    for x, y in (path[0], path[-1]):
        draw.ellipse((x - 2, y - 2, x + 2, y + 2), outline=(255, 255, 0), width=1)
    return np.asarray(image)


def save_result(folder, data, points, p, result, exclusions=None,
                upper_points=None, lower_points=None):
    from PIL import Image
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    folder = Path(folder)
    folder.mkdir(parents=True, exist_ok=False)
    g = result.geometry
    labels = source_labels(data, result)
    annotated = source_annotation(data, result)
    write_tiff(folder / "source_raw.tif", data.image, data.sx, data.sy)
    write_tiff(folder / "source_labels.tif", labels, data.sx, data.sy)
    write_tiff(folder / "source_shaft_mask.tif", (labels == 1).astype(np.uint8) * 255, data.sx, data.sy)
    write_tiff(folder / "source_spine_mask.tif", (labels == 2).astype(np.uint8) * 255, data.sx, data.sy)
    for name, mask in (("source_auto_spine_mask", result.automatic_spines if result.automatic_spines is not None else result.spines),
                       ("source_manual_added_mask", result.manual_added), ("source_manual_removed_mask", result.manual_removed)):
        full = np.zeros(data.image.shape, np.uint8)
        if mask is not None:
            x0, y0, x1, y1 = g.bbox
            full[y0:y1, x0:x1] = mask.astype(np.uint8) * 255
        write_tiff(folder / (name + ".tif"), full, data.sx, data.sy)
    tifffile.imwrite(folder / "source_annotated.tif", annotated, imagej=True, photometric="rgb",
                     resolution=(1 / data.sx, 1 / data.sy), metadata={"axes": "YXS", "unit": "um"})
    Image.fromarray(annotated).save(folder / "source_annotated.png")
    masks = {"03_corridor": g.corridor, "04_positive": result.positive,
             "05_attached": result.attached, "06_shaft_domain": result.shaft_domain,
             "07_shaft_positive": result.shaft, "08_spine_mask": result.spines,
             "09_other_positive": result.positive & ~result.shaft & ~result.spines,
             "10_excluded": result.excluded, "11_candidates_before_filters": result.candidates}
    write_tiff(folder / "01_raw_crop.tif", result.raw, data.sx, data.sy)
    write_tiff(folder / "02_corrected_crop.tif", result.corrected, data.sx, data.sy)
    for name, mask in masks.items():
        write_tiff(folder / (name + ".tif"), mask.astype(np.uint8) * 255, data.sx, data.sy)
    ds = g.length / len(g.s_axis)
    dt = g.t_axis[1] - g.t_axis[0]
    write_tiff(folder / "straightened_intensity_VIEW_ONLY.tif",
               g.straighten(data.image.astype(np.float32)), ds, dt)
    write_tiff(folder / "straightened_spines_VIEW_ONLY.tif",
               g.straighten(result.spines.astype(np.uint8) * 255, order=0, crop=True), ds, dt)
    session = {"program_version": VERSION, "created_utc": datetime.now(timezone.utc).isoformat(),
               "source": data.source, "source_selection": data.selection,
               "source_shape_yx": list(data.image.shape), "source_dtype": str(data.image.dtype),
               "projection_sha256": image_digest(data.image),
               "sx_um_per_pixel": data.sx, "sy_um_per_pixel": data.sy,
               "settings": asdict(p), "trace_xy_source_px": np.asarray(points).tolist(),
               "selected_centerline_xy_source_px": g.path.tolist(),
               "crop_bbox_source_px_xyxy": list(g.bbox),
               "selection_box_source_px_xyxy": g.selection_box,
               "exclusion_polygons_xy_source_px": exclusions or [],
               "upper_boundary_points_s_t_um": upper_points or [],
               "lower_boundary_points_s_t_um": lower_points or [],
               "brush_strokes": result.brush_strokes,
               "brush_rule": "last stroke wins; add requires threshold-positive signal inside corridor, outside shaft/exclusions; manual edits override automatic connectivity/area/reach filters",
               "components_csv_scope": "automatic candidates before brush edits; final_components.csv describes the final mask",
               "summary": result.summary, "warnings": result.warnings,
               "source_annotation": {"shape_yx": list(data.image.shape),
                                     "labels": {"0": "unlabeled", "1": "shaft positive", "2": "spine candidate"},
                                     "colors": {"shaft": "green", "spines": "magenta", "centerline": "cyan", "ends": "yellow"},
                                     "scope": "selected segment in the analyzed channel/projection; original pixel coordinates",
                                     "display": "red fluorescence with green/magenta contours; quantitative pixels remain in source_raw.tif"},
               "area_definition": "count original source pixels; inclusive threshold; 8-connectivity; no hole filling",
               "runtime": {"python": platform.python_version(), "numpy": np.__version__,
                           "scipy": __import__("scipy").__version__,
                           "skimage": __import__("skimage").__version__,
                           "czifile": getattr(sys.modules.get("czifile"), "__version__", "not loaded")}}
    (folder / "session.json").write_text(json.dumps(session, ensure_ascii=False, indent=2), encoding="utf-8")
    if data.metadata_xml:
        (folder / "czi_metadata.xml").write_text(data.metadata_xml, encoding="utf-8")
    with (folder / "results.csv").open("w", newline="", encoding="utf-8-sig") as f:
        fields = ["source", *result.summary]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow({"source": data.source, **result.summary})
    with (folder / "components.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.DictWriter(f, fieldnames=["component_id", "area_um2", "centroid_x_source_px",
                                             "centroid_y_source_px", "status"])
        writer.writeheader()
        writer.writerows(result.components)
    final_labels, _ = ndi.label(result.spines, STRUCTURE)
    with (folder / "final_components.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["component_id", "area_um2", "centroid_x_source_px", "centroid_y_source_px"])
        for i, location in enumerate(ndi.find_objects(final_labels), 1):
            if location is not None:
                yy, xx = np.nonzero(final_labels[location] == i)
                writer.writerow([i, len(xx)*data.sx*data.sy, float(xx.mean()+location[1].start+g.bbox[0]),
                                 float(yy.mean()+location[0].start+g.bbox[1])])
    with (folder / "shaft_boundaries.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["s_um", "upper_t_um", "lower_t_um"])
        writer.writerows(zip(g.s_axis, -result.upper, result.lower))
    with (folder / "centerline.csv").open("w", newline="", encoding="utf-8-sig") as f:
        writer = csv.writer(f)
        writer.writerow(["x_source_px", "y_source_px", "x_crop_px", "y_crop_px"])
        writer.writerows((x, y, x - g.bbox[0], y - g.bbox[1]) for x, y in g.path)
    np.savetxt(folder / "centerline_xy.txt", g.path, fmt="%.10f", delimiter="\t")
    fig = Figure(figsize=(12, 7), layout="constrained")
    FigureCanvasAgg(fig)
    a1, a2 = fig.subplots(2, 1)
    bright = np.percentile(result.raw, 99.8)
    a1.imshow(result.raw, cmap="gray", vmin=0, vmax=max(1, bright))
    rgba = overlay(result)
    a1.imshow(rgba)
    a1.plot(g.path[:, 0] - g.bbox[0], g.path[:, 1] - g.bbox[1], color="cyan", lw=1)
    a1.set_title("Source pixels: shaft green; protrusions magenta; other positive amber")
    plot_horizontal_view(a2, data, result)
    a2.set(xlabel="Centerline distance (um)", ylabel="Normal offset (um)",
           title="Horizontal view: green shaft / magenta spines\n"
                 f"Area {result.summary['spine_area_um2']:.4f} um2 / {g.length:.4f} um = "
                 f"{result.summary['area_per_length_um2_per_um']:.5f} um2/um; VIEW ONLY")
    fig.savefig(folder / "review_overlay.png", dpi=160)
    fig.clear()
    fig.set_size_inches(9, 5)
    horizontal_ax = fig.subplots()
    for filled, name in ((True, "horizontal_filled_VIEW_ONLY.png"),
                         (False, "horizontal_outlines_VIEW_ONLY.png")):
        horizontal_ax.clear()
        plot_horizontal_view(horizontal_ax, data, result, filled=filled)
        fig.savefig(folder / name, dpi=160)
    return folder


def plot_horizontal_view(ax, data, result, high=None, filled=True):
    g = result.geometry
    high = max(1., float(np.percentile(data.image, 99.7)) if high is None else high)
    rgb = np.zeros((*g.grid_x.shape, 3), np.float32)
    rgb[:, :, 0] = np.clip(g.straighten(data.image.astype(np.float32)) / high, 0, 1)
    extent = (0, g.length, g.t_axis[-1], g.t_axis[0])
    ax.imshow(rgb, extent=extent, aspect="equal", interpolation="nearest", origin="upper")
    rgba = np.zeros((*g.grid_x.shape, 4), np.float32)
    # Nearest-neighbor sampling preserves mask classes; areas still use source pixels.
    for mask, color in ((result.shaft, (.1, 1, .2, .7)), (result.spines, (1, 0, 1, .85))):
        straight = g.straighten(mask.astype(np.uint8), order=0, crop=True).astype(bool)
        if filled:
            rgba[straight] = color
        elif straight.any():
            # Outline the same accepted pixels without covering fluorescence inside them.
            rows, cols = straight.shape
            xs = (np.arange(-1, cols + 1) + .5) * g.length / cols
            ys = g.t_axis[0] + (np.arange(-1, rows + 1) + .5) * np.ptp(g.t_axis) / rows
            ax.contour(xs, ys, np.pad(straight, 1), levels=[.5], colors=[color[:3]], linewidths=.8)
    if filled:
        ax.imshow(rgba, extent=extent, aspect="equal", interpolation="nearest", origin="upper")
    ax.plot(g.s_axis, -result.upper, color="lime", lw=1)
    ax.plot(g.s_axis, result.lower, color="lime", lw=1)
    style = "filled" if filled else "outlines only"
    ax.set(title=f"Horizontal view: {style} / green shaft / magenta spines",
           xlabel="Centerline distance (um)", ylabel="Normal offset (um)",
           xlim=(0, g.length), ylim=(g.t_axis[-1], g.t_axis[0]))


def overlay(r: Result):
    out = np.zeros((*r.raw.shape, 4), np.float32)
    out[r.positive & ~r.shaft & ~r.spines] = [1, .65, 0, .55]
    out[r.shaft] = [.1, 1, .2, .5]
    out[r.spines] = [1, 0, 1, .8]
    out[r.excluded] = [.2, .6, 1, .4]
    return out


def demo():
    """A deterministic synthetic dendrite, unrelated to biological validation."""
    rng = np.random.default_rng(12)
    yy, xx = np.mgrid[:180, :380]
    trunk = (np.abs(yy - 90) <= 6) & (xx >= 20) & (xx <= 360)
    spine = np.zeros(trunk.shape, bool)
    for x, direction in ((65, -1), (110, 1), (160, -1), (210, 1), (260, -1), (305, 1)):
        spine |= (np.abs(xx - x) <= 2) & (direction * (yy - 90) > 6) & (direction * (yy - 90) < 16)
        spine |= (xx - x) ** 2 + (yy - (90 + direction * 17)) ** 2 <= 5 ** 2
    a = rng.normal(100, 15, trunk.shape)
    a[trunk] += 6000
    a[spine] += 4000
    data = ImageData(np.clip(a, 0, 65535).astype(np.uint16), .1, .1, "SYNTHETIC_DEMO",
                     {"calibration_source": "synthetic 0.1 um/pixel; not measured CZI"})
    points = [[40., 90.], [140., 90.], [240., 90.], [340., 90.]]
    p = Settings(length_um=25, start_offset_um=0, threshold=1000)
    return data, points, p
