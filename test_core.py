"""Small deterministic checks for geometry, calibration, segmentation and export."""
from pathlib import Path
import csv
import json
import tempfile
import unittest

import numpy as np
import tifffile

from core import (ImageData, Settings, analyze, demo, geometry, preprocess,
                  save_result, scaling_from_xml, trim_path, source_labels, source_annotation)


class MeasurementTests(unittest.TestCase):
    def simple(self):
        a = np.zeros((100, 160), np.uint16)
        a[48:53, 20:130] = 5000
        a[40:48, 45:49] = 3000  # 32 pixels, connected protrusion
        a[53:61, 75:79] = 3000  # 32 pixels
        a[30:33, 100:103] = 3000  # detached neighbour
        d = ImageData(a, .1, .2, "synthetic_test")
        p = Settings(length_um=10, start_offset_um=0, threshold=1000,
                     half_width_um=5, shaft_halfwidth_um=.5, min_area_um2=0,
                     max_area_um2=0, spine_reach_um=3)
        points = [[20, 50], [120, 50]]
        g = geometry(d, points, p)
        r = analyze(d, preprocess(d, p), g, p)
        return d, p, points, r

    def test_exact_length_offset_and_anisotropic_scale(self):
        out = trim_path([[0, 0], [30, 0], [30, 40]], .2, .1, 7, 1)
        np.testing.assert_allclose(out, [[5, 0], [30, 0], [30, 20]])
        physical = out * [.2, .1]
        self.assertAlmostEqual(np.linalg.norm(np.diff(physical, axis=0), axis=1).sum(), 7)
        with self.assertRaises(ValueError):
            trim_path([[0, 0], [1, 0]], .1, .1, 20)

    def test_counts_source_pixels_not_resampled_pixels(self):
        d, p, points, r = self.simple()
        self.assertEqual(r.summary["spine_pixels"], 64)
        self.assertAlmostEqual(r.summary["spine_area_um2"], 1.28)
        self.assertAlmostEqual(r.summary["area_per_length_um2_per_um"], .128)
        self.assertEqual(r.summary["accepted_components_not_spine_count"], 2)
        self.assertTrue(np.any(r.positive & ~r.attached))
        self.assertFalse(np.any(r.shaft & r.spines))

    def test_threshold_inclusive_and_noise_excluded(self):
        d, p, points, r = self.simple()
        p.threshold = 3000
        r2 = analyze(d, preprocess(d, p), r.geometry, p)
        self.assertEqual(r2.summary["spine_pixels"], 64)
        p.threshold = 3001
        r3 = analyze(d, preprocess(d, p), r.geometry, p)
        self.assertEqual(r3.summary["spine_pixels"], 0)

    def test_independent_shaft_threshold_and_component_coordinates(self):
        d, p, points, _ = self.simple()
        p.shaft_halfwidth_um = 0
        p.shaft_threshold = 4000
        g = geometry(d, points, p)
        low = analyze(d, preprocess(d, p), g, p)
        p.threshold = 3001
        high = analyze(d, preprocess(d, p), g, p)
        np.testing.assert_array_equal(low.upper, high.upper)
        np.testing.assert_array_equal(low.lower, high.lower)
        self.assertEqual(low.summary["spine_pixels"], 64)
        self.assertEqual(high.summary["spine_pixels"], 0)
        centers = [(c["centroid_x_source_px"], c["centroid_y_source_px"]) for c in low.components]
        np.testing.assert_allclose(centers, [(46.5, 43.5), (76.5, 56.5)])
        self.assertEqual(low.summary["shaft_threshold"], 4000)

    def test_rotation_preserves_source_pixel_area(self):
        d, p, points, r = self.simple()
        rotated = ImageData(np.rot90(d.image).copy(), d.sy, d.sx, "rotated")
        xy = [[y, d.image.shape[1] - 1 - x] for x, y in points]
        g = geometry(rotated, xy, p)
        rr = analyze(rotated, preprocess(rotated, p), g, p)
        self.assertEqual(r.summary["spine_pixels"], rr.summary["spine_pixels"])
        self.assertAlmostEqual(r.summary["area_per_length_um2_per_um"], rr.summary["area_per_length_um2_per_um"])

    def test_long_branch_rejected_whole_not_clipped(self):
        d, p, points, r = self.simple()
        d.image[27:48, 45:49] = 3000
        p.spine_reach_um = 2
        r2 = analyze(d, preprocess(d, p), r.geometry, p)
        self.assertEqual(r2.summary["spine_pixels"], 32)
        self.assertIn("exceeds_reach", [c["status"] for c in r2.components])

    def test_size_filter_and_polygon(self):
        d, p, points, r = self.simple()
        poly = [[[43, 38], [51, 38], [51, 48], [43, 48]]]
        r2 = analyze(d, preprocess(d, p), r.geometry, p, poly)
        self.assertEqual(r2.summary["spine_pixels"], 32)
        p.max_area_um2 = .6
        r3 = analyze(d, preprocess(d, p), r.geometry, p)
        self.assertEqual(r3.summary["spine_pixels"], 0)

    def test_auto_shaft_and_manual_edges(self):
        d, points, p = demo()
        g = geometry(d, points, p)
        r = analyze(d, preprocess(d, p), g, p)
        self.assertGreater(r.summary["spine_area_um2"], 1)
        np.testing.assert_allclose(np.median(r.upper), .65, atol=.05)
        up, down = [[0, -.65], [25, -.65]], [[0, .65], [25, .65]]
        manual = analyze(d, preprocess(d, p), g, p, [], up, down)
        self.assertEqual(r.summary["spine_pixels"], manual.summary["spine_pixels"])

    def test_xml_scaling_namespace_and_missing(self):
        xml = '<ImageDocument xmlns="test"><Metadata><Scaling><Items>' \
              '<Distance Id="X"><Value>1e-7</Value></Distance>' \
              '<Distance Id="Y"><Value>2e-7</Value></Distance>' \
              '</Items></Scaling></Metadata></ImageDocument>'
        np.testing.assert_allclose(scaling_from_xml(xml), [.1, .2])
        with self.assertRaises(ValueError):
            scaling_from_xml("<ImageDocument/>")

    def test_export_calibration_and_mask_area(self):
        d, p, points, r = self.simple()
        test_root = Path(__file__).resolve().parent
        with tempfile.TemporaryDirectory(prefix="test-export-", dir=test_root) as tmp:
            self.assertTrue(Path(tmp).resolve().is_relative_to(test_root))
            out = save_result(Path(tmp) / "result", d, points, p, r)
            with tifffile.TiffFile(out / "08_spine_mask.tif") as tif:
                mask = tif.asarray()
                self.assertEqual(np.count_nonzero(mask), 64)
                tx = tif.pages[0].tags["XResolution"].value
                ty = tif.pages[0].tags["YResolution"].value
                self.assertAlmostEqual(tx[0] / tx[1], 10)
                self.assertAlmostEqual(ty[0] / ty[1], 5)
                self.assertEqual(tif.imagej_metadata["unit"], "um")
            np.testing.assert_array_equal(tifffile.imread(out / "source_raw.tif"), d.image)
            labels = tifffile.imread(out / "source_labels.tif")
            x0, y0, x1, y1 = r.geometry.bbox
            expected = np.zeros(d.image.shape, np.uint8)
            expected[y0:y1, x0:x1][r.shaft] = 1
            expected[y0:y1, x0:x1][r.spines] = 2
            np.testing.assert_array_equal(labels, expected)
            self.assertEqual(np.count_nonzero(labels == 2), 64)
            for value, name in ((1, "shaft"), (2, "spine")):
                mask = tifffile.imread(out / f"source_{name}_mask.tif")
                np.testing.assert_array_equal(mask, (expected == value).astype(np.uint8) * 255)
            from PIL import Image
            annotated = np.asarray(Image.open(out / "source_annotated.png"))
            self.assertEqual(annotated.shape, (*d.image.shape, 3))
            np.testing.assert_array_equal(annotated, tifffile.imread(out / "source_annotated.tif"))
            with tifffile.TiffFile(out / "source_annotated.tif") as tif:
                self.assertEqual(tif.imagej_metadata["unit"], "um")
                self.assertEqual(tif.pages[0].tags["PhotometricInterpretation"].value, 2)
            np.testing.assert_allclose(np.loadtxt(out / "centerline_xy.txt"), r.geometry.path)
            self.assertFalse(list(out.glob("*.ijm")))
            with self.assertRaises(FileExistsError):
                save_result(out, d, points, p, r)

    def test_source_annotation_alignment_after_rotation(self):
        d, p, points, r = self.simple()
        original = d.image.copy()
        rotated = ImageData(np.rot90(d.image).copy(), d.sy, d.sx, "rotated")
        xy = [[y, d.image.shape[1] - 1 - x] for x, y in points]
        rr = analyze(rotated, preprocess(rotated, p), geometry(rotated, xy, p), p)
        np.testing.assert_array_equal(source_labels(rotated, rr), np.rot90(source_labels(d, r)))
        image = source_annotation(d, r)
        outside = np.ones(d.image.shape, bool)
        x0, y0, x1, y1 = r.geometry.bbox
        outside[y0:y1, x0:x1] = False
        self.assertFalse(np.any(image[outside, 1:]))
        np.testing.assert_array_equal(d.image, original)

    def test_rectangle_limits_source_pixel_measurement(self):
        d, p, points, r = self.simple()
        box = (15, 45, 125, 56)
        g = geometry(d, points, p, box)
        selected = analyze(d, preprocess(d, p), g, p)
        self.assertEqual(selected.summary["spine_pixels"], 24)
        labels = source_labels(d, selected)
        self.assertFalse(labels[:45].any())
        self.assertFalse(labels[56:].any())
        self.assertEqual(g.selection_box, box)
        with self.assertRaises(ValueError):
            geometry(d, points, p, (30, 45, 125, 56))

    def test_reject_invalid_parameters_and_self_crossing(self):
        p = Settings(threshold=0)
        with self.assertRaises(ValueError):
            p.validate()
        with self.assertRaises(ValueError):
            trim_path([[10, 10], [90, 90], [10, 90], [90, 10]], .1, .1, 25)

    def test_brush_calibration_addition_and_guards(self):
        d, p, points, r = self.simple()
        original = d.image.copy()
        remove = {"mode":"remove", "radius_um":.21, "points_xy_source_px":[[46,44]]}
        edited = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[remove])
        self.assertEqual(edited.summary["spine_pixels"], 58)
        self.assertEqual(edited.manual_removed.sum(), 6)
        np.testing.assert_array_equal(edited.shaft, r.shaft)
        add = {"mode":"add", "radius_um":.3, "points_xy_source_px":[[101,31]]}
        edited = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[add])
        self.assertEqual(edited.summary["spine_pixels"], 73)
        self.assertEqual(edited.manual_added.sum(), 9)
        for center in ([90,50], [30,30], [5,5]):
            paint = dict(add, points_xy_source_px=[center])
            guarded = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[paint])
            np.testing.assert_array_equal(guarded.spines, r.spines)
        polygon = [[[99,29],[104,29],[104,34],[99,34]]]
        excluded = analyze(d, preprocess(d,p), r.geometry, p, exclusions=polygon, brush_strokes=[add])
        self.assertEqual(excluded.manual_added.sum(), 0)
        p.threshold = 3001
        dim = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[add])
        self.assertEqual(dim.summary["spine_pixels"], 0)
        np.testing.assert_array_equal(d.image, original)
        with self.assertRaises(ValueError):
            analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[dict(add, radius_um=0)])

    def test_brush_continuous_strokes_last_action_and_export_replay(self):
        d, p, points, r = self.simple()
        remove = {"mode":"remove", "radius_um":.06, "points_xy_source_px":[[46,40],[46,47]]}
        erased = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=[remove])
        self.assertEqual(erased.summary["spine_pixels"], 56)
        restored = analyze(d, preprocess(d,p), r.geometry, p,
                           brush_strokes=[remove, dict(remove, mode="add")])
        np.testing.assert_array_equal(restored.spines, r.spines)
        strokes = [remove, {"mode":"add", "radius_um":.3, "points_xy_source_px":[[101,31]]}]
        edited = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=strokes)
        self.assertAlmostEqual(edited.summary["spine_area_um2"], 1.3)
        with tempfile.TemporaryDirectory(prefix="test-brush-", dir=Path(__file__).resolve().parent) as tmp:
            out = save_result(Path(tmp)/"result", d, points, p, edited)
            saved = json.loads((out/"session.json").read_text(encoding="utf-8"))
            self.assertEqual(saved["brush_strokes"], strokes)
            replay = analyze(d, preprocess(d,p), r.geometry, p, brush_strokes=saved["brush_strokes"])
            np.testing.assert_array_equal(replay.spines, edited.spines)
            for name, count in (("source_auto_spine_mask",64), ("source_spine_mask",65),
                                ("source_manual_added_mask",9), ("source_manual_removed_mask",8)):
                mask = tifffile.imread(out/(name+".tif"))
                self.assertEqual(mask.shape, d.image.shape)
                self.assertEqual(np.count_nonzero(mask), count)
            with (out/"final_components.csv").open(encoding="utf-8-sig") as f:
                rows = list(csv.DictReader(f))
            self.assertAlmostEqual(sum(float(row["area_um2"]) for row in rows), 1.3)
            self.assertEqual(len(rows), edited.summary["accepted_components_not_spine_count"])


if __name__ == "__main__":
    unittest.main(verbosity=2)
