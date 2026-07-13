"""Tab 5 "Simulation cell": preview and configure the elongated cell the
engine builds for a run. The whole periodic box is drawn: the crystal slab
tiled from the saved unit cell, the field zones as background bands, and
underneath, sharing the same x axis, the damping profile lambda(x) and the
external potential u_ext(x), so it is immediately visible which molecules
fall into which field while the parameters are tuned.

The tab edits the cell-shaping parameters (uc_in_x, uc_in_y, free_space,
temperature_in_transition_zone, lambdam) and stores them in the project.
The temperature and u_m loops belong to the Run tab; the local T and u_m
spins here only drive the preview curves.
"""

import math

from PySide6.QtCore import Qt, Signal
from PySide6.QtGui import QBrush, QColor, QPainterPath, QPen, QTransform
from PySide6.QtWidgets import (QDoubleSpinBox, QFrame, QGraphicsItem,
                               QGraphicsPathItem, QGraphicsRectItem,
                               QGraphicsSimpleTextItem, QGraphicsView,
                               QGridLayout, QGroupBox, QHBoxLayout, QLabel,
                               QMessageBox, QPushButton, QSpinBox, QSplitter,
                               QVBoxLayout, QWidget)

from .. import theme
from ..glyph import fallback_glyph, model_glyph
from ..gridview import GridView
from ..project import Project
from ..simcell import (ENGINE_BUFFER, KSI_CRYSTAL, KSI_EDGE, KSI_GAS,
                       KSI_PLATEAU, KSI_PLATEAU_END, SimCellSpec,
                       lambda_profile, tile_unit_cell, u_profile)

N_A = 6.02214076e23
# past this many glyph disks the preview falls back to one dot per molecule
MAX_DISKS = 40000


class SimCellCanvas(GridView):
    """Read-only, pannable view of the elongated cell with the two field
    profiles drawn in the same scene right below the box, so they stay
    aligned with the molecules at any zoom. The angstrom grid is off: at
    this scale it only drowns the molecules."""

    MIN_SCALE = 0.2
    SHOW_GRID = False

    def __init__(self, parent=None):
        super().__init__(parent)
        self.setDragMode(QGraphicsView.ScrollHandDrag)
        self.setHorizontalScrollBarPolicy(Qt.ScrollBarAlwaysOff)
        self.setVerticalScrollBarPolicy(Qt.ScrollBarAlwaysOff)

    def clear_data(self) -> None:
        self.scene().clear()

    def set_data(self, layout, placements, temp: float, um: float,
                 glyph) -> None:
        scene = self.scene()
        scene.clear()
        lx, ly = layout.lx, layout.ly
        spec = layout.spec

        # strip geometry below the cell (scene y is up, strips are negative)
        h = max(0.22 * ly, 10.0)
        gap = max(0.08 * ly, 4.0)
        lam_top, lam_bot = -gap, -gap - h
        u_top, u_bot = lam_bot - gap, lam_bot - gap - h

        # field zones as vertical bands through the cell and both strips
        bands = ((0.0, KSI_CRYSTAL, theme.ACCENT, 34),
                 (KSI_CRYSTAL, KSI_PLATEAU, theme.GOLD, 16),
                 (KSI_PLATEAU, KSI_PLATEAU_END, theme.GOLD, 36),
                 (KSI_PLATEAU_END, KSI_GAS, theme.TEXT_DIM, 14),
                 (KSI_GAS, KSI_EDGE, theme.TEXT_DIM, 26))
        for k0, k1, color, alpha in bands:
            col = QColor(color)
            col.setAlpha(alpha)
            for side in (-1, 1):
                xa = lx / 2.0 + side * k0 * lx / 32.0
                xb = lx / 2.0 + side * k1 * lx / 32.0
                band = QGraphicsRectItem(min(xa, xb), u_bot,
                                         abs(xb - xa), ly - u_bot)
                band.setPen(QPen(Qt.NoPen))
                band.setBrush(QBrush(col))
                band.setZValue(-2)
                scene.addItem(band)
        self._text("crystal", lx / 2.0, ly + 0.5, center=True, above=True)
        self._text("transition", layout.zone_x(
            (KSI_PLATEAU + KSI_PLATEAU_END) / 2.0)[1], ly + 0.5,
            center=True, above=True)
        self._text("gas", layout.zone_x(
            (KSI_GAS + KSI_EDGE) / 2.0)[1], ly + 0.5, center=True, above=True)

        # the cell box and the molecules
        cell = QGraphicsRectItem(0, 0, lx, ly)
        cell.setPen(QPen(QColor(theme.ACCENT), 0))
        cell.setZValue(-1)
        scene.addItem(cell)
        glyph = list(glyph) if glyph else fallback_glyph()
        if len(placements) * len(glyph) > MAX_DISKS:
            glyph = fallback_glyph()
        for px, py, phi in placements:
            rad = math.radians(phi)
            c, s = math.cos(rad), math.sin(rad)
            for gx, gy, color, r in glyph:
                x = px + gx * c - gy * s
                y = py + gx * s + gy * c
                dot = scene.addEllipse(x - r, y - r, 2 * r, 2 * r,
                                       QPen(Qt.NoPen), QBrush(QColor(color)))
                dot.setZValue(1)

        # the two profiles, sampled over the full cell width
        lam = lambda x: lambda_profile(x, lx, temp, spec)
        lam_y = lambda v: lam_bot + v * h        # lambda is 0..1 by design
        self._strip(lx, lam_top, lam_bot)
        self._curve(lx, lambda x: lam_y(lam(x)), theme.ACCENT_HOVER)
        self._text("λ(x)", 1, lam_top)
        self._text("1", lx, lam_top, vcenter=True)
        self._text("0", lx, lam_bot, vcenter=True)

        u_lo, u_hi = min(0.0, um), max(0.0, um)
        span = (u_hi - u_lo) or 1.0
        u_y = lambda v: u_bot + (v - u_lo) / span * h
        self._strip(lx, u_top, u_bot)
        self._curve(lx, lambda x: u_y(u_profile(x, lx, um)), theme.GOLD)
        self._text("u_ext(x)", 1, u_top)
        self._text("0", lx, u_y(0.0), vcenter=True)
        if um:
            self._text(f"{um / 1000.0:g} kJ/mol", lx, u_y(um), vcenter=True)

        self.fit_points([(-2, u_bot - 2), (lx + 2, ly + 2)], pad=4)
        # second pass: the zone captions above the box are device-sized,
        # so convert their pixel height into scene units at the found zoom
        scale = abs(self.transform().m11())
        self.fit_points([(-2, u_bot - 2), (lx + 2, ly + 2 + 24.0 / scale)],
                        pad=4)

    # -- drawing helpers -----------------------------------------------------

    def _strip(self, lx: float, top: float, bot: float) -> None:
        frame = QGraphicsRectItem(0, bot, lx, top - bot)
        pen = QPen(QColor(theme.BORDER), 0)
        pen.setCosmetic(True)
        frame.setPen(pen)
        self.scene().addItem(frame)

    def _curve(self, lx: float, f, color: str, samples: int = 512) -> None:
        path = QPainterPath()
        for i in range(samples + 1):
            x = lx * i / samples
            if i == 0:
                path.moveTo(x, f(x))
            else:
                path.lineTo(x, f(x))
        item = QGraphicsPathItem(path)
        pen = QPen(QColor(color))
        pen.setWidthF(1.6)
        pen.setCosmetic(True)
        item.setPen(pen)
        item.setZValue(2)
        self.scene().addItem(item)

    def _text(self, s: str, x: float, y: float, center: bool = False,
              vcenter: bool = False, above: bool = False) -> None:
        t = QGraphicsSimpleTextItem(s)
        t.setBrush(QColor(theme.TEXT_DIM))
        # keep labels upright and screen-sized regardless of the y-flip/zoom;
        # the item transform below is applied in device pixels
        t.setFlag(QGraphicsItem.ItemIgnoresTransformations)
        t.setPos(x, y)
        rect = t.boundingRect()
        dx = -rect.width() / 2.0 if center else (4.0 if vcenter else 0.0)
        dy = -rect.height() if above else (-rect.height() / 2.0 if vcenter
                                           else 0.0)
        if dx or dy:
            t.setTransform(QTransform.fromTranslate(dx, dy))
        self.scene().addItem(t)


class SimulationCellTab(QWidget):
    statusMessage = Signal(str)
    projectCellChanged = Signal()

    def __init__(self, project: Project, parent=None):
        super().__init__(parent)
        self.project = project
        self._dirty = False

        root = QVBoxLayout(self)
        root.setContentsMargins(12, 12, 12, 12)

        bar = QHBoxLayout()
        fit = QPushButton("Fit view")
        fit.clicked.connect(self._recompute)
        bar.addWidget(fit)
        bar.addStretch(1)
        attach = QPushButton("Use in project")
        attach.setProperty("primary", True)
        attach.clicked.connect(self.use_in_project)
        bar.addWidget(attach)
        root.addLayout(bar)

        # the cell is long in x: the canvas takes the full width, the
        # parameters sit in a strip below it
        splitter = QSplitter(Qt.Vertical)
        root.addWidget(splitter, 1)
        self.canvas = SimCellCanvas()
        self.canvas.cursorMoved.connect(
            lambda x, y: self.statusMessage.emit(f"x = {x:.2f} Å,  y = {y:.2f} Å"))
        splitter.addWidget(self.canvas)
        splitter.addWidget(self._build_bottom())
        splitter.setStretchFactor(0, 1)
        splitter.setStretchFactor(1, 0)
        splitter.setSizes([620, 240])

    # -- construction --------------------------------------------------------

    def _spin(self, lo, hi, val, step, dec=0, suffix="") -> QDoubleSpinBox:
        s = QDoubleSpinBox() if dec else QSpinBox()
        s.setRange(lo, hi)
        if dec:
            s.setDecimals(dec)
        s.setSingleStep(step)
        s.setValue(val)
        if suffix:
            s.setSuffix(suffix)
        s.valueChanged.connect(self._on_changed)
        return s

    def _build_bottom(self) -> QWidget:
        panel = QWidget()
        row = QHBoxLayout(panel)
        row.setContentsMargins(0, 8, 0, 0)
        row.setSpacing(16)

        left = QVBoxLayout()
        self.info = QLabel()
        self.info.setWordWrap(True)
        self.info.setTextFormat(Qt.RichText)
        left.addWidget(self.info)
        left.addStretch(1)
        caption = QLabel("Simulation cell in project")
        caption.setProperty("dim", True)
        left.addWidget(caption)
        frame = QFrame()
        frame.setProperty("panel", True)
        fl = QVBoxLayout(frame)
        self.slot_label = QLabel()
        self.slot_label.setWordWrap(True)
        fl.addWidget(self.slot_label)
        left.addWidget(frame)
        row.addLayout(left, 1)

        self.cell_box = QGroupBox("Elongated cell")
        grid = QGridLayout(self.cell_box)
        self.uc_in_x = self._spin(1, 200, 22, 1)
        self.uc_in_y = self._spin(1, 100, 6, 1)
        self.free_space = self._spin(0.0, 0.45, 0.24, 0.01, dec=3)
        self.temp_transition = self._spin(1, 5000, 900, 50, dec=1, suffix=" K")
        self.lambdam = self._spin(0.0, 1.0, 0.0, 0.05, dec=2)
        # two columns of parameters keep the strip low
        for pos, (label, w, tip) in enumerate((
                ("unit cells in x", self.uc_in_x, "uc_in_x"),
                ("unit cells in y", self.uc_in_y, "uc_in_y"),
                ("free space", self.free_space,
                 "free_space: fraction of the box left empty on each side"),
                ("T transition zone", self.temp_transition,
                 "temperature_in_transition_zone"),
                ("λ_m (gas)", self.lambdam,
                 "lambdam: sqrt of the damping in the gas zone; 0 = ideal gas"))):
            lab = QLabel(label)
            lab.setToolTip(tip)
            w.setToolTip(tip)
            grid.addWidget(lab, pos % 3, 2 * (pos // 3))
            grid.addWidget(w, pos % 3, 2 * (pos // 3) + 1)
        row.addWidget(self.cell_box)

        self.preview_box = QGroupBox("Field preview")
        grid2 = QGridLayout(self.preview_box)
        self.temp = self._spin(1, 3000, 300, 10, dec=1, suffix=" K")
        self.um = self._spin(-100000, 100000, -5000, 500, dec=0, suffix=" J/mol")
        grid2.addWidget(QLabel("T"), 0, 0)
        grid2.addWidget(self.temp, 0, 1)
        grid2.addWidget(QLabel("u_m"), 1, 0)
        grid2.addWidget(self.um, 1, 1)
        note = QLabel("The T and u_m loops of a run are set on the Run "
                      "tab; these values only drive the preview curves.")
        note.setWordWrap(True)
        note.setProperty("dim", True)
        note.setMaximumWidth(240)
        grid2.addWidget(note, 2, 0, 1, 2)
        row.addWidget(self.preview_box)
        return panel

    # -- state ----------------------------------------------------------------

    def _spec(self) -> SimCellSpec:
        return SimCellSpec(int(self.uc_in_x.value()), int(self.uc_in_y.value()),
                           self.free_space.value(), self.lambdam.value(),
                           self.temp_transition.value())

    def refresh(self) -> None:
        """Called when the tab becomes current: pick up the unit cell, the
        glyph and (unless edited here) the saved settings."""
        saved = self.project.simulation_cell
        if saved is not None and not self._dirty:
            for w, key in ((self.uc_in_x, "uc_in_x"), (self.uc_in_y, "uc_in_y"),
                           (self.free_space, "free_space"),
                           (self.lambdam, "lambdam"),
                           (self.temp_transition, "temperature_in_transition_zone"),
                           (self.temp, "preview_temp"), (self.um, "preview_um")):
                if key in saved:
                    w.blockSignals(True)
                    w.setValue(saved[key])
                    w.blockSignals(False)
        self._refresh_slot()
        self._recompute()

    def _on_changed(self) -> None:
        self._dirty = True
        self._recompute()

    def _recompute(self) -> None:
        uc = self.project.unit_cell
        if uc is None or not uc["molecules"]:
            self.canvas.clear_data()
            self.cell_box.setEnabled(False)
            self.preview_box.setEnabled(False)
            self.info.setText("No unit cell yet. Build and save one on the "
                              "Unit cell tab first.")
            return
        self.cell_box.setEnabled(True)
        self.preview_box.setEnabled(True)
        layout, placements = tile_unit_cell(uc, self._spec())
        self.canvas.set_data(layout, placements, self.temp.value(),
                             self.um.value(), model_glyph(self.project))
        n = layout.n_molecules
        crystal = (len(uc["molecules"]) / (uc["cell_x"] * uc["cell_y"])
                   * 1.0e26 / N_A)
        overall = n / (layout.lx * layout.ly) * 1.0e26 / N_A
        warn = ("" if n <= ENGINE_BUFFER else
                f"  <span style='color:{theme.DANGER}'>&gt; engine buffer "
                f"{ENGINE_BUFFER}</span>")
        self.info.setText(
            f"{n} molecules{warn}<br>"
            f"box {layout.lx:.1f} × {layout.ly:.1f} Å, crystal slab "
            f"{layout.slab_x:.1f} Å<br>"
            f"density {overall:.4f} µmol/m² (crystal {crystal:.4f})")

    # -- project integration ---------------------------------------------------

    def _refresh_slot(self) -> None:
        saved = self.project.simulation_cell
        if saved is None:
            self.slot_label.setText("Not saved yet. Tune the cell and press "
                                    "'Use in project'.")
        else:
            self.slot_label.setText(
                f"{saved['uc_in_x']} × {saved['uc_in_y']} unit cells, "
                f"free space {saved['free_space']:g}, "
                f"T_tz {saved['temperature_in_transition_zone']:g} K, "
                f"λ_m {saved['lambdam']:g}")

    def use_in_project(self) -> None:
        if self.project.unit_cell is None:
            QMessageBox.information(self, "No unit cell",
                                    "Save a unit cell on the Unit cell tab first.")
            return
        spec = self._spec()
        self.project.set_simulation_cell({
            "uc_in_x": spec.uc_in_x, "uc_in_y": spec.uc_in_y,
            "free_space": spec.free_space, "lambdam": spec.lambdam,
            "temperature_in_transition_zone": spec.temp_transition,
            "preview_temp": self.temp.value(), "preview_um": self.um.value()})
        self._dirty = False
        self._refresh_slot()
        self.projectCellChanged.emit()
        self.statusMessage.emit("Simulation cell saved to the project")
