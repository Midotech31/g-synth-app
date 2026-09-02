"""Render the G-Synth double-helix and leaf brand mark."""
import math
import pathlib

CX, AMP = 40.0, 21.0
TOP, BOT = 14.0, 126.0
TURNS = 1.0
W_MAX, W_MIN = 9.0, 2.4
STEPS = 110

VIEWBOX = (12, 9, 80, 121)

NAVY, MID, TEAL = "#0b2545", "#1b5c86", "#34a0bd"


def strand(phase: float) -> str:
    """One ribbon, as a closed outline whose width breathes along the sweep.

    A helix seen side-on shows its ribbon broadside at the extremes and
    edge-on where the strands cross, so the width follows |sin|.
    """
    left, right = [], []
    for i in range(STEPS + 1):
        f = i / STEPS
        t = phase + f * TURNS * 2 * math.pi
        y = TOP + f * (BOT - TOP)
        x = CX + AMP * math.sin(t)

        half = (W_MIN + (W_MAX - W_MIN) * abs(math.sin(t))) / 2
        half *= 0.42 + 0.58 * min(1.0, f / 0.05, (1 - f) / 0.05)

        # Offset along the normal to preserve ribbon width through curves.
        dx = AMP * math.cos(t) * TURNS * 2 * math.pi
        dy = BOT - TOP
        n = math.hypot(dx, dy)

        left.append((x + dy / n * half, y - dx / n * half))
        right.append((x - dy / n * half, y + dx / n * half))

    pts = left + right[::-1]
    return (f"M{pts[0][0]:.1f} {pts[0][1]:.1f}"
            + "".join(f"L{x:.1f} {y:.1f}" for x, y in pts[1:]) + "Z")


def rungs() -> list[str]:
    """Return legible base-pair rungs away from helix crossings."""
    out = []
    for i in range(1, 10):
        f = i / 10
        t = f * TURNS * 2 * math.pi
        sep = abs(math.sin(t))
        if sep < 0.50:
            continue
        y = TOP + f * (BOT - TOP)
        span = AMP * sep
        inset = 3.4
        out.append(
            f'<line x1="{CX - span + inset:.2f}" y1="{y:.2f}" '
            f'x2="{CX + span - inset:.2f}" y2="{y:.2f}" '
            f'stroke-width="{1.5 + 1.6 * sep:.2f}"/>'
        )
    return out


def leaf(base, tip, bulge, back=0.82):
    """An almond leaf: two arcs bowed to opposite sides of the base–tip axis."""
    (bx, by), (tx, ty) = base, tip
    mx, my = (tx + bx) / 2, (ty + by) / 2
    dx, dy = tx - bx, ty - by
    n = math.hypot(dx, dy)
    px, py = -dy / n * bulge, dx / n * bulge
    return (f"M{bx:.1f} {by:.1f} "
            f"Q{mx + px:.1f} {my + py:.1f} {tx:.1f} {ty:.1f} "
            f"Q{mx - px * back:.1f} {my - py * back:.1f} {bx:.1f} {by:.1f} Z")


def midrib(base, tip, bulge):
    """The vein, bowed a little less than the leaf's own upper edge."""
    (bx, by), (tx, ty) = base, tip
    mx, my = (tx + bx) / 2, (ty + by) / 2
    dx, dy = tx - bx, ty - by
    n = math.hypot(dx, dy)
    px, py = -dy / n * bulge * 0.28, dx / n * bulge * 0.28
    return f"M{bx:.1f} {by:.1f} Q{mx + px:.1f} {my + py:.1f} {tx:.1f} {ty:.1f}"


INNER = ((49.5, 41.0), (57.0, 12.0), 13.0)
OUTER = ((51.0, 43.0), (89.0, 19.0), 17.5)

svg = f'''<svg xmlns="http://www.w3.org/2000/svg" viewBox="12 9 80 121" fill="none" role="img">
  <defs>
    <linearGradient id="gsHelix" x1="16" y1="126" x2="76" y2="14" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{NAVY}"/>
      <stop offset="0.62" stop-color="{MID}"/>
      <stop offset="1" stop-color="{TEAL}"/>
    </linearGradient>
    <linearGradient id="gsLeaf" x1="54" y1="44" x2="88" y2="18" gradientUnits="userSpaceOnUse">
      <stop offset="0" stop-color="{MID}"/>
      <stop offset="1" stop-color="{TEAL}"/>
    </linearGradient>
  </defs>

  <g fill="url(#gsHelix)">
    <path d="{strand(0.0)}"/>
    <path d="{strand(math.pi)}"/>
  </g>

  <g stroke="url(#gsHelix)" stroke-linecap="round" opacity="0.5">
    {chr(10).join("    " + r for r in rungs()).strip()}
  </g>

  <path d="{leaf(*INNER)}" fill="{NAVY}"/>
  <path d="{leaf(*OUTER)}" fill="url(#gsLeaf)"/>
  <path d="{midrib(*OUTER)}" stroke="#ffffff" stroke-width="1.1" opacity="0.26"/>
</svg>
'''

out = pathlib.Path(__file__).with_name("mark.svg")
out.write_text(svg)
print("wrote", out, len(svg), "bytes")
