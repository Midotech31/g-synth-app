import { useMemo, useState } from "react";

import type { GelSimulation as GelSimulationData } from "../api/client";

export function gelBandPosition(size: number, minimum: number, maximum: number): number {
  if (maximum <= minimum || size <= 0) return 50;
  const clamped = Math.min(maximum, Math.max(minimum, size));
  const fraction = (Math.log10(maximum) - Math.log10(clamped))
    / (Math.log10(maximum) - Math.log10(minimum));
  return 6 + fraction * 88;
}

export default function GelSimulation({ simulation }: { simulation: GelSimulationData }) {
  const [ladderKey, setLadderKey] = useState(simulation.recommended_ladder);
  const ladder = simulation.ladders.find((item) => item.key === ladderKey)
    ?? simulation.ladders[0];
  const plottedLanes = useMemo(() => [
    {
      name: "M",
      description: ladder?.name ?? "DNA ladder",
      bands: (ladder?.bands ?? []).map((size) => ({ size_bp: size, label: `${size.toLocaleString()} bp` })),
      marker: true,
    },
    ...simulation.lanes.map((lane) => ({ ...lane, marker: false })),
  ], [ladder, simulation.lanes]);
  const sizes = plottedLanes.flatMap((lane) => lane.bands.map((band) => band.size_bp));
  const minimum = Math.max(20, Math.min(...sizes, 100) * 0.8);
  const maximum = Math.max(...sizes, 1500) * 1.15;

  return (
    <div className="gel-simulation">
      <div className="gel-toolbar">
        <div>
          <span className="pill pill-warn">In silico prediction</span>
          <p>{simulation.notice}</p>
        </div>
        <div className="field gel-ladder-picker">
          <label htmlFor={`gel-ladder-${simulation.title.replace(/\W+/g, "-")}`}>DNA size marker</label>
          <select
            id={`gel-ladder-${simulation.title.replace(/\W+/g, "-")}`}
            value={ladder?.key}
            onChange={(event) => setLadderKey(event.target.value)}
          >
            {simulation.ladders.map((item) => (
              <option key={item.key} value={item.key}>{item.name} · {item.range}</option>
            ))}
          </select>
        </div>
      </div>

      {simulation.lanes.length ? (
        <div className="gel-layout">
          <div className="virtual-gel" role="img" aria-label={`${simulation.title}. Predicted DNA bands by size.`}>
            {plottedLanes.map((lane, laneIndex) => (
              <div className="gel-lane" key={`${lane.name}-${laneIndex}`}>
                <div className="gel-well" />
                {lane.bands.map((band, bandIndex) => (
                  <span
                    className={`gel-band ${lane.marker ? "marker" : "sample"}`}
                    key={`${band.size_bp}-${bandIndex}`}
                    style={{ top: `${gelBandPosition(band.size_bp, minimum, maximum)}%` }}
                    title={`${lane.name}: ${band.label}`}
                  >
                    {lane.marker && <i>{band.size_bp.toLocaleString()}</i>}
                  </span>
                ))}
                <strong>{lane.name}</strong>
              </div>
            ))}
          </div>
          <div className="gel-lane-key">
            {plottedLanes.map((lane, index) => (
              <div key={`${lane.name}-key-${index}`}>
                <strong>{lane.name}</strong>
                <span>{lane.description}</span>
                <small>
                  {lane.bands.length
                    ? lane.bands.map((band) => band.label).join(" · ")
                    : "No specific band expected"}
                </small>
              </div>
            ))}
          </div>
        </div>
      ) : (
        <div className="notice notice-info">{simulation.notice}</div>
      )}

      {ladder && (
        <p className="gel-marker-list">
          <strong>{ladder.name} markers:</strong>{" "}
          {ladder.bands.map((size) => size.toLocaleString()).join(", ")} bp.
        </p>
      )}
    </div>
  );
}
