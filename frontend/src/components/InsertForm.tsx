import type { Catalogue, DesignParams } from "../api/client";
import EnzymePicker from "./EnzymePicker";


type Props = {
  params: DesignParams;
  catalogue: Catalogue | null;
  onChange: <K extends keyof DesignParams>(key: K, value: DesignParams[K]) => void;

  showFragmentation?: boolean;

  expert?: boolean;
  idPrefix?: string;
};

export default function InsertForm({
  params,
  catalogue,
  onChange,
  showFragmentation = true,
  expert = true,
  idPrefix = "",
}: Props) {
  const id = (name: string) => `${idPrefix}${name}`;

  const enzymes = catalogue?.enzymes ?? [];
  const byName = new Map(enzymes.map((e) => [e.name, e]));
  const left = byName.get(params.left_enzyme);
  const right = byName.get(params.right_enzyme);

  return (
    <div style={{ display: "flex", flexDirection: "column", gap: "0.85rem" }}>
      <div className="field">
        <label htmlFor={id("name")}>Construct name</label>
        <input
          id={id("name")}
          type="text"
          value={params.name ?? ""}
          onChange={(e) => onChange("name", e.target.value)}
          placeholder="pGS-EntA"
        />
      </div>

      <div className="field">
        <label htmlFor={id("sequence")}>Sequence (A/C/G/T)</label>
        <textarea
          id={id("sequence")}
          value={params.sequence}
          onChange={(e) => onChange("sequence", e.target.value)}
          rows={6}
          className="mono"
          style={{ fontSize: "0.8rem" }}
        />
        <span className="label">
          {params.sequence.replace(/[^ACGTacgt]/g, "").length} nt entered
        </span>
      </div>

      <div className="row-2">
        <div>
          <EnzymePicker
            id={id("left")}
            label="5' enzyme"
            enzymes={enzymes}
            value={params.left_enzyme}
            onChange={(value) => onChange("left_enzyme", value)}
          />
          {left && (
            <span className="label">
              {left.overhang ? `${left.overhang_type} ${left.overhang}` : "blunt"}
              {left.supplies_start_codon && " · supplies ATG"}
            </span>
          )}
        </div>
        <div>
          <EnzymePicker
            id={id("right")}
            label="3' enzyme"
            enzymes={enzymes}
            value={params.right_enzyme}
            onChange={(value) => onChange("right_enzyme", value)}
          />
          {right && (
            <span className="label">
              {right.overhang ? `${right.overhang_type} ${right.overhang}` : "blunt"}
            </span>
          )}
        </div>
      </div>

      {expert && <div className="field">
        <label htmlFor={id("cleavage")}>Protease site</label>
        <select
          id={id("cleavage")}
          value={params.cleavage_site ?? ""}
          onChange={(e) => onChange("cleavage_site", e.target.value || null)}
        >
          <option value="">None</option>
          {catalogue?.cleavage_sites.map((c) => (
            <option key={c.name} value={c.name}>
              {c.name}
            </option>
          ))}
        </select>
      </div>}

      <div className="checks">
        {expert && <label>
          <input
            type="checkbox"
            checked={params.include_his_tag}
            onChange={(e) => onChange("include_his_tag", e.target.checked)}
          />
          6×His tag
        </label>}
        {expert && <label>
          <input
            type="checkbox"
            checked={params.include_linkers}
            onChange={(e) => onChange("include_linkers", e.target.checked)}
          />
          Flexible linkers
        </label>}
        <label>
          <input
            type="checkbox"
            checked={params.is_coding}
            onChange={(e) => onChange("is_coding", e.target.checked)}
          />
          Standalone ORF (starts with ATG)
        </label>
        {expert && params.is_coding && (
          <label>
            <input
              type="checkbox"
              checked={params.remove_stop}
              onChange={(e) => onChange("remove_stop", e.target.checked)}
            />
            Remove the stop codon
          </label>
        )}
      </div>

      {expert && showFragmentation && (
        <div className="row-2">
          <div className="field">
            <label htmlFor={id("oligo")}>Oligo length (nt)</label>
            <input
              id={id("oligo")}
              type="number"
              min={20}
              max={300}
              value={params.target_oligo_length}
              onChange={(e) => onChange("target_oligo_length", Number(e.target.value))}
            />
          </div>
          <div className="field">
            <label htmlFor={id("overhang")}>Junction overhang (nt)</label>
            <input
              id={id("overhang")}
              type="number"
              min={4}
              max={8}
              value={params.overhang_length}
              onChange={(e) => onChange("overhang_length", Number(e.target.value))}
            />
          </div>
        </div>
      )}
    </div>
  );
}
