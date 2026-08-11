import type { Catalogue, DesignParams } from "../api/client";

/**
 * The insert controls, shared by the design and cloning pages.
 *
 * Both pages design the same cassette from the same inputs; keeping one form
 * means they cannot drift into offering different options for what is
 * supposed to be the same operation.
 */

type Props = {
  params: DesignParams;
  catalogue: Catalogue | null;
  onChange: <K extends keyof DesignParams>(key: K, value: DesignParams[K]) => void;
  /** Hide the fragmentation controls where they do not apply. */
  showFragmentation?: boolean;
  idPrefix?: string;
};

export default function InsertForm({
  params,
  catalogue,
  onChange,
  showFragmentation = true,
  idPrefix = "",
}: Props) {
  const id = (name: string) => `${idPrefix}${name}`;

  const enzymes = catalogue?.enzymes ?? [];
  const byName = new Map(enzymes.map((e) => [e.name, e]));
  // A hundred names in one flat list is worse than nineteen. The set this
  // lab keeps in the freezer goes first; the rest stay selectable, because
  // an enzyme in your vector's polylinker should not need a code change.
  const common = enzymes.filter((e) => e.common !== false);
  const rest = enzymes.filter((e) => e.common === false);
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
        <div className="field">
          <label htmlFor={id("left")}>5' enzyme</label>
          <select
            id={id("left")}
            value={params.left_enzyme}
            onChange={(e) => onChange("left_enzyme", e.target.value)}
          >
            <optgroup label="In the freezer">
              {common.map((e) => (
                <option key={e.name} value={e.name}>
                  {e.name} · {e.recognition}
                </option>
              ))}
            </optgroup>
            {rest.length > 0 && (
              <optgroup label={`Others (${rest.length})`}>
                {rest.map((e) => (
                  <option key={e.name} value={e.name}>
                    {e.name} · {e.recognition}
                  </option>
                ))}
              </optgroup>
            )}
          </select>
          {left && (
            <span className="label">
              {left.overhang ? `${left.overhang_type} ${left.overhang}` : "blunt"}
              {left.supplies_start_codon && " · supplies ATG"}
            </span>
          )}
        </div>
        <div className="field">
          <label htmlFor={id("right")}>3' enzyme</label>
          <select
            id={id("right")}
            value={params.right_enzyme}
            onChange={(e) => onChange("right_enzyme", e.target.value)}
          >
            <optgroup label="In the freezer">
              {common.map((e) => (
                <option key={e.name} value={e.name}>
                  {e.name} · {e.recognition}
                </option>
              ))}
            </optgroup>
            {rest.length > 0 && (
              <optgroup label={`Others (${rest.length})`}>
                {rest.map((e) => (
                  <option key={e.name} value={e.name}>
                    {e.name} · {e.recognition}
                  </option>
                ))}
              </optgroup>
            )}
          </select>
          {right && (
            <span className="label">
              {right.overhang ? `${right.overhang_type} ${right.overhang}` : "blunt"}
            </span>
          )}
        </div>
      </div>

      <div className="field">
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
      </div>

      <div className="checks">
        <label>
          <input
            type="checkbox"
            checked={params.include_his_tag}
            onChange={(e) => onChange("include_his_tag", e.target.checked)}
          />
          6×His tag
        </label>
        <label>
          <input
            type="checkbox"
            checked={params.include_linkers}
            onChange={(e) => onChange("include_linkers", e.target.checked)}
          />
          Flexible linkers
        </label>
        <label>
          <input
            type="checkbox"
            checked={params.is_coding}
            onChange={(e) => onChange("is_coding", e.target.checked)}
          />
          Insert already has its own ATG
        </label>
        {params.is_coding && (
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

      {showFragmentation && (
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
