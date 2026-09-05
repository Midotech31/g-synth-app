import { useMemo, useState } from "react";

import type { Enzyme } from "../api/client";

type Props = {
  id: string;
  label: string;
  enzymes: Enzyme[];
  value: string;
  onChange: (value: string) => void;
  disabled?: boolean;
};

function searchableText(enzyme: Enzyme): string {
  return [
    enzyme.name,
    ...(enzyme.aliases ?? []),
    enzyme.recognition,
    enzyme.overhang,
    enzyme.overhang_type,
  ].join(" ").toLocaleLowerCase();
}

function optionLabel(enzyme: Enzyme): string {
  const geometry = enzyme.overhang
    ? `${enzyme.overhang_type} ${enzyme.overhang}`
    : "blunt";
  return `${[enzyme.name, ...(enzyme.aliases ?? [])].join(" / ")} · ${enzyme.recognition} · ${geometry}`;
}

export default function EnzymePicker({ id, label, enzymes, value, onChange, disabled = false }: Props) {
  const [query, setQuery] = useState("");
  const filtered = useMemo(() => {
    const needle = query.trim().toLocaleLowerCase();
    if (!needle) return enzymes;
    const matches = enzymes.filter((enzyme) => searchableText(enzyme).includes(needle));
    const selected = enzymes.find((enzyme) => enzyme.name === value);
    if (selected && !matches.includes(selected)) return [selected, ...matches];
    return matches;
  }, [enzymes, query, value]);
  const common = filtered.filter((enzyme) => enzyme.common !== false);
  const additional = filtered.filter((enzyme) => enzyme.common === false);
  const actualMatches = filtered.filter((enzyme) => enzyme.name !== value || searchableText(enzyme).includes(query.trim().toLocaleLowerCase())).length;

  return (
    <div className="field enzyme-picker">
      <label htmlFor={id}>{label}</label>
      <input
        type="search"
        value={query}
        onChange={(event) => setQuery(event.target.value)}
        aria-label={`Search ${label.toLocaleLowerCase()}`}
        aria-controls={id}
        placeholder="Search name, site or overhang"
        autoComplete="off"
        disabled={disabled}
      />
      <select id={id} value={value} onChange={(event) => onChange(event.target.value)} disabled={disabled}>
        {enzymes.length === 0 && <option value={value}>{value}</option>}
        {common.length > 0 && (
          <optgroup label="Common cloning enzymes">
            {common.map((enzyme) => (
              <option key={enzyme.name} value={enzyme.name}>{optionLabel(enzyme)}</option>
            ))}
          </optgroup>
        )}
        {additional.length > 0 && (
          <optgroup label={`Additional verified enzymes (${additional.length})`}>
            {additional.map((enzyme) => (
              <option key={enzyme.name} value={enzyme.name}>{optionLabel(enzyme)}</option>
            ))}
          </optgroup>
        )}
      </select>
      {query && <span className="label" role="status">{actualMatches} matching enzyme{actualMatches === 1 ? "" : "s"}</span>}
    </div>
  );
}
