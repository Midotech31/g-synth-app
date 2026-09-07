import type { DesignParams } from "../api/client";

export default function InsertSettingsSummary({ params, fragment = true }: { params: DesignParams; fragment?: boolean }) {
  return <div className="notice notice-info compact">
    Current insert settings: {params.include_his_tag ? "6×His tag enabled" : "no added His tag"};{" "}
    {params.include_linkers ? "flexible linkers enabled" : "no added linkers"};{" "}
    {params.cleavage_site ? `${params.cleavage_site} site` : "no protease site"}.{" "}
    {fragment ? `${params.target_oligo_length ?? 90} nt target oligos with ${params.overhang_length ?? 4} nt assembly junctions.` : "Use the SSD duplex without fragmentation."}
  </div>;
}
